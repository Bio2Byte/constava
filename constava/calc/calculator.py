"""constava.calc.calculator contains the calculator class that calculates the
conformational state propensities and conformational state variability from
a protein ensemble
"""

import os
import multiprocessing
import logging
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from typing import List
import numpy as np

import tqdm
from .subsampling import SubsamplingABC, SubsamplingMethodError
from .csmodels import ConfStateModelABC
from ..utils.ensembles import ProteinEnsemble
from ..utils.results import ConstavaResults, ConstavaResultsEntry

# The logger for the wrapper
logger = logging.getLogger("Constava")

_SELF_WORKER_CSMODEL: ConfStateModelABC = None
_SELF_WORKER_METHODS: List[SubsamplingABC] = None


def _init_self_worker(csmodels, methods):
    global _SELF_WORKER_CSMODEL, _SELF_WORKER_METHODS
    _SELF_WORKER_CSMODEL = csmodels
    _SELF_WORKER_METHODS = methods


def _self_compute_logpdf_worker(phipsi: np.ndarray):
    """
    Worker function executed in a child process.

    Performance change:
    - Returns two lists aligned to method order instead of a nested defaultdict,
      greatly reducing pickle/serialization overhead.
    """
    csmodel = _SELF_WORKER_CSMODEL
    methods = _SELF_WORKER_METHODS

    logpdf = csmodel.get_logpdf(phipsi)

    # Return plain lists to minimize serialization costs.
    # props_list[i] and var_list[i] correspond to methods[i].
    props_list = [None] * len(methods)
    var_list = [None] * len(methods)

    for method_idx, method in enumerate(methods):
        state_propensities, state_variability = method.calculate(logpdf)

        props_list[method_idx] = state_propensities
        var_list[method_idx] = state_variability

    return props_list, var_list


def check_subsampling_methods(func):
    """Decorator for ConfStateCalculator that checks, if appropriate subsampling
    methods have been defined.
    """

    def __inner(self, *args, **kwargs):
        check = [isinstance(m, SubsamplingABC) for m in self.methods]
        if len(check) < 1:
            raise SubsamplingMethodError("No subsampling methods specified.")
        elif not all(check):
            raise SubsamplingMethodError(
                "Incompatible subsampling method found. All subsampling methods should inherit from SubsamplingABC."
            )
        return func(self, *args, **kwargs)

    return __inner


class ConfStateCalculator:

    def __init__(
        self, csmodels: ConfStateModelABC, methods: List[SubsamplingABC] = None
    ):
        """Initializes the calculator class with given conformational state
        models (csmodels) and zero or more sub-sampling methods.

        Parameters:
        -----------
            csmodels : ConfStateModelABC
                Probabilistic model of conformational states used for the calculation.

            methods : List[SubsamplingABC] = None
                (Optional) A list of sub-sampling methods to use in the calculation.
        """
        self.csmodels = csmodels
        self.methods = methods or []

        # Improvement: keep a reusable process pool to avoid paying startup cost
        # on every calculate() call.
        self._executor = None
        self._executor_max_workers = None
        self._mp_context = None

    def add_method(self, new_method: SubsamplingABC):
        """Adds a new sub-sampling methods to the calculator."""
        self.methods.append(new_method)

    def close(self):
        """Shut down the internal reusable pool, if any."""
        if self._executor is not None:
            try:
                self._executor.shutdown(wait=True, cancel_futures=False)
            finally:
                self._executor = None
                self._executor_max_workers = None
                self._mp_context = None

    def __del__(self):
        # Best-effort cleanup (don't raise during GC)
        try:
            self.close()
        except BaseException:
            pass

    def _get_mp_context(self):
        """
        Try to choose a fast start method.
        - 'fork' is typically fastest on Linux (low startup overhead).
        - On platforms where 'fork' isn't available/appropriate, fall back to default.
        """
        if self._mp_context is not None:
            return self._mp_context

        try:
            # "fork" is usually fastest on Linux.
            ctx = multiprocessing.get_context("fork")
        except (ValueError, RuntimeError, AttributeError):
            # Fallback to default context ("spawn" on Windows, often "spawn" on macOS).
            ctx = multiprocessing.get_context()

        self._mp_context = ctx
        return ctx

    def _get_or_create_executor(self, max_workers):
        """
        Reuse a cached executor if possible; otherwise create a new one.

        Note: Reusing requires that csmodels/methods are intended to stay the same.
        If you mutate methods frequently, consider calling close() and rebuilding.
        """
        if self._executor is not None and self._executor_max_workers == max_workers:
            return self._executor

        # If something changed, rebuild the pool.
        self.close()

        ctx = self._get_mp_context()

        self._executor = ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_init_self_worker,
            initargs=(self.csmodels, self.methods),
            mp_context=ctx,
        )
        self._executor_max_workers = max_workers
        return self._executor

    @check_subsampling_methods
    def calculate(self, ensemble: ProteinEnsemble) -> List[ConstavaResults]:
        """Calculates conformational state propensities and conformational
        state variabilities based on the given model and methods for the
        ProteinEnsemble.

        Parameters:
        -----------
            ensemble : ProteinEnsemble
                A ProteinEnsemble object with the relevant backbone
                dihedral information of the ensemble.

        Returns:
        --------
            results : List[ConstavaResults]
                A list of ConstavaResults, each representing the same
                ProteinEnsemble but calculated with one or more sub-sampling methods.
        """

        n_residues = ensemble.n_residues
        residues_propensities_variabilities = dict()
        future_to_idx = dict()

        try:
            max_workers = os.process_cpu_count()
        except AttributeError:
            max_workers = multiprocessing.cpu_count()

        max_in_flight = max_workers * 2

        logger.debug(
            "Starting inference of log-probability densities & propensities/variability"
        )

        # Improvement: reuse pool across calls to avoid process startup cost
        process_pool_executor = self._get_or_create_executor(max_workers=max_workers)

        it = iter(enumerate(ensemble.get_residues()))
        in_flight = set()

        with tqdm.tqdm(
            total=n_residues,
            desc="Residues",
            unit="residue",
            bar_format="{l_bar}{bar} | {n_fmt}/{total_fmt} "
            "[{rate_fmt}, elapsed: {elapsed}, remaining: {remaining}]",
        ) as progress_bar:
            # Prime the pipeline
            for _ in range(min(max_in_flight, n_residues)):
                idx, residue = next(it)

                fut = process_pool_executor.submit(
                    _self_compute_logpdf_worker, residue.phipsi
                )

                future_to_idx[fut] = idx
                in_flight.add(fut)

            completed = 0
            while in_flight:
                done, in_flight = wait(in_flight, return_when=FIRST_COMPLETED)

                for fut in done:
                    idx = future_to_idx.pop(fut)
                    residues_propensities_variabilities[idx] = fut.result()

                    progress_bar.update(1)
                    completed += 1

                    # refill
                    try:
                        idx, residue = next(it)
                    except StopIteration:
                        continue

                    nfut = process_pool_executor.submit(
                        _self_compute_logpdf_worker, residue.phipsi
                    )

                    future_to_idx[nfut] = idx
                    in_flight.add(nfut)

        results = [
            ConstavaResults(
                method=method.getShortName(),
                protein=ensemble,
                state_labels=self.csmodels.get_labels(),
            )
            for method in self.methods
        ]

        for idx, res in enumerate(ensemble.get_residues(sorted_list=True)):
            props_list, var_list = residues_propensities_variabilities[idx]

            for method_idx, result in enumerate(results):
                state_propensities = props_list[method_idx]
                state_variability = var_list[method_idx]

                result.add_entry(
                    ConstavaResultsEntry(res, state_propensities, state_variability),
                    sorted_insertion=True,
                )

        return results
