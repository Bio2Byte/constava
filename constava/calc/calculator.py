"""constava.calculator contains the calculator class that calculates the
conformational state propensities and conformational state variability from
a protein ensemble
"""

from collections import defaultdict
import os
import multiprocessing
import logging
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from typing import List
import numpy as np

import tqdm
from tqdm.contrib.logging import logging_redirect_tqdm
from .subsampling import SubsamplingABC, SubsamplingMethodError
from .csmodels import ConfStateModelABC
from ..utils.ensembles import ProteinEnsemble
from ..utils.results import ConstavaResults, ConstavaResultsEntry

# The logger for the wrapper
logger = logging.getLogger("Constava")

_SELF_WORKER_CSMODEL = None
_SELF_WORKER_METHODS = None

def _init_self_worker(csmodels, methods):
    global _SELF_WORKER_CSMODEL, _SELF_WORKER_METHODS
    _SELF_WORKER_CSMODEL = csmodels
    _SELF_WORKER_METHODS = methods

def _self_compute_logpdf_worker(phipsi: np.ndarray):
    csmodel = _SELF_WORKER_CSMODEL
    methods = _SELF_WORKER_METHODS

    logpdf = csmodel.get_logpdf(phipsi)

    n_methods = len(methods)
    n_states = len(csmodel.get_labels())

    matrix_results = np.empty((n_methods, n_states + 1))

    for method_idx, method in enumerate(methods):
        state_propensities, state_variability = method.calculate(logpdf)
        matrix_results[method_idx, :n_states] = state_propensities
        matrix_results[method_idx, n_states] = state_variability

    return matrix_results


def check_subsampling_methods(func):
    """Decorator for ConfStateCalculator that checks, if appropriate subsampling
    methods have been defined.
    """
    def __inner(self, *args, **kwargs):
        check = [isinstance(m, SubsamplingABC) for m in self.methods]
        if len(check) < 1:
            raise SubsamplingMethodError("No subsampling methods specified.")
        elif not all(check):
            raise SubsamplingMethodError("Incompatible subsampling method found. All subsampling methods should inherit from SubsamplingABC.")
        return func(self, *args, **kwargs)
    return __inner


class ConfStateCalculator:

    def __init__(self, csmodels: ConfStateModelABC, methods: List[SubsamplingABC] = None):
        """Initializes the calcualtor class with given conformational state
        models (csmodels) and zero or more subsampling methods.

        Parameters:
        -----------
            csmodels : ConfStateModelABC
                Probabilistic model of conformational states used for the calculation.

            methods : List[SubsamplingABC] = None
                (Optional) A list of subsampling methods to use in the calculation.
        """
        self.csmodels = csmodels
        self.methods = methods or []

    def add_method(self, new_method: SubsamplingABC):
        """Adds a new subsampling methods to the calculator."""
        self.methods.append(new_method)

    @check_subsampling_methods
    def calculate(self, ensemble: ProteinEnsemble) -> List[ConstavaResults]:
        """Calculates conformational state propensities and conformational
        state variabilites based on the given model and methods for the
        ProteinEnsemble.

        Parameters:
        -----------
            ensemble : ProteinEnsemble
                A ProteinEnsemble object with the relevant backbone
                dihedal information of the ensemble.

        Returns:
        --------
            results : List[ConstavaResults]
                A list of ConstavaResults, each representing the same
                ProteinEnsemble but calculated with one or more subsampling methods.
        """

        logger.debug(f"Instantiating Constava results for each method ({len(self.methods)} methods)...")
        results = [
            ConstavaResults(
                method=method.getShortName(),
                protein=ensemble,
                state_labels=self.csmodels.get_labels()
            )
            for method in self.methods
        ]

        n_residues = ensemble.n_residues

        logpdf_results = np.empty(
            (n_residues, len(self.methods), len(self.csmodels.get_labels()) + 1)
        )

        try:
            max_workers = os.process_cpu_count()
        except AttributeError:
            max_workers = multiprocessing.cpu_count()

        max_in_flight = max_workers * 3
        future_to_idx = defaultdict()
        
        logger.debug(f"Starting inference of log-probability densities & propens/var")
        
        with ProcessPoolExecutor(
            max_workers=max_workers,
            initializer=_init_self_worker,
            initargs=(self.csmodels, self.methods)
        ) as process_pool_executor:

            it = iter(enumerate(ensemble.get_residues()))
            in_flight = set()
            
            logger.debug(f"Max_in_flight={max_in_flight}, Max_workers={max_workers}")

            with tqdm.tqdm(
                total=n_residues, 
                desc="Residues", 
                unit="residue",
                bar_format="{l_bar}{bar} | {n_fmt}/{total_fmt} "
               "[{rate_fmt}, elapsed: {elapsed}, remaining: {remaining}]",
            ) as progress_bar:
                with logging_redirect_tqdm():
                    # Prime the pipeline
                    for _ in range(min(max_in_flight, n_residues)):
                        idx, residue = next(it)
                        phi_psi_angles = np.ascontiguousarray(residue.phipsi)

                        fut = process_pool_executor.submit(
                            _self_compute_logpdf_worker, phi_psi_angles
                        )
                            
                        future_to_idx[fut] = idx

                        in_flight.add(fut)

                    completed = 0
                    while in_flight:
                        done, in_flight = wait(in_flight, return_when=FIRST_COMPLETED)

                        for fut in done:
                            idx = future_to_idx.pop(fut)
                            logpdf_results[idx] = fut.result()

                            progress_bar.update(1)
                            completed += 1

                            # refill
                            try:
                                idx, residue = next(it)
                                phi_psi_angles = np.ascontiguousarray(residue.phipsi)
                            except StopIteration:
                                continue

                            nfut = process_pool_executor.submit(
                                _self_compute_logpdf_worker, phi_psi_angles
                            )
                            
                            # logger.debug(f"Parallel task submitted ({idx}: {residue})")
                            
                            future_to_idx[nfut] = idx
                            
                            in_flight.add(nfut)

        logger.debug("Building the results objects...")
        for idx, residue in enumerate(ensemble.get_residues()):
            res_logpdf_results = logpdf_results[idx]

            for result_idx, result in enumerate(results):
                state_propensities = res_logpdf_results[result_idx][:-1]
                state_variability = res_logpdf_results[result_idx][-1]

                result.add_entry(
                    ConstavaResultsEntry(
                        residue,
                        state_propensities,
                        state_variability
                    )
                )

        logger.debug("All the calculations have been done with success!")

        return results
