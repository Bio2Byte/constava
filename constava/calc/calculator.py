"""constava.calculator contains the calculator class that calculates the
conformational state propensities and conformational state variability from
a protein ensemble
"""

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing
from typing import List
import tqdm
from .subsampling import SubsamplingABC, SubsamplingMethodError
from .csmodels import ConfStateModelABC
from ..utils.ensembles import ProteinEnsemble
from ..utils.results import ConstavaResults, ConstavaResultsEntry

# The logger for the wrapper
logger = logging.getLogger("Constava")

_LOGPDF_WORKER_CSMODEL = None

def _init_logpdf_worker(csmodel: ConfStateModelABC, methods: List[SubsamplingABC]):
    """Initializer to share the conformational state model with worker processes."""
    global _LOGPDF_WORKER_CSMODEL
    global _LOGPDF_WORKER_METHODS

    _LOGPDF_WORKER_CSMODEL = csmodel
    _LOGPDF_WORKER_METHODS = methods

def _compute_logpdf_worker(phipsi):
    """Worker function executed in separate processes to compute logpdf."""
    logpdf = _LOGPDF_WORKER_CSMODEL.get_logpdf(phipsi)
    
    result = { "logpdf": logpdf }
    
    for method in _LOGPDF_WORKER_METHODS:
        state_propensities, state_variability = method.calculate(logpdf)

        result[method.getShortName()] = { 
            "state_propensities": state_propensities, 
            "state_variability": state_variability 
        }

    return result


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

        residues = list(ensemble.get_residues())
        n_residues = ensemble.n_residues
        logpdf_workers = min(4, multiprocessing.cpu_count())
        logpdf_results = {}

        with tqdm.tqdm(total=n_residues, unit="residue") as pbar:
            with ProcessPoolExecutor(
                max_workers=logpdf_workers,
                initializer=_init_logpdf_worker,
                initargs=(self.csmodels, self.methods)
            ) as executor:
                futures = {
                    executor.submit(_compute_logpdf_worker, residue.phipsi): idx
                    for idx, residue in enumerate(residues)
                }

                for future in as_completed(futures):
                    idx = futures[future]
                    logpdf_results[idx] = future.result()
                    pbar.update(1)

        for idx, res in enumerate(residues):
            res_logpdf_results = logpdf_results[idx]

            for result in results:
                current_method = result.method
                
                result.add_entry(
                    ConstavaResultsEntry(
                        res,
                        res_logpdf_results[current_method]["state_propensities"], 
                        res_logpdf_results[current_method]["state_variability"]
                    )
                )

        return results
