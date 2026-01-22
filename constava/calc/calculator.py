"""constava.calculator contains the calculator class that calculates the 
conformational state propensities and conformational state variability from
a protein ensemble
"""

import logging
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List
import tqdm
from .subsampling import SubsamplingABC, SubsamplingMethodError
from .csmodels import ConfStateModelABC
from ..utils.ensembles import ProteinEnsemble
from ..utils.results import ConstavaResults, ConstavaResultsEntry

# The logger for the wrapper
logger = logging.getLogger("Constava")

_LOGPDF_WORKER_CSMODEL = None

def _init_logpdf_worker(csmodel: ConfStateModelABC):
    """Initializer to share the conformational state model with worker processes."""
    global _LOGPDF_WORKER_CSMODEL
    _LOGPDF_WORKER_CSMODEL = csmodel

def _compute_logpdf_worker(phipsi):
    """Worker function executed in separate processes to compute logpdf."""
    return _LOGPDF_WORKER_CSMODEL.get_logpdf(phipsi)

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
        logger.debug(f"Calculating the state propensities and variability for each residue ({len(residues)}) for each method ({len(self.methods)} methods)...")

        logpdf_workers = min(len(residues), 8) or 1
        logger.debug(f"Calculating log-probability densities of ({len(residues)}) residues using {logpdf_workers} process workers ...")

        logpdfs = [None] * len(residues)
        
        with tqdm.tqdm(total=len(residues), unit="residue") as pbar:
            with ProcessPoolExecutor(
                max_workers=logpdf_workers,
                initializer=_init_logpdf_worker,
                initargs=(self.csmodels,)
            ) as executor:
                futures = {
                    executor.submit(_compute_logpdf_worker, residue.phipsi): idx
                    for idx, residue in enumerate(residues)
                }
                for future in as_completed(futures):
                    idx = futures[future]
                    logpdfs[idx] = future.result()
                    pbar.update(1)

        logger.debug("Finished computing log-probability densities, aggregating per-method results...")

        for res, logpdf in tqdm.tqdm(zip(residues, logpdfs), total=len(residues), unit='residues'):
            for method, result in zip(self.methods, results):
                state_propensities, state_variability = method.calculate(logpdf)

                result.add_entry(ConstavaResultsEntry(
                    res, state_propensities, state_variability))
        
        logger.debug(f"Returning the state propensities and variability results ({len(results)} results)...")  
        
        return results
