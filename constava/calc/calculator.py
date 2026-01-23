"""constava.calculator contains the calculator class that calculates the
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

def _self_compute_logpdf_worker(phipsi: np.ndarray, csmodel: ConfStateModelABC, methods: List[SubsamplingABC]):
    """Worker function executed in separate processes to compute logpdf."""
    logpdf = csmodel.get_logpdf(phipsi)
    
    result = dict()
    
    for method in methods:
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

        logpdf_results = [None] * n_residues
        # import ipdb; ipdb.set_trace()
        
        try:
            max_workers = os.process_cpu_count()
        except AttributeError:
            max_workers = multiprocessing.cpu_count()

        max_in_flight = max_workers * 4  # tweak: 1–4x workers is usually sane
        
        logger.debug(f"Starting the parallel inference of log-probability densities and calculations for propensities & variability...")
        with ProcessPoolExecutor(max_workers=max_workers) as ex:    
            it = iter(enumerate(residues))
            in_flight = set()
            
            logger.debug(f"Process Pool Executor is ready to start: max_in_flight={max_in_flight}, max_workers={max_workers}")
            
            with tqdm.tqdm(total=n_residues, unit="residue") as progress_bar:                
                
                # Prime the pipeline
                for _ in range(min(max_in_flight, n_residues)):
                    idx, residue = next(it)
                    phipsi = np.ascontiguousarray(residue.phipsi, dtype=np.float32)
                    
                    fut = ex.submit(_self_compute_logpdf_worker, phipsi, self.csmodels, self.methods)
                    
                    fut.idx = idx  # attach index (simple and cheap)
                    in_flight.add(fut)

                completed = 0
                while in_flight:
                    done, in_flight = wait(in_flight, return_when=FIRST_COMPLETED)
                    
                    for fut in done:
                        idx = fut.idx
                        logpdf_results[idx] = fut.result()
                        
                        progress_bar.update(1)
                        completed += 1

                        # refill
                        try:
                            idx, residue = next(it)
                        except StopIteration:
                            continue
                        
                        phipsi = np.ascontiguousarray(residue.phipsi, dtype=np.float32)
                        
                        nfut = ex.submit(_self_compute_logpdf_worker, phipsi, self.csmodels, self.methods)
                        
                        nfut.idx = idx
                        in_flight.add(nfut)

        logger.debug(f"Building the results objects...")
        for idx, residue in enumerate(residues):
            res_logpdf_results = logpdf_results[idx]

            for result in results:
                current_method = result.method
                
                result.add_entry(
                    ConstavaResultsEntry(
                        residue,
                        res_logpdf_results[current_method]["state_propensities"], 
                        res_logpdf_results[current_method]["state_variability"]
                    )
                )

        return results
