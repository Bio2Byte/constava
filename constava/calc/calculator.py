"""constava.calculator contains the calculator class that calculates the 
conformational state propensities and conformational state variability from
a protein ensemble
"""

from typing import List
import tqdm
from .subsampling import SubsamplingABC
from .pdfestimators import StatePdfABC
from ..datautils.ensembles import ProteinEnsemble
from ..datautils.results import ConfStateResults, ConfStateResultsEntry


class ConfStateCalculator:

    def __init__(self, pdfestimator: StatePdfABC, methods: List[SubsamplingABC] = None):
        """Initializes the calcualtor class with given conformational state 
        models (pdfestimator) and zero or more subsampling methods.

        Parameters:
        -----------
            pdfestimator : StatePdfABC
                Probabilistic model of conformational states used for the calculation.

            methods : List[SubsamplingABC] = None
                (Optional) A list of subsampling methods to use in the calculation.
        """
        self.pdfestimator = pdfestimator
        self.methods = methods or []

    def add_method(self, new_method: SubsamplingABC):
        """Adds a new subsampling methods to the calculator."""
        self.methods.append(new_method)

    def calculate(self, ensemble: ProteinEnsemble) -> List[ConfStateResults]:
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
            results : List[ConfStateResults]
                A list of ConfStateResults, each representing the same 
                ProteinEnsemble but calculated with one or more subsampling methods.
        """
        results = [
            ConfStateResults(method, ensemble, self.pdfestimator.labels) 
            for method in self.methods
        ]

        for res in tqdm.tqdm(ensemble.get_residues(), total=ensemble.n_residues, unit='residues'):
            logpdf = self.pdfestimator.get_logpdf(res.phipsi)

            for meth, result in zip(self.methods, results):
                state_propensities, state_variability = meth.calculate(logpdf)
                result.add_entry(ConfStateResultsEntry(
                    res, state_propensities, state_variability))
                
        return results

