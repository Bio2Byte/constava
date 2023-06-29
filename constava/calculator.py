from typing import List
import tqdm
from .ensembles import ProteinEnsemble
from .methods import ConstavaABC
from .pdfestimators import StatePdfABC
from .results import ConfStateResults, ConfStateResultsEntry


class ConfStateCalculator:

    def __init__(self, pdfestimator: StatePdfABC, methods: List[ConstavaABC] = None):
        self.pdfestimator = pdfestimator
        self.methods = methods or []

    def add_method(self, new_method: ConstavaABC):
        self.methods.append(new_method)

    def calculate(self, ensemble: ProteinEnsemble) -> List[ConfStateResults]:
        
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

