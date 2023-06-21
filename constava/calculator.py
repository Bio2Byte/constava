from dataclasses import dataclass
from typing import List
import numpy as np
import tqdm
from constava.ensembles import ProteinEnsemble, ResidueEnsemble
from constava.methods import ConstavaABC
from constava.pdfestimators import StatePdfABC


@dataclass
class ConfStateResultsEntry:
    residue: ResidueEnsemble = None
    state_propensities: np.ndarray = None
    state_variability: float = None


class ConfStateResults:
    """ A result class for ConfStateCalculator """
    def __init__(self, method: ConstavaABC, state_labels: List[str], entries: List[ConfStateResultsEntry] = None):
        self.method = method
        self.state_labels = state_labels
        self.entries = entries or []

    def add_entry(self, new_entry: ConfStateResultsEntry):
        self.entries.append(new_entry)
        self.entries.sort(key=lambda x: x.residue.respos)


class ConfStateCalculator:

    def __init__(self, pdfestimator: StatePdfABC, methods: List[ConstavaABC] = None):
        self.pdfestimator = pdfestimator
        self.methods = methods or []

    def add_method(self, new_method: ConstavaABC):
        self.methods.append(new_method)

    def calculate(self, ensemble: ProteinEnsemble) -> List[ConfStateResults]:
        
        results = [
            ConfStateResults(method, self.pdfestimator.labels) 
            for method in self.methods
        ]

        for res in tqdm.tqdm(ensemble.get_residues()):
            logpdf = self.pdfestimator.get_logpdf(res.phipsi)

            for meth, result in zip(self.methods, results):
                state_propensities, state_variability = meth.calculate(logpdf)
                result.add_entry(ConfStateResultsEntry(
                    res, state_propensities, state_variability))
                
        return results

