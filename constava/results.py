from dataclasses import dataclass
from typing import List
import numpy as np
from .ensembles import ProteinEnsemble, ResidueEnsemble
from .methods import ConstavaABC

@dataclass
class ConfStateResultsEntry:
    residue: ResidueEnsemble = None
    state_propensities: np.ndarray = None
    state_variability: float = None


class ConfStateResults:
    """ A result class for ConfStateCalculator """
    def __init__(self, method: ConstavaABC, protein: ProteinEnsemble, state_labels: List[str], entries: List[ConfStateResultsEntry] = None):
        self.method = method
        self.protein = protein
        self.state_labels = state_labels
        self.entries = entries or []

    def add_entry(self, new_entry: ConfStateResultsEntry):
        self.entries.append(new_entry)
        self.entries.sort(key=lambda x: x.residue.respos)