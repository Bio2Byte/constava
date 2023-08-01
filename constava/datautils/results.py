""" constava.results contains classes that hold output information.

ConfStateResultsEntry stores state_propensities and state_varaibility for 
a single ResidueEnsemble.

ConfStateResults hold multiple ConfStateResultsEntry objects, thus describing 
a while ProteinEnsemble.
"""


from dataclasses import dataclass
from typing import List
import numpy as np
from .ensembles import ProteinEnsemble, ResidueEnsemble
from ..calc.subsampling import SubsamplingABC

@dataclass
class ConfStateResultsEntry:
    """ Results for a single ResidueEnsemble"""
    residue: ResidueEnsemble = None
    state_propensities: np.ndarray = None
    state_variability: float = None


class ConfStateResults:
    """ Results from ConfStateCalculator for a given ProteinEnsemble and Method """
    def __init__(self, method: SubsamplingABC, protein: ProteinEnsemble, state_labels: List[str], entries: List[ConfStateResultsEntry] = None):
        self.method = method
        self.protein = protein
        self.state_labels = state_labels
        self.entries = entries or []

    def add_entry(self, new_entry: ConfStateResultsEntry):
        self.entries.append(new_entry)
        self.entries.sort(key=lambda x: x.residue.respos)