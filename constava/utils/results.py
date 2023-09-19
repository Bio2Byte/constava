""" constava.results contains classes that hold output information.

ConstavaResultsEntry stores state_propensities and state_varaibility for 
a single ResidueEnsemble.

ConstavaResults hold multiple ConstavaResultsEntry objects, thus describing 
a while ProteinEnsemble.
"""


from dataclasses import dataclass
from typing import List
import numpy as np
from .ensembles import ProteinEnsemble, ResidueEnsemble, EnsembleMismatchError
from ..calc.subsampling import SubsamplingABC

@dataclass
class ConstavaResultsEntry:
    """Results for a single ResidueEnsemble.
    
    Attributes:
    -----------
        residue : ResidueEnsemble
            The residue ensemble for which the results were calculated

        state_propensities : Array[N] or Array[M,N]
            The calculated conformational state propensities for N states. If 
            time series information is calculated, this is a 2-dimensional array
            of M samples * N conformational states.
        
        state_variability : float or Array[M]
            The calculated conformational state variability across the models.
            If time series information is calculated, this is an array of 
            M samples.
    """
    residue: ResidueEnsemble = None
    state_propensities: np.ndarray = None
    state_variability: float or np.ndarray = None


class ConstavaResults:
    """ Results from ConfStateCalculator for a given ProteinEnsemble and Method 

    Attributes:
    -----------
        method : str
            A string indicating the subsampling method with which the results 
            were calculated

        protein : ProteinEnsemble
            The conformational ensemble for which all results were calculated

        state_labels : List[str]
            The labels of the conformational states, which should be human-readable.
        
        entries : List[ConstavaResultsEntry]
            A list of result entries, each of which contains the information
            for a single residue in the conformational ensemble.

    Methods:
    --------
        add_entry(new_entry)
            Adds an additional entry to the result class
    """
    def __init__(self, method: str, protein: ProteinEnsemble, state_labels: List[str], entries: List[ConstavaResultsEntry] = None):
        self.method = method
        self.protein = protein
        self.state_labels = state_labels
        self.entries = entries or []
    
    def __repr__(self):
        return f"Results(method={self.method}, {len(self.entries)} entries)"

    def add_entry(self, new_entry: ConstavaResultsEntry):
        """Adds an additional entry to the result class"""
        if new_entry.residue.protein is not self.protein:
            raise EnsembleMismatchError(f"Result for residue {new_entry.residue} does not belong to the given Protein")
        self.entries.append(new_entry)
        self.entries.sort(key=lambda x: x.residue.respos)