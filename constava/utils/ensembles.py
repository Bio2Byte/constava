"""constava.ensembles contains classes describing the input data from the 
conformational ensemble.
"""

from dataclasses import dataclass
from typing import List, Generator
import numpy as np

from .constants import aminoacids3to1


class EnsembleMismatchError(ValueError):
    """Error raised when ResidueEnsembles from different ProteinEnsembles are mixed"""
    pass


@dataclass
class ResidueEnsemble:
    """ 
    Dataclass to hold information on a given residue in a conformational ensemble:

    Attributes:
    -----------
        restype: str            Residue name of the residue
        respos: int             Index of the residue in the sequence
        phipsi: array[N,2]      Array of the phi/psi angles of the residue
                                for all conformations in the ensemble
        proteinensemble: ProteinEnsemble
                                Refernece to the ProteinEnsemble-object, the 
                                residue belongs to.
    """
    restype: str = ""
    respos: int = None
    phipsi: np.ndarray = None
    protein = None

    @property
    def restype1(self):
        """ Converts the restype attribute to one-letter code """
        return aminoacids3to1.get(self.restype, "X")

    def __repr__(self):
        """ Short representation of the object """
        return f"<{self.restype}:{self.respos}>"
    
    def __lt__(self, other):
        """ Implemented to allow for easy sorting of residues """
        if not isinstance(other, self.__class__):
            raise TypeError(f"Cannot compare {self.__class__.__name__} with {other.__class__.__name__}")
        return self.respos < other.respos
    
    def to_dict(self):
        """ Returns the residue data as a dictionary (e.g., for later conversion to JSON) """
        _data = {
            "restype":  self.restype,
            "respos":   self.respos,
            "phipsi":   (self.phipsi.tolist() if isinstance(self.phipsi, np.ndarray) else None),
        }
        return _data


class ProteinEnsemble:
    """ 
    Class to hold information on conformational ensemble. It contains multiple
    ResidueEnsemble objects which hold the actual data per residue
    """
    def __init__(self, residues: List[ResidueEnsemble] = None):
        """ Constructor for the ProteinEnsemble class 
        
        Parameters:
        -----------
            residues: List[ResidueEnsemble]     A list of ResidueEnsemble objects 
                                                to be added to the ensemble. 
        """
        self._residues = []
        self.add_residues(*residues)

    def __repr__(self):
        """ Short string representation of a class-object """
        return f"<ProteinEnsemble: {self.n_residues} residues>"
    
    @property
    def n_residues(self):
        """ Returns the number of residues form the first to last residue 
        this might include gaps (residues without data) """
        return 1 + self._residues[-1].respos - self._residues[0].respos

    @property
    def sequence(self):
        """ Returns the protein sequence as a string """
        seq = self._getPropertyFromResidues("restype1", fillvalue="-")
        return "".join(seq)
    
    def _getPropertyFromResidues(self, resattr: str, *, fillvalue=np.nan, dtype=None):
        """ Private helper function, to retrieve properties from the individual
        residues across the whole protein. """
        offset = self._residues[0].respos
        result = None
        for i, value in ((res.respos-offset, getattr(res, resattr)) for res in self.get_residues()):
            if value is None:
                continue
            elif result is None and isinstance(value, np.ndarray):
                result = np.full((self.n_residues,)+value.shape, fillvalue)
            elif result is None:
                result = np.full((self.n_residues,), fillvalue)
            result[i] = value
        return result
    
    def get_residues(self) -> Generator[ResidueEnsemble, None, None]:
        """ Returns a generator for all residues in the class """
        return (res for res in self._residues)
    
    def add_residues(self, *new_residues: ResidueEnsemble):
        """ Adds a new residue to the ensemble, and sorts the residue list 
        according to their indices """
        for res in new_residues:
            res.protein = self
            self._residues.append(res)
        self._residues.sort()
        # self.__dict__.pop('sequence', None)
        # self.__dict__.pop('dynamics', None)
        # self.__dict__.pop('conformation', None)

    def to_dict(self):
        _data = {
            "residues": [res.to_dict() for res in self.get_residues()]
        }
        return _data

