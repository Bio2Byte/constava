from dataclasses import dataclass, field
from typing import List, Generator
import numpy as np

from .constants import aminoacids3to1


@dataclass
class ResidueEnsemble:

    restype: str = ""
    respos: int = None
    phipsi: np.ndarray = None
    proteinensemble = None

    @property
    def restype1(self):
        return aminoacids3to1.get(self.restype, "X")

    def __repr__(self):
        return f"<{self.restype}:{self.respos}>"
    
    def to_dict(self):
        _data = {
            "restype":  self.restype,
            "respos":   self.respos,
            "phipsi":   (self.phipsi.tolist() if isinstance(self.phipsi, np.ndarray) else None),
        }
        return _data


class ProteinEnsemble:

    def __init__(self, residues: List[ResidueEnsemble] = None):
        self._residues = []
        self.add_residues(*residues)

    def __repr__(self):
        return f"<ProteinEnsemble: {self.n_residues} residues>"
    
    @property
    def n_residues(self):
        return 1 + self._residues[-1].respos - self._residues[0].respos

    @property
    def sequence(self):
        seq = self._getPropertyFromResidues("restype1", fillvalue="-")
        return "".join(seq)
    
    # @property
    # def conformational_state_variablity(self):
    #     return self._getPropertyFromResidues("conformational_state_variablity")
    
    # @property
    # def conformational_state_propensities(self):
    #     return self._getPropertyFromResidues("conformational_state_propensities")
    
    def _getPropertyFromResidues(self, resattr: str, *, fillvalue=np.nan, dtype=None):
        offset = self._residues[0].respos
        result = None
        for i, value in ((res.respos-offset, getattr(res, resattr)) for res in self.get_residues()):
            if value is None:
                continue
            elif result is None and isinstance(value, np.ndarray):
                result = np.full((self.n_residues,)+value.shape, fillvalue)
            else:
                result = np.full((self.n_residues,), fillvalue)
            result[i] = value
        return result
    
    def get_residues(self) -> Generator[ResidueEnsemble, None, None]:
        return (res for res in self._residues)
    
    def add_residues(self, *new_residues: ResidueEnsemble):
        for res in new_residues:
            res.owner = self
            self._residues.append(res)
        self._residues.sort(key=lambda res: res.respos)
        # self.__dict__.pop('sequence', None)
        # self.__dict__.pop('dynamics', None)
        # self.__dict__.pop('conformation', None)

    def to_dict(self):
        _data = {
            "residues": [res.to_dict() for res in self.get_residues()]
        }
        return _data

