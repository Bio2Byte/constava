import abc, os, re
from typing import List
from warnings import warn
import numpy as np
import pandas as pd
from .ensembles import ProteinEnsemble, ResidueEnsemble


class UnknownFileStructureError(ValueError):
    pass

class DihedralRangeError(ValueError):
    pass

class DihedralRangeWarning(UserWarning):
    pass

def check_dihedral_range(arr: np.ndarray) -> bool:
    vmin, vmax = np.min(arr), np.max(arr)
    if vmin < -np.pi or vmax > np.pi:
        raise DihedralRangeError(f"Dihedrals outside the range [-pi, pi] detected: [{vmin:.3f}, {vmax:.3f}]")
    elif vmin >= np.radians(np.pi) and vmax <= np.radians(np.pi):
        warn(("Provided dihedrals a very small: [{vmin:.3f}, {vmax:.3f}]. "
            "Please check that convertion to radians was only applied once."), 
            DihedralRangeWarning())

class ReaderABC(metaclass=abc.ABCMeta):

    def __init__(self, degrees2radians: bool = False):
        self.degrees2radians = degrees2radians

    @abc.abstractclassmethod
    def checkFileFormat(self, *input_files: str) -> bool:
        """ 
        Checks if the files provided adhere to the given format. Returns True 
        if all provided files adhere to the file format, else False is returned.
        """
        pass

    @abc.abstractmethod
    def readFiles(self, *input_files: str) -> ProteinEnsemble:
        """ Reads the input dihedrals and returns them as a ProteinEnsemble object """
        pass


class DihedralCsvReader(ReaderABC):
    """
    """
    @classmethod
    def checkFileFormat(cls, *input_files: str) -> bool:
        for infile in input_files:
            with open(infile, "r") as f:
                first_line = f.readline()
                if not first_line.startswith("#Frame,ResIndex,ResName,Phi[rad],Psi[rad]"):
                    return False
        return True

    def readFiles(self, *input_files) -> ProteinEnsemble:
        # Read all input files
        dihedrals = pd.DataFrame()
        for infile in input_files:
            dihedrals = pd.concat([dihedrals, pd.read_csv(infile)], axis=0)
        unique_residues = dihedrals[["ResIndex", "ResName"]].drop_duplicates(inplace=False)
        residue_list = []
        for i, (resid, resname) in unique_residues.iterrows():
            phipsi = dihedrals.loc[
                (dihedrals["ResIndex"] == resid) & (dihedrals["ResName"] == resname),
                ["Phi[rad]", "Psi[rad]"]].to_numpy()
            if self.degrees2radians:
                phipsi = np.radians(phipsi)
            check_dihedral_range(phipsi) # Throws errors/warnings if data is not in radians
            residue_list.append(
                ResidueEnsemble(restype=resname, respos=resid, phipsi=phipsi))
        return ProteinEnsemble(residue_list)


class GmxChiReader(ReaderABC):
    """
    A reader that to parse the output of:
        gmx chi -s <structure> -f <trajectory> -rama
    """

    FILENAME_REGEX = re.compile("ramaPhiPsi([A-Z][A-Z0-9][A-Z0-9])([0-9]+).xvg")

    @classmethod
    def checkFileFormat(cls, *input_files: str) -> bool:
        for filename in (os.path.basename(fp) for fp in input_files):
            if not cls.FILENAME_REGEX.match(filename):
                return False
        return True
    
    def readFiles(self, *input_files) -> ProteinEnsemble:
        residue_list = []
        for infile in input_files:
            filename = os.path.basename(infile)
            m = self.FILENAME_REGEX.match(filename)
            if m is None:
                raise UnknownFileStructureError((
                    f"File does not match the structure of `gmx chi` outputs: {filename}"))
            restype = m.group(1)
            respos = int(m.group(2))
            phipsi = np.loadtxt(infile, comments=["#", "@"])
            if self.degrees2radians:
                phipsi = np.radians(phipsi)
            check_dihedral_range(phipsi) # Throws errors/warnings if data is not in radians
            residue_list.append(
                ResidueEnsemble(restype=restype, respos=respos, phipsi=phipsi))
        return ProteinEnsemble(residue_list)


class EnsembleReader:

    def __init__(self, filetype_str: str = "auto", degrees2radians: bool = False):
        self.filetype_str = filetype_str
        self.degrees2radians = degrees2radians

    def readFiles(self, *input_files: str) -> ProteinEnsemble:
        strategy = self.get_strategy(self.filetype_str, input_files)
        return strategy.readFiles(*input_files)
    
    def get_strategy(self, filetype: str, input_files: List[str]) -> ReaderABC:
        if filetype.lower() == "xvg":
            return GmxChiReader(degrees2radians=self.degrees2radians)
        elif filetype.lower() == "csv":
            return DihedralCsvReader(degrees2radians=self.degrees2radians)
        elif filetype.lower() == "auto":
            return self.guess_strategy(input_files)(degrees2radians=self.degrees2radians)
        else:
            raise ValueError(f"Unknown argument for --input-format flag: `{filetype}`")

    @classmethod
    def guess_strategy(cls, input_files: List[str]) -> ReaderABC:
        # Check for the internal CSV-format
        if DihedralCsvReader.checkFileFormat(*input_files):
            return DihedralCsvReader
        # Check if the files correspond to the output of `gmx chi`
        elif GmxChiReader.checkFileFormat(*input_files):
            return GmxChiReader
        # IF all checks fail, raise an UnknownFileStructureError
        else:
            raise UnknownFileStructureError()