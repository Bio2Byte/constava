import abc, os, re
from warnings import warn
import numpy as np
import pandas as pd

from ..utils.ensembles import ProteinEnsemble, ResidueEnsemble
from ..utils.utils import check_dihedral_range


class UnknownFileStructureError(ValueError):
    """Raised if the file structure is not recognized"""
    pass


class ReaderABC(metaclass=abc.ABCMeta):
    """Base class for all file reader strategies"""

    def __init__(self, degrees: bool = False):
        self.degrees = degrees

    @abc.abstractmethod
    def checkFileFormat(self, *input_files: str) -> bool:
        """Checks if the files provided adhere to the given format."""
        pass

    @abc.abstractmethod
    def readFiles(self, *input_files: str) -> ProteinEnsemble:
        """Reads the dihedral angles from one or more input files, converts them
        into radians (if degrees == True) and returns them as a 
        ProteinEnsemble object
        """
        pass


class DihedralCsvReader(ReaderABC):
    """A reader strategy for csv files, as provided by constava.dihedrals.

    Attributes:
    -----------
        degrees : bool
            Sets if read data should be converted from degrees to radians.

    Methods:
    --------
        checkFileFormat(*input_files)
            Checks if the files provided adhere to the given format.
        readFiles(*input_files)
            Reads the dihedral angles from one or more input files.
    """

    def checkFileFormat(self, *input_files: str) -> bool:
        """Checks if the files provided adhere to the given format.
        
        Parameters:
        -----------
            *input_files : str
                One or more files to be read.

        Returns:
        --------
            check_ok : bool
                True if all provided files adhere to the file format, else False.
        """
        phicol = "Phi[{0}]".format("deg" if self.degrees else "rad")
        psicol = "Psi[{0}]".format("deg" if self.degrees else "rad")
        expected_columns = {"#Frame", "ResIndex", "ResName", phicol, psicol}
        for infile in input_files:
            with open(infile, "r") as f:
                first_line = f.readline()
            columns = set(first_line.strip().split(","))
            if expected_columns.intersection(columns) != expected_columns:
                return False
        return True

    def readFiles(self, *input_files) -> ProteinEnsemble:
        """Reads the dihedral angles from one or more input files, converts them
        into radians (if degrees == True) and returns them as a 
        ProteinEnsemble object.
        
        Parameters:
        -----------
            *input_files : str
                One or more files to be read.

        Returns:
        --------
            prot : ProteinEnsemble
                Object that stores the dihedral angles for all the residues.
        """
        phicol = "Phi[{0}]".format("deg" if self.degrees else "rad")
        psicol = "Psi[{0}]".format("deg" if self.degrees else "rad")
        data = pd.DataFrame()
        for infile in input_files:
            data = pd.concat([data, pd.read_csv(infile)], axis=0)
        unique_residues = data[["ResIndex", "ResName"]].drop_duplicates(inplace=False)
        residue_list = []
        for i, (resid, resname) in unique_residues.iterrows():
            phipsi = data.loc[(data["ResIndex"] == resid) & (data["ResName"] == resname), [phicol, psicol]].to_numpy()
            if self.degrees:
                phipsi = np.radians(phipsi)
            check_dihedral_range(phipsi) # Throws errors/warnings if data is not in radians
            residue_list.append(
                ResidueEnsemble(restype=resname, respos=resid, phipsi=phipsi))
        return ProteinEnsemble(residue_list)


class GmxChiReader(ReaderABC):
    """A reader strategy designed to read the output of GROMACS' chi module 
    command: `gmx chi -s <structure> -f <trajectory> -rama`. The files are
    generally named `ramaPhiPsi[RESNAME][RESINDEX].xvg`. Both RESNAME and
    RESINDEX are directly extracted from the filename.

    Attributes:
    -----------
        FILENAME_REGEX : re.Pattern
            Regular expression to extract RESNAME and RESINDEX from the filename.
        degrees : bool
            Sets if read data should be converted from degrees to radians.
    
    Methods:
    --------
        checkFileFormat(*input_files)
            Checks if the files provided adhere to the given format.
        readFiles(*input_files)
            Reads the dihedral angles from one or more input files.
    """
    FILENAME_REGEX = re.compile("ramaPhiPsi([A-Z][A-Z0-9][A-Z0-9])([0-9]+).xvg")

    @classmethod
    def checkFileFormat(cls, *input_files: str) -> bool:
        """Checks if the files provided adhere to the given format.
        
        Parameters:
        -----------
            *input_files : str
                One or more files to be read.

        Returns:
        --------
            check_ok : bool
                True if all provided files adhere to the file format, else False.
        """
        for filename in (os.path.basename(fp) for fp in input_files):
            if not cls.FILENAME_REGEX.match(filename):
                return False
        return True
    
    def readFiles(self, *input_files) -> ProteinEnsemble:
        """Reads the dihedral angles from one or more input files, converts them
        into radians (if degrees == True) and returns them as a 
        ProteinEnsemble object.
        
        Parameters:
        -----------
            *input_files : str
                One or more files to be read.

        Returns:
        --------
            prot : ProteinEnsemble
                Object that stores the dihedral angles for all the residues.
        """
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
            if self.degrees:
                phipsi = np.radians(phipsi)
            check_dihedral_range(phipsi) # Throws errors/warnings if data is not in radians
            residue_list.append(
                ResidueEnsemble(restype=restype, respos=respos, phipsi=phipsi))
        return ProteinEnsemble(residue_list)