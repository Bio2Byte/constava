"""constava.ensemblereader contains classes to read input data in different file
formats."""

import abc, os, re
from typing import List
from warnings import warn
import numpy as np
import pandas as pd
from ..datautils.ensembles import ProteinEnsemble, ResidueEnsemble


class UnknownFileStructureError(ValueError):
    """Raised if the file structure is not recognized"""
    pass

class DihedralRangeError(ValueError):
    """Raised if any dihedral angles are not correctly in radians"""
    pass

class DihedralRangeWarning(UserWarning):
    """Raised on the suspicion that  dihedral angles are not correctly in radians"""
    pass

def check_dihedral_range(arr: np.ndarray):
    """Helper method that checks if the dihedral angles are correctly in 
    radians.
    
    Parameters:
    -----------
        arr : Array[N,2]
            An array of N (phi, psi) pairs.

    Raises:
    -------
        DihedralRangeError
            If any dihedrals fall outside the range [-pi, pi]

        DihedralRangeWarning
            If all dihedrals fall in the range of [-(pi*pi/180), (pi*pi/180)], 
            as this suggests that angles were converted to radians twice.
    """
    vmin, vmax = np.min(arr), np.max(arr)
    if vmin < -np.pi or vmax > np.pi:
        raise DihedralRangeError(f"Dihedrals outside the range [-pi, pi] detected: [{vmin:.3f}, {vmax:.3f}]")
    elif vmin >= np.radians(np.pi) and vmax <= np.radians(np.pi):
        warn(("Provided dihedrals a very small: [{vmin:.3f}, {vmax:.3f}]. "
            "Please check that convertion to radians was only applied once."), 
            DihedralRangeWarning())

class ReaderABC(metaclass=abc.ABCMeta):
    """Base class for all file reader strategies"""

    def __init__(self, degrees2radians: bool = False):
        """Initializes the reader.
        
        Parameters:
        -----------
            degrees2radians : bool = False
                Sets if the read data should be converted from degrees to radians
                (default: Do not convert.)
        """
        self.degrees2radians = degrees2radians

    @abc.abstractclassmethod
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
        pass

    @abc.abstractmethod
    def readFiles(self, *input_files: str) -> ProteinEnsemble:
        """Reads the dihedral angles from one or more input files, converts them
        into radians (if degrees2radians == True) and returns them as a 
        ProteinEnsemble object
        
        Parameters:
        -----------
            *input_files : str
                One or more files to be read.

        Returns:
        --------
            prot : ProteinEnsemble
                Object that stores the dihedral angles for all the residues.
        """
        pass


class DihedralCsvReader(ReaderABC):
    """A reader strategy for csv files, as provided by constava.dihedrals.

    Attributes:
    -----------
        degrees2radians : bool
            Sets if read data should be converted from degrees to radians.

    Methods:
    --------
        checkFileFormat(*input_files)
            Checks if the files provided adhere to the given format.
        readFiles(*input_files)
            Reads the dihedral angles from one or more input files.
    """
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
        for infile in input_files:
            with open(infile, "r") as f:
                first_line = f.readline()
                if not first_line.startswith("#Frame,ResIndex,ResName,Phi[rad],Psi[rad]"):
                    return False
        return True

    def readFiles(self, *input_files) -> ProteinEnsemble:
        """Reads the dihedral angles from one or more input files, converts them
        into radians (if degrees2radians == True) and returns them as a 
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
    """A reader strategy designed to read the output of GROMACS' chi module 
    command: `gmx chi -s <structure> -f <trajectory> -rama`. The files are
    generally named `ramaPhiPsi[RESNAME][RESINDEX].xvg`. Both RESNAME and
    RESINDEX are directly extracted from the filename.

    Attributes:
    -----------
        FILENAME_REGEX : re.Pattern
            Regular expression to extract RESNAME and RESINDEX from the filename.
        degrees2radians : bool
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
        into radians (if degrees2radians == True) and returns them as a 
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
            if self.degrees2radians:
                phipsi = np.radians(phipsi)
            check_dihedral_range(phipsi) # Throws errors/warnings if data is not in radians
            residue_list.append(
                ResidueEnsemble(restype=restype, respos=respos, phipsi=phipsi))
        return ProteinEnsemble(residue_list)


class EnsembleReader:
    """Interface class for the file reader strategies. This is the only class
    a user should interact with.
    
    Attributes:
    -----------
        filetype_str : str {`auto`, `xvg`, `csv`}
            A string to indicate which reader strategy should be used.
        degrees2radians : bool
            Sets if read data should be converted from degrees to radians.

    Methods:
    --------
        checkFileFormat(*input_files)
            Checks if the files provided adhere to the given format.
        readFiles(*input_files)
            Reads the dihedral angles from one or more input files.
        get_strategy(self, filetype, input_files)
            Returns a reader strategy based on the `filetype_str` or guessing.
        guess_strategy(input_files)
            Guesses a reader strategy based on the structure/filenames of the `input_files`.
    """

    def __init__(self, filetype_str: str = "auto", degrees2radians: bool = False):
        self.filetype_str = filetype_str
        self.degrees2radians = degrees2radians

    def readFiles(self, *input_files: str) -> ProteinEnsemble:
        """Reads the dihedral angles from one or more input files using an 
        appropriate strategy. It further converts the dihedrals into radians 
        (if degrees2radians == True) and returns them as a ProteinEnsemble object.
        
        Parameters:
        -----------
            *input_files : str
                One or more files to be read.

        Returns:
        --------
            prot : ProteinEnsemble
                Object that stores the dihedral angles for all the residues.

        Raises:
        -------
            UnknownFileStructureError
                If no appropriate reader strategy was found.
        """
        strategy = self.get_strategy(self.filetype_str, input_files)
        return strategy.readFiles(*input_files)
    
    def get_strategy(self, filetype: str, input_files: List[str]) -> ReaderABC:
        """Method to return an appropriate reader strategy. This is primarily
        based on the `filetype_str`. If `auto` was provieded the reader strategy
        is guessed based on file structures and filenames.

        Parameters:
        -----------
            filetype : str
                String indicating the reader strategy to be used.
            input_files : List[str]
                Lst of one or more files to be read (used for guessing filetype).
        
        Returns:
        --------
            strategy : ReaderABC
                Reader strategy that is to be used to parse the files.
        """
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
        """Method to try to guess an appropriate reader strategy, based on file
        structures and file names.

        Parameters:
        -----------
            *input_files : str
                One or more files to be read.

        Returns:
        --------
            strategy : ReaderABC
                Reader strategy that is to be used to parse the files.

        Raises:
        -------
            UnknownFileStructureError
                If no appropriate reader strategy was found.
        """
        # Check for the internal CSV-format
        if DihedralCsvReader.checkFileFormat(*input_files):
            return DihedralCsvReader
        # Check if the files correspond to the output of `gmx chi`
        elif GmxChiReader.checkFileFormat(*input_files):
            return GmxChiReader
        # IF all checks fail, raise an UnknownFileStructureError
        else:
            raise UnknownFileStructureError()