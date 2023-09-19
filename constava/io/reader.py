"""constava.reader contains the reader interface to read input data"""

from typing import List
from .rstrategies import ReaderABC, DihedralCsvReader, GmxChiReader, UnknownFileStructureError
from ..utils.ensembles import ProteinEnsemble


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
    STRATEGIES = {
        "csv": DihedralCsvReader,
        "xvg": GmxChiReader,
    }

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
        filetype = (filetype or "").lower()
        if filetype in self.STRATEGIES:
            return self.STRATEGIES[filetype](degrees=self.degrees2radians)
        elif filetype in ("auto", "guess", ""):
            return self.guess_strategy(input_files, degrees=self.degrees2radians)
        else:
            raise ValueError(f"Unknown argument for --input-format flag: `{filetype}`")

    @classmethod
    def guess_strategy(cls, input_files: List[str], degrees=False) -> ReaderABC:
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
        for Reader in cls.STRATEGIES.values():
            reader = Reader(degrees=degrees)
            if reader.checkFileFormat(*input_files):
                return reader
        else:
            # IF all checks fail, raise an UnknownFileStructureError
            raise UnknownFileStructureError("Dihedral input corresponds to no known input format.")