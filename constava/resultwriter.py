import abc
import csv
from .ensembles import ProteinEnsemble
from .calculator import ConfStateResults


class WriterABC(metaclass=abc.ABCMeta):

    def __init__(self, float_precision: int = 4):
        self.float_precision = float_precision

    def floatToString(self, x):
        return "{0:.{1:d}f}".format(x, self.float_precision)

    @abc.abstractmethod
    def writeToFile(self, results: ConfStateResults, output_file: str):
        pass


class CsvWriter(WriterABC):

    def writeToFile(self, results: ConfStateResults, output_file: str):
        with open(output_file, "w") as fhandle:
            writer = csv.writer(fhandle)
            writer.writerow([
                "#ResIndex", 
                "ResName",
                "ConStaVa",
                *results.state_labels])
            for entry in results.entries:
                writer.writerow([
                    entry.residue.respos,
                    entry.residue.restype,
                    self.floatToString(entry.state_variability),
                    *map(self.floatToString, entry.state_propensities)
                ])


class ResultWriter:

    def __init__(self, filetype_str: str = "csv", float_precision: int = 4):
        self.filetype_str = filetype_str
        self.float_precision = float_precision

    def writeToFile(self, results: ConfStateResults, output_file: str):
        strategy = self.get_strategy(self.filetype_str)
        strategy.writeToFile(results, output_file)
    
    def get_strategy(self, filetype: str) -> WriterABC:
        if filetype.lower() == "csv":
            return CsvWriter(self.float_precision)
        else:
            raise ValueError(f"Unknown argument for --input-format flag: `{filetype}`")
