import abc
import csv
import json
import os 
from typing import List
import numpy as np
from .results import ConfStateResults, ConfStateResultsEntry


class WriterABC(metaclass=abc.ABCMeta):

    def __init__(self, float_precision: int = 4):
        self.float_precision = float_precision

    def floatToString(self, x):
        return "{0:.{1:d}f}".format(x, self.float_precision)

    @abc.abstractmethod
    def writeToFile(self, results: ConfStateResults, output_file: str):
        pass


class CsvWriter(WriterABC):

    def writeToFile(self, results_list: List[ConfStateResults], output_file: str):
        with open(output_file, "w") as fhandle:
            writer = csv.writer(fhandle)
            # Write the header line
            writer.writerow(["#ResIndex", "ResName", "ConStaVa", 
                             *results_list[0].state_labels, "Method"])
            # Iterate over results (methods) and write one line per reside
            for results in results_list:
                self._writeEntries(writer, 
                    results.entries, method=results.method.getShortName())
    
    def _writeEntries(self, writer: csv.writer, entries_list: List[ConfStateResultsEntry], method: str = "N/A"):
        for entry in entries_list:
            writer.writerow([
                entry.residue.respos,
                entry.residue.restype,
                self.floatToString(entry.state_variability),
                *map(self.floatToString, entry.state_propensities),
                method
            ])


class JsonWriter(WriterABC):

    def writeToFile(self, results_list: List[ConfStateResults], output_file: str):
        output_data = []
        for results in results_list:
            output_data.append({
                "method":  results.method.getShortName(),
                "sequence": results.protein.sequence,
                "results": self._reportResults(results)
            })
        with open(output_file, "w") as fhandle:
            json.dump(output_data, fhandle)

    def _reportResults(self, results: ConfStateResults):
        datakeys = ["ResIndex", "ResName", "ConStaVa", *results.state_labels]
        result_dict = {key: [] for key in datakeys}
        for entry in results.entries:
            result_dict["ResIndex"].append(entry.residue.respos)
            result_dict["ResIndex"].append(entry.residue.restype)
            result_dict["ConStaVa"].append(entry.state_variability)
            for key, value in zip(results.state_labels, np.round(entry.state_propensities, decimals=self.float_precision)):
                result_dict[key].append(value)
        return result_dict


class ResultWriter:

    def __init__(self, filetype_str: str = "auto", float_precision: int = 4):
        self.filetype_str = filetype_str
        self.float_precision = float_precision

    def writeToFile(self, results: ConfStateResults, output_file: str):
        strategy = self.get_strategy(self.filetype_str, output_file)
        strategy.writeToFile(results, output_file)
    
    def get_strategy(self, filetype: str, output_file: str) -> WriterABC:
        if filetype.lower() == "csv":
            return CsvWriter(self.float_precision)
        elif filetype.lower() == "json":
            return JsonWriter(self.float_precision)
        elif filetype.lower() == "auto":
            return self.guess_strategy(output_file)
        else:
            raise ValueError(f"Unknown argument for --input-format flag: `{filetype}`")

    def guess_strategy(self, output_file: str):
        filename, fileext = os.path.splitext(output_file)
        if fileext.lower() == ".csv":
            return CsvWriter(self.float_precision)
        elif fileext.lower() == ".json":
            return JsonWriter(self.float_precision)
        else:
            raise ValueError(f"Cannot guess output format for file: `{output_file}`")