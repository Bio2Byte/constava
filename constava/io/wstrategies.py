import abc
import csv
import json
from datetime import datetime
from typing import List, TextIO 
import numpy as np
from ..utils.constants import CONSTAVA_NAME, CONSTAVA_VERSION
from ..utils.results import ConstavaResults, ConstavaResultsEntry


class WriterABC(metaclass=abc.ABCMeta):
    """Abstract base class for file writers"""
    def __init__(self, float_precision: int = 4):
        self.float_precision = float_precision

    @abc.abstractmethod
    def write_results(self, fhandle: TextIO, results: List[ConstavaResults]):
        pass


class JsonWriter(WriterABC):
    """Writer strategy for writing Json files"""

    def write_results(self, fhandle: TextIO, results: List[ConstavaResults]):
        __dict = {
            "tool": CONSTAVA_NAME,
            "version": CONSTAVA_VERSION,
            "creation_date": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
            "results": {
                "ResIndex": [entry.residue.respos for entry in results[0].entries],
                "ResName": [entry.residue.restype for entry in results[0].entries],
            }
        }
        for result in results:
            cspropensites = np.round(np.stack([entry.state_propensities for entry in result.entries]), decimals=self.float_precision)
            csvariability = np.round(np.stack([entry.state_variability for entry in result.entries]), decimals=self.float_precision)
            __dict["results"][result.method] = {
                state_label: cspropensites[:,i].tolist() 
                for i, state_label in enumerate(result.state_labels)
            }
            __dict["results"][result.method]["Variability"] = csvariability.tolist()

        json.dump(__dict, fhandle)
        
        
class CsvWriter(WriterABC):
    """Writer strategy for writing CSV files"""

    csvdialect = {}
    
    def write_results(self, fhandle: TextIO, results: List[ConstavaResults]):
        # Initialize csv writer
        writer = csv.writer(fhandle, **self.csvdialect)
        # Write the header line
        writer.writerow([
            "#Method", "SeriesIndex", "ResIndex", "ResName", 
            *results[0].state_labels, "Variability"])
        # Iterate over results (methods) and entries (residues)
        for result in results:
            for entry in result.entries:
                writer.writerows(self._write_entry(entry, result.method))

    def _write_entry(self, entry: ConstavaResultsEntry, method: str) -> List:
        flt2str = lambda x: "{0:.{1:d}f}".format(x, self.float_precision)
        respos, restype = entry.residue.respos, entry.residue.restype
        # If average state propensities are reported, one row is written
        if len(entry.state_propensities.shape) == 1:
            values = np.concatenate([entry.state_propensities, [entry.state_variability]])
            yield (method, None, respos, restype, *map(flt2str, values))
        # If state propensities are reported as series, multiple rows are written
        else:
            arr = np.concatenate([entry.state_propensities, [entry.state_variability]]).T
            for i, values in enumerate(arr):
                yield (method, i, respos, restype, *map(flt2str, values))


class TsvWriter(CsvWriter):
    """Writer strategy for writing TSV files (based on CsvWriter)"""
    csvdialect = {"delimiter": "\t"}