"""constava.writer contains classes the write the output to different 
data formats."""

import os 
from typing import Any
from .wstrategies import JsonWriter, CsvWriter, TsvWriter
from ..utils.results import ConstavaResults


class ResultsWriter:

    STRATEGIES = {
        "json": JsonWriter,
        "csv": CsvWriter,
        "tsv": TsvWriter,
    }

    def __init__(self, filename: str, format: str = "auto", float_precision: int = 4):
        self.filename = filename
        self.format = format
        self.float_precision =  float_precision
        self._strategy = self.get_strategy()

    def __setattr__(self, __name: str, __value: Any) -> None:
        super().__setattr__(__name, __value)
        # Update the strategy
        if __name in ("filename", "format", "float_precision") and hasattr(self, "_strategy"):
            super().__setattr__(self._strategy, self.get_strategy())
        
    def get_strategy(self):
        if (self.format or "").lower() in ["auto", "guess", ""]:
            __fmt = os.path.splitext(self.filename)[1].lstrip(".")
        else:
            __fmt = self.format
        if __fmt in self.STRATEGIES:
            return self.STRATEGIES[__fmt](self.float_precision)
        else:
            raise ValueError(f"Unknown output file format: `{__fmt}`")

    def write_results(self, results: ConstavaResults):
        with open(self.filename, "w") as fhandle:
            self._strategy.write_results(fhandle, results)