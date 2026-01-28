"""constava.writer contains classes the write the output to different
data formats."""

import os
from typing import Any
from .wstrategies import JsonWriter, CsvWriter, TsvWriter
from ..utils.results import ConstavaResults


class ResultsWriter:
    """
    ResultsWriter writes the ConstavaResults to a file using different strategies.
    """

    STRATEGIES = {
        "json": JsonWriter,
        "csv": CsvWriter,
        "tsv": TsvWriter,
    }

    def __init__(
        self,
        filename: str,
        file_format: str = "auto",
        float_precision: int = 4,
        indent_size: int = 0,
    ):
        self.filename = filename
        self.file_format = file_format
        self.float_precision = float_precision
        self.indent_size = indent_size
        self._strategy = self.get_strategy()

    def __setattr__(self, __name: str, __value: Any) -> None:
        super().__setattr__(__name, __value)
        # Update the strategy
        if __name in ("filename", "file_format", "float_precision", "indent_size") and hasattr(
            self, "_strategy"
        ):
            super().__setattr__(self._strategy, self.get_strategy())

    def get_strategy(self):
        """
        Getter of the Writer strategy based on the file_format attribute (file_format).
        """

        if (self.file_format or "").lower() in ["auto", "guess", ""]:
            __fmt = os.path.splitext(self.filename)[1].lstrip(".")
        else:
            __fmt = self.file_format
        if __fmt in self.STRATEGIES:
            return self.STRATEGIES[__fmt](self.float_precision, self.indent_size)
        else:
            raise ValueError(f"Unknown output file format: `{__fmt}`")

    def write_results(self, results: ConstavaResults):
        """
        Write the results to the output file using the `filename` attribute as filename.

        :param results: The results to write to a file.
        :type results: ConstavaResults
        """

        with open(self.filename, "w", encoding="utf-8") as fhandle:
            self._strategy.write_results(fhandle, results)
