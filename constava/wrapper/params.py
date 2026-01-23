from dataclasses import dataclass, field
import typing
from ..utils.logging import logging

logger = logging.getLogger("Constava")

def set_logger_level(func):
    """Set logger.level while modifying verbose-attribute"""
    def _inner_(self, __attr, __value):
        rvalue = func(self, __attr, __value)
        if __attr == "verbose":
            if __value == 0:
                logger.setLevel(logging.WARNING)
            elif __value == 1:
                logger.setLevel(logging.INFO)
            else:
                logger.setLevel(logging.DEBUG)
        return rvalue
    return _inner_

def set_single_as_list(func):
    """Allows list-attributes (e.g., window) to be set with a single value"""
    def _inner_(self, __attr, __value):
        # ADIAZ: It works for <=3.13
        # dtype = typing.get_type_hints(self).get(__attr, None)
        hints = typing.get_type_hints(type(self))      # or self.__class__
        dtype = hints.get(__attr, None)
        
        if typing.get_origin(dtype) is list:
            if __value is None:
                __value = []
            elif not isinstance(__value, typing.Iterable) or isinstance(__value, str):
                __value = [__value]
        return func(self, __attr, __value)
    return _inner_


@dataclass
class ConstavaParameters:
    """The parameters that govern the function of Constava
    
    Parameters:
    -----------
        input_files : List[str] or str
            Input file(s) that contain the dihedral angles.
        input_format : str
            Format of the input file. Options include:
            - 'auto': Automatically detect the file format (default).
            - 'csv': Comma-separated values format.
            - 'xvg': XVG format used by GROMACS for graphing.
        output_file : str
            The file to write the output to.
        output_format : str
            Format of the output file. Options include:
            - 'auto': Automatically select the output format based on the input format or other criteria (default).
            - 'csv': Comma-separated values format, suitable for spreadsheets and simple data analyses.
            - 'json': JSON format, which is lightweight and easy for humans to read and write, and easy for machines
            to parse and generate.
            - 'tsv': Tab-separated values format, useful for tabular data that is less complex than CSV data.
        model_type : str
            Specifies the probabilistic conformational state model used. Options include:
            - 'kde': Kernel Density Estimator (default).
            - 'grid': A grid-based approximation of the KDE. Runs significantly faster with minor sacrifice to accuracy.
        model_load : str
            Load a conformational state model from the given pickled file.
        model_data : str
            Fit conformational state models to data provided in the given file.
        model_dump : str
            Write the generated model to a pickled file, that can be loaded
            again using `model_load`.
        window : List[int] or int
            Do inference using a moving reading-frame of <int> consecutive samples.
            Multiple values can be given as a list.
        window_series : List[int] or int
            Do inference using a moving reading-frame of <int> consecutive samples.
            Return the results for every window rather than the average. Multiple
            values can be given as a list.
        bootstrap : List[int] or int
            Do inference using <Int> samples obtained through bootstrapping.
            Multiple values can be given as a list.
        bootstrap_series : List[int] or int
            Do inference using <Int> samples obtained through bootstrapping.
            Return the results for every subsample rather than the average. Multiple
            values can be given as a list.
        bootstrap_samples : int
            When bootstrapping, sample <Int> times from the input data.
        input_degrees : bool
            Set `True` if input files are in degrees.
        model_data_degrees : bool
            Set `True` if the data given under `model_data` to is given in degrees.
        precision : int
            Sets the number of decimals in the output files. By default, 4 decimal.
        indent_size : int
            Sets the number of spaces used to indent the output document. By default, 0.
        kde_bandwidth : float
            This controls the bandwidth of the Gaussian kernel density estimator.
        grid_points : int
            When `model_type` == 'grid', this controls how many grid points
            are used to describe the probability density function.
        seed : int
            Set the random seed especially for bootstrapping.
    """
    # Input/Output Options
    input_files : typing.List[str] = field(default_factory=list)
    input_format : str = "auto"
    output_file : str = None
    output_format : str = "auto"

    # Conformational State Model Options
    model_type : str = "kde"
    model_load : str = None
    model_data : str = None
    model_dump : str = None

    # Subsampling Options
    window : typing.List[int] = field(default_factory=list)
    bootstrap : typing.List[int] = field(default_factory=list)
    window_series : typing.List[int] = field(default_factory=list)
    bootstrap_series : typing.List[int] = field(default_factory=list)
    bootstrap_samples : int = 500

    # Miscellaneous Options
    input_degrees : bool = False
    model_data_degrees : bool = False
    precision : int = 4
    indent_size : int = 0
    kde_bandwidth : float = .13
    grid_points : int = 10_000
    seed : int = None
    verbose : int = 0

    # References to the owning Constava object, not to be set by user
    _constava = None

    @set_logger_level
    @set_single_as_list
    def __setattr__(self, __attr, __value) -> None:
        """Custom function to set attributes, to catch certain special behaviours"""
        super().__setattr__(__attr, __value)
