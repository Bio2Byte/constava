from dataclasses import dataclass, field
from typing import List, ClassVar
from .logging import logging

logger = logging.getLogger("Constava")


@dataclass
class ConstavaParameters:
    """The parameters that govern the function of Constava
    
    Parameters:
    -----------
        input_files : List[str]
            Input file(s) that contain the dihedral angles.
        input_format : str
            Format of the input file: {'auto', 'csv', 'xvg'}
        output_file : str
            The file to write the output to.
        output_format : str
            Format of output file: {'auto', 'csv', 'json', 'tsv'}

        model_type : str
            The probabilistic conformational state model used. Default is `kde`.
            The alternative `grid` runs significantly faster while slightly 
            sacrificing accuracy: {'kde', 'grid'}
        model_load : str
            Load a conformational state model from the given pickled file.
        model_data : str
            Fit conformational state models to data provided in the given file.
        model_dump : str
            Write the generated model to a pickled file, that can be loaded 
            again using `model_load`.

        window : List[int]
            Do inference using a moving reading-frame of <int> consecutive samples.
            Multiple values can be given as a list.
        window_series : List[int]
            Do inference using a moving reading-frame of <int> consecutive samples.
            Return the results for every window rather than the average. Multiple 
            values can be given as a list.
        bootstrap : List[int]
            Do inference using <Int> samples obtained through bootstrapping.
            Multiple values can be given as a list.
        bootstrap_samples : int
            When bootstrapping, sample <Int> times from the input data.

        input_degrees : bool
            Set `True` if input files are in degrees.
        model_data_degrees : bool
            Set `True` if the data given under `model_data` to is given in degrees.
        precision : int
            Sets the number of decimals in the output files. By default, 4 decimal.
        kde_bandwidth : float
            This controls the bandwidth of the Gaussian kernel density estimator.
        grid_points : int
            When `model_type` == 'grid', this controls how many grid points
            are used to describe the probability density function.
        seed : int
            Set the random seed especially for bootstrapping.
    """
    # Input/Output Options
    input_files : List[str] = None
    input_format : str = "auto"
    output_file : str = None
    output_format : str = "auto"

    # Conformational State Model Options
    model_type : str = "kde"
    model_load : str = None
    model_data : str = None
    model_dump : str = None

    # Subsampling Options
    window : List[int] = field(default_factory=list)
    bootstrap : List[int] = field(default_factory=list)
    window_series : List[int] = field(default_factory=list)
    bootstrap_samples : int = 500

    # Miscellaneous Options
    input_degrees : bool = False
    model_data_degrees : bool = False
    precision : int = 4
    kde_bandwidth : float = .13
    grid_points : int = 10_000
    seed : int = None
    verbose : int = 0

    # References to the owning Constava object, not to be set by user
    _constava = None

    def __setattr__(self, __attr, __value) -> None:
        """Custom function to set attributes, to catch certain special behaviours"""
        # Set the actual value
        super(ConstavaParameters, self).__setattr__(__attr, __value)
        # Some special behaviours for selected cases
        # ... "model_load" and "model_data" are mutually exclusive
        if __attr == "model_load" and __value is not None and self.model_data is not None:
            logger.warning("`model_load` takes precedence over `model_data`. The latter will be ignored.")
        # ... "model_load" and "model_data" are mutually exclusive
        elif __attr == "model_data" and __value is not None and __value and self.model_load is not None:
            logger.warning("`model_load` is not None. While `model_load` is set `model_data` will be ignored")
        # ... if verbosity is changed, apply to logger
        elif __attr == "verbose":
            if __value == 0:
                logger.setLevel(logging.WARNING)
            elif __value == 1:
                logger.setLevel(logging.INFO)
            else:
                logger.setLevel(logging.DEBUG)