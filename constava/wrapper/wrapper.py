import os
from typing import List

from ..utils.constants import DEFAULT_TRAINING_DATA_PATH
from .logging import logging
from .params import ConstavaParameters
from ..io import ResultsWriter, EnsembleReader
from ..calc.calculator import ConfStateCalculator
from ..calc.subsampling import SubsamplingBootstrap, SubsamplingWindow, SubsamplingWindowSeries
from ..calc.csmodels import ConfStateModelABC, ConfStateModelKDE, ConfStateModelGrid

# The logger for the wrapper
logger = logging.getLogger("Constava")

def cache_csmodel(func):
    """Decorator for caching conformational state models"""
    def __inner(self: "Constava", *args, **kwargs):
        if func.__name__ == "load_csmodel":
            new_cshash = hash((
                args[0] if len(args) > 0 else kwargs.get("pickled_csmodel", None),))
        elif func.__name__ == "fit_csmodel":
            new_cshash = hash((
                args[0] if len(args) > 0 else kwargs.get("model_type", "kde"),
                args[1] if len(args) > 1 else kwargs.get("model_data", None),
                args[2] if len(args) > 2 else kwargs.get("kde_bandwidth", .13),
                args[3] if len(args) > 3 else kwargs.get("grid_points", 10_000),
                args[4] if len(args) > 4 else kwargs.get("model_data_degrees", False),))
        else:
            raise TypeError("Decorator `cache_csmodel` only works with `load_csmodel` and `fit_csmodel`")
        
        if self._cshash == new_cshash:
            logger.info("No change to model parameters. Using preloaded model.")
        else:
            csmodel = func(self, *args, **kwargs)
            self._csmodel, self._cshash = csmodel, new_cshash
        return self._csmodel
    return __inner


class Constava:
    """Interface class for all functionalities of Constava.

    Methods:
    --------
        set_param(parameter, value)
            Sets a parameters to a given value.
        get_param(parameter) -> value
            Returns the current value of the given parameter
        show_params() -> str
            Returns the current set of parameters as a string
        run()
            Runs Constava with the current parameters
    """

    def __init__(self, parameters: ConstavaParameters = None, **kwargs):
        """Initializes the python interface for Constava. Parameters can be 
        provided as a ConstavaParameters class
        
        Parameters:
        -----------
            parameters : ConstavaParameters
                ConstavaParameters object ontaining all parameters (if provided
                kwargs will be ignored)
            **kwargs :
                To only set individual parameters, those parameters can be 
                provided as keyword arguments. For all other parameters
                default values are used. For a full list of available settings
                and their defaults, check: `help(ConstavaParameters)`
        """
        logger.info("Constava: Initializing python interface...")
        if parameters is None:
            parameters = ConstavaParameters(**kwargs)
        parameters._constava = self
        self.parameters = parameters
        self.results = None
        self._csmodel = None    # Preloaded conformational state models
        self._cshash = None     # Hashed parameters of the models

    def get_param(self, parameter: str):
        """Returns the current value of the given parameter"""
        return getattr(self.parameters, parameter)
        
    def set_param(self, parameter: str, value):
        """Sets a parameter to a given value"""
        logger.info(f"Setting `{parameter} = {value}`")
        setattr(self.parameters, parameter, value)
        logger.debug(f"New parameters: {self.show_params()}")

    def unset_param(self, parameter: str):
        """Sets a parameter to None"""
        setattr(self.parameters, parameter, None)

    def show_params(self) -> str:
        """Returns a string with all currently set parameters"""
        return repr(self.parameters)
    
    def run(self) -> None:
        """Calculate conformational state variabilities and conformational
        state variabilites with the given parameters."""
        # Reset results
        self.results = None

        # Initialize an reader for input file(s)
        reader = self.initialize_reader(
                format = self.get_param("input_format"), 
                in_degrees = self.get_param("input_degrees"))
        
        # Initialize writer for results
        writer = self.initialize_writer(
                outfile = self.get_param("output_file"),
                format = self.get_param("output_format"),
                float_precision = self.get_param("precision"))
        
        # Fit or load a conformational state model
        if os.path.isfile(self.get_param("model_load") or ""):
            csmodel = self.load_csmodel(pickled_csmodel = self.get_param("model_load"))
        else:
            csmodel = self.fit_csmodel(
                model_type = self.get_param("model_type"),
                model_data = self.get_param("model_data"), 
                kde_bandwidth = self.get_param("kde_bandwidth"),
                grid_points = self.get_param("grid_points"),
                model_data_degrees = self.get_param("model_data_degrees"))

        # Initialize a calculator (logged inside function)
        calculator = self.initialize_calculator(
                csmodel = csmodel,
                window = self.get_param("window"),
                window_series = self.get_param("window_series"),
                bootstrap = self.get_param("bootstrap"),
                bootstrap_samples = self.get_param("bootstrap_samples"),
                bootstrap_seed  = self.get_param("seed"))

        # Read input files
        input_files = self.get_param("input_files")
        logger.info(f"Reading dihedrals from {len(input_files)} files...")
        logger.debug("\n\t*  ".join(["... input file list:", *input_files]))
        ensemble = reader.readFiles(*input_files)
        
        # Do the inference
        logger.info("Starting inference...")
        self.results = calculator.calculate(ensemble)
        
        # Write results
        logger.info(f"Writing results to file: {writer.filename}")
        writer.write_results(self.results)
    
    def initialize_reader(self, format: str = "auto", in_degrees: bool = False) -> EnsembleReader:
        """Initializes an EnsembleReader.
        
        Parameters:
        -----------
            format: str
                File format of the files to be read. {'auto', 'csv', 'xvg'}
            in_degrees: bool
                Set `True` if input files are in degrees,
        
        Returns:
        --------
            reader
                An EnsembleReader object.
        """
        logger.info("Initializing reader for input file(s)...")
        logger.debug(f"... setting reader parameters: {format=}, {in_degrees=}")
        reader = EnsembleReader(filetype_str=format, degrees2radians=in_degrees)
        return reader
    
    def initialize_writer(self, outfile, format: str = "auto", float_precision: int = 4) -> ResultsWriter:
        """Initializes a ResultsWriter.
        
        Parameters:
        -----------
            outfile: str
                The file path to write results to. File extension is used to infer
                the output format, if not provided explicitly.
            format: str
                Format in which the results should be written out. {'auto', 'csv', 'json'}
            float_precision: int
                Sets de number of decimals in the output files. By default, 4 decimal.
        
        Returns:
        --------
            writer
                A ResultsWriter object.
        """
        if outfile is None:
            return None
        logger.info("Initializing writer for results...")
        logger.debug(f"... setting writer parameters: {outfile=}, {format=}, {float_precision=}")
        writer = ResultsWriter(outfile, format=format, float_precision=float_precision)
        return writer
    
    @cache_csmodel
    def fit_csmodel(self, model_type: str = "kde", model_data: str = None, 
            kde_bandwidth: float = .13, grid_points: int = 10_000,
            model_data_degrees: bool = False) -> ConfStateModelABC:
        """Fits a conformational state model to the provided data.
        
        Parameters:
        -----------
            model_type : str
                The probabilistic model used: {'kde', 'grid'}. (default: 'kde')
            model_data : str
                File in JSON format with the fitting data for the conformational 
                state models. If not provided, the default data from the 
                publication is used.
            kde_bandwidth : float
                This controls the bandwidth of the Gaussian kernel density 
                estimator. (default: 0.13)
            grid_points : int
                When `model_type` == 'grid', this controls how many grid points
                are used to describe the probability density function (default: 10000).
            model_data_degrees : bool
                Set `True` if the data given under `model_data` to is given in 
                degrees. (default: False)

        Returns:
        --------
            csmodel : ConfStateModelABC
                Probabilistic model describing the conformational states
        """
        PdfModel = {"kde": ConfStateModelKDE, "grid": ConfStateModelGrid}[model_type]
        model_data = (model_data or DEFAULT_TRAINING_DATA_PATH)
        logger.info(f"Fitting model to data in: {model_data}")
        csmodel = PdfModel.from_fitting(
            model_data, 
            in_degrees = model_data_degrees,
            bandwidth = kde_bandwidth,
            grid_points = grid_points)
        logger.info(f"... model fitted: {csmodel}")
        return csmodel

    @cache_csmodel
    def load_csmodel(self, pickled_csmodel: str) -> ConfStateModelABC:
        """Load a previously fitted conformational state model from a pickled
        file.
        
        Parameters:
        -----------
            pickled_csmodel : str
                Path to a pickled file from which the model is to be loaded.

        Returns:
        --------
            csmodel : ConfStateModelABC
                Probabilistic model describing the conformational states
        """
        logger.info(f"Loading conformational state models from file: {pickled_csmodel}")
        csmodel = ConfStateModelABC.from_pickle(pickled_csmodel)
        logger.info(f"... model loaded: {csmodel}")
        return csmodel

    def initialize_calculator(self, csmodel: ConfStateModelABC = None, 
            window: List[int] = None, window_series: List[int] = None, 
            bootstrap: List[int] = None, bootstrap_samples: int = 500, 
            bootstrap_seed: int = None) -> ConfStateCalculator:
        """Initializes a ConfStateCalculator.

        Parameters:
        -----------
            csmodel : ConfStateModel
                A conformational state model (as returned by the fit_csmodel 
                method). If None, the standard model from the publication will
                be used.
            window : List[int]
                Subsampling using a moving reading-frame of size <int>. Multiple
                values can be given as a list.
            window_series : List[int]
                Subsampling using a moving reading-frame of size <int>. Returns 
                the results for every window rather than the average. Multiple
                values can be given as a list.
            bootstrap : List[int]
                Subsampling using by bootstrapping <int> datapoints. Multiple
                values can be given as a list.
            bootstrap_samples : int
                When bootstrapping, sample <int> times from the input data.
            bootstrap_seed : int
                Random seed used when bootstrapping.

        Returns:
        --------
            calculator: ConfStateCalculator
                A ConfStateCalculator object
        """
        # Quickly generate a csmodel if not done before
        if csmodel is None:
            csmodel = self.fit_csmodel()
        # Initialize the calculator
        logger.info(f"Initializing calculator with {csmodel}...")
        calculator = ConfStateCalculator(csmodel)
        # Add subsampling methods to calculator
        for window_size in (window or []):
            new_method = SubsamplingWindow(window_size)
            logger.info(f"... adding subsampling method: {new_method.getShortName()}")
            calculator.add_method(new_method)
        for sample_size in (bootstrap or []):
            new_method = SubsamplingBootstrap(sample_size, bootstrap_samples, seed=bootstrap_seed)
            logger.info(f"... adding subsampling method: {new_method.getShortName()}")
            calculator.add_method(new_method)
        for window_size in (window_series or []):
            new_method = SubsamplingWindowSeries(window_size)
            logger.info(f"... adding subsampling method: {new_method.getShortName()}")
            calculator.add_method(new_method)
        return calculator
