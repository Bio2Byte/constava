# import sys
# import argparse
import os
import glob
from typing import List, Union, Optional, Literal

from constava.constants import DEFAULT_KDE_PATH, DEFAULT_TRAINING_DATA_PATH
from constava.calc.calculator import ConfStateCalculator
from constava.calc.subsampling import SubsamplingBootstrap, SubsamplingWindow
from constava.calc.pdfestimators import KDEStatePdf, GridStatePdf
from constava.io.ensemblereader import EnsembleReader
from constava.io.resultswriter import ResultWriter



class ConStaVa:
	"""
	Initializes all the variables in required for the fitting and calculation of the Conformational State Variability.
	It also stores all relevant information to train, save and load Kernel Density Estimators (KDE) as well as to
	calculate and save Conformational State Variability for different windows and bootstrap parameters.

	Parameters:
		None
	"""

	def __init__(self):
		"""
		Constructor for the ConStaVa class. Initializes all variables that the different methods can require.
		Does not return anything.

		Parameters:
			None

		Returns:
			None
		"""
		self.bootstrap_seed = None
		self.results = None
		self.pdfestimator = None
		self.ensemble = None
		self.kdes = None
		self.kde_from_data = None
		self.kde_from_degrees = False
		self.dump_kdes = None
		self.bootstrap = None
		self.bootstrap_samples = 500
		self.window = None
		self.input_file = None
		self.input_format = 'auto'
		self.output_file = None
		self.output_format = 'auto'
		self.input_degrees = False
		self.quick = False
		self.precision = 4
		self.kde_bandwidth = 0.13

	def train_kde(self, training_data: Optional[str] = None, kde_degrees: bool = False):
		"""
		It trains KDEs for all conformational states according given a json file with the names of conformational
		states as keys and their corresponding phi and psi angles in a list of lists as values {"conf0":[[phi, psi],
		[phi, psi], ...], "conf1": [...], ...}

		Parameters:
			training_data (str, optional): path to the json file which contains the dihedrals for KDE fitting. If None,
			a default set of angles will be used to fit the KDEs.

			kde_degrees (bool, optional): Indicate with True if the dihedrals in training_data are in degrees rather
			than radians. Leave empty or use False if the training data is in radians.

		"""
		self.kde_from_data = training_data
		self.kde_from_degrees = kde_degrees

		if self.kde_from_data is None:
			self.kde_from_data = DEFAULT_TRAINING_DATA_PATH

		self._fit_kde()

	def _fit_kde(self):
		"""
		Performs the fitting of the KDEs once the right parameters have been fetched in self.train_kde().
		Does not return anything.

		Parameters:
			None

		Returns:
			None
		"""
		PDFEstimator = GridStatePdf if self.quick else KDEStatePdf 
		self.pdfestimator = PDFEstimator.from_fitting(
			self.kde_from_data,
			bandwidth=self.kde_bandwidth,
			degrees2radians=self.input_degrees)

	def load_kde(self, kdes_path: Optional[str] = None):
		"""
		Loads pre-fitted KDEs which are stored in a path. Does not return anything.

		Parameters:
			kdes_path (Optional[str], default is None): Path where the KDEs file is (pickle file). Will use default
			path if not provided.

		Returns:
			None
		"""
		self.kdes = kdes_path  # Needs to be pickle (I think)
		PDFEstimator = GridStatePdf if self.quick else KDEStatePdf
		if self.kdes is None:
			if os.path.isfile(DEFAULT_KDE_PATH):
				self.kdes = DEFAULT_KDE_PATH
		self.pdfestimator = PDFEstimator.from_pickle(self.kdes)

	def save_kde(self, dump_kdes: Optional[str] = None):
		"""
		Saves the trained KDEs, by default as "constava_default_kdes.pkl". Does not return anything.

		Parameters:
			dump_kdes (Optional[str], default is None): Save path for the KDE file. Must be a .pkl file.

		Returns:
			None
		"""

		if dump_kdes is not None:
			self.dump_kdes = dump_kdes
		self.pdfestimator.dump_pickle(self.dump_kdes)

	# Read input files
	def read_input_files(self, infile: Union[str, List[str]], degrees: bool = False,
	                     input_format: Literal['auto', 'csv', 'json'] = 'auto'):
		"""
		Reads the input files to calculate ConStaVa. Does not return anything.

		Parameters:
			infile (Union[str, List[str]]): A string or list of strings where each string is a path to a protein structure ensemble file.
			If a string is provided, it should represent the path in bash syntax (e.g. path/to/files/*.xvg).

			degrees (bool, default is False): True if the dihedrals of the ensembles are in degrees rather than radians.

			input_format (Literal['auto', 'csv', 'json'], default is 'auto'): Specifies the format of the input file.

		Returns:
			None
		"""

		if type(infile) == list:
			self.input_file = infile
		elif type(infile) == str:
			self.input_file = glob.glob(infile)
		self.input_degrees = degrees
		self.input_format = input_format

		reader = EnsembleReader(filetype_str=self.input_format,
		                        degrees2radians=self.input_degrees)
		self.ensemble = reader.readFiles(*self.input_file)

	def calculate_results(self, window: Optional[List[int]] = None, bootstrap: Optional[List[int]] = None,
	                      bootstrap_samples: int = 500, bootstrap_seed: Optional[int] = None):
		"""
		Calculates ConStaVa given that the KDEs have been trained or loaded and the protein ensembles are also loaded.
		Does not return anything.

		Parameters:
			window (Optional[List[int]], default is None): Specifies window sizes for calculating the ConStaVa. Each
			element in the list represents a different window size.

			bootstrap (Optional[List[int]], default is None): Specifies the bootstrap sample sizes for calculating the ConStaVa. Each
			element in the list represents a different bootstrap sample size. (By default a run with 3 and 25 is
			performed if no window or bootstrap is provided.)

			bootstrap_samples (int, default is 500): The number of bootstrap samples to generate for the calculation of the ConStaVa.

			bootstrap_seed (Optional[int], default is None): The seed for random number generator in bootstrap sampling.

		Returns:
			None
		"""

		self.window = window
		self.bootstrap = bootstrap
		self.bootstrap_samples = bootstrap_samples
		self.bootstrap_seed = bootstrap_seed

		# Load the calculation methods
		cscalc = ConfStateCalculator(self.pdfestimator)

		if self.window is None and self.bootstrap is None:
			self.bootstrap = [3, 25]

		if self.window is not None:
			for window_size in self.window:
				cscalc.add_method(SubsamplingWindow(window_size))

		if self.bootstrap is not None:
			for sample_size in self.bootstrap:
				cscalc.add_method(SubsamplingBootstrap(sample_size, self.bootstrap_samples, seed=self.bootstrap_seed))

		# Do the inference
		self.results = cscalc.calculate(self.ensemble)

	def save_results(self, output_file: str, output_format: Literal['auto', 'csv', 'json'] = 'auto',
	                 precision: int = 4):
		"""
		Saves the results of the calculations. Does not return anything.

		Parameters:
			output_file (str): Path to output file. Format must be .csv or .json. CSV follows a standard dataframe format for easy data handling.

			output_format (Literal['auto', 'csv', 'json'], default is 'auto'): Specifies save format.

			precision (int, default is 4): Sets de number of decimals in the output files. By default, 4 decimal

		Returns:
			None
		"""

		self.output_file = output_file
		self.output_format = output_format
		self.precision = precision

		# Write output
		writer = ResultWriter(self.output_format, self.precision)
		writer.writeToFile(self.results, self.output_file)
