# import sys
# import argparse
import os
import glob

from constava.calculator import ConfStateCalculator
from constava.ensemblereader import EnsembleReader
from constava.methods import ConstavaBootstrap, ConstavaWindow
from constava.resultswriter import ResultWriter
from constava.pdfestimators import KDEStatePdf


class ConStaVa:
	def __init__(self):
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
		# self.quick = False
		self.precision = 4
		self.kde_bandwidth = 0.13

	def train_kde(self, training_data=None, kde_degrees=False):
		self.kde_from_data = training_data
		self.kde_from_degrees = kde_degrees
		self.fit_kde()

	# Read input files
	def read_input_files(self, infile=None, degrees=False, input_format='auto'):
		# self.infile = [os.path.join(infile, the_file) for the_file in os.listdir(infile)]
		self.input_file = glob.glob(infile)
		self.input_degrees = degrees
		self.input_format = input_format

		reader = EnsembleReader(filetype_str=self.input_format,
		                        degrees2radians=self.input_degrees)
		self.ensemble = reader.readFiles(*self.input_file)

	def load_kde(self):
		self.pdfestimator = KDEStatePdf.from_pickle(self.kdes)

	def fit_kde(self):
		self.pdfestimator = KDEStatePdf.from_fitting(
			self.kde_from_data,
			bandwidth=self.kde_bandwidth,
			degrees2radians=self.input_degrees)

	def save_kde(self):
		self.pdfestimator.dump_pickle(self.dump_kdes)

	def calculate_results(self):
		# Load the calculation methods
		cscalc = ConfStateCalculator(self.pdfestimator)

		if self.window is None and self.bootstrap is None:
			self.bootstrap = [3, 5]

		if self.window is not None:
			for window_size in self.window:
				cscalc.add_method(ConstavaWindow(window_size))

		if self.bootstrap is not None:
			for sample_size in self.bootstrap:
				cscalc.add_method(ConstavaBootstrap(sample_size, self.bootstrap_samples))

		# Do the inference
		self.results = cscalc.calculate(self.ensemble)

	def save_results(self, output_file, output_format='auto', precision=4):
		self.output_file = output_file
		self.output_format = output_format
		self.precision = precision

		# Write output
		writer = ResultWriter(self.output_format, self.precision)
		writer.writeToFile(self.results, self.output_file)
