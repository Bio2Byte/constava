"""constava.constants contains general information for all modules"""
from importlib.metadata import version
import os

CONSTAVA_NAME = __name__.split(".")[0]
CONSTAVA_VERSION = version(CONSTAVA_NAME)
CONSTAVA_SOURCE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONSTAVA_DATA_DIR = os.path.join(CONSTAVA_SOURCE_DIR, "data")
DEFAULT_TRAINING_DATA_PATH = os.path.join(CONSTAVA_DATA_DIR, "constava_csdata.json")

aminoacids1to3 = dict(
	A="ALA", C="CYS", D="ASP", E="GLU", F="PHE",
	G="GLY", H="HIS", I="ILE", K="LYS", L="LEU",
	M="MET", N="ASN", P="PRO", Q="GLN", R="ARG",
	S="SER", T="THR", V="VAL", W="TRP", Y="TYR",
)

aminoacids3to1 = {j: i for i, j in aminoacids1to3.items()}
