"""constava.__main__ is the executable for the command line functionality of 
the tool."""

import sys, os
import argparse
from warnings import warn

from constava.constants import DEFAULT_KDE_PATH, DEFAULT_TRAINING_DATA_PATH
from constava.calc.calculator import ConfStateCalculator
from constava.calc.subsampling import SubsamplingBootstrap, SubsamplingWindow
from constava.calc.pdfestimators import KDEStatePdf, GridStatePdf
from constava.io.ensemblereader import EnsembleReader
from constava.io.resultswriter import ResultWriter



def parse_commandline_arguments(arguments):

    # Required for the description to consider return characters \n
    class CustomFormatter(argparse.RawDescriptionHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
        pass

    parser = argparse.ArgumentParser(description=(
            "This software is used to calculate Conformational State Variability "
            "(ConStaVa) from a protein structure ensemble. This is done by "
            "calculating the propensities for each conformational state for each "
            "residue in a protein ensemble. Then, ConStaVa is calculated from the "
            "change among Conformational States, inferred from trained Kernel "
            "Density Estimators (KDEs).\n\nBy default, this code retrains the "
            "Conformational States KDEs with the provided data set, as described in "
            "the associated publication. This will generate KDEs that are compatible "
            "with your current SciKit-learn version, and can be stored with the "
            "--kde-dump flag, to be later loaded with the --kde flag.\n\n"
            "Alternatively, the same pre-trained KDEs used for the publication can "
            "be used with the flag --use-publication-kdes. Due to SciKit-learn "
            "limitations, this flag must only be used with SciKit-learn version "
            "x.x.xx. If the  user wishes to train KDEs with a different set of "
            "dihedrals, this can be indicated with the flag --training-data. This "
            "must be provided in a json file, with the name of the conformational "
            "states as keys and a list of lists [[phi, psi], [phi, psi], ...] as values."),
        formatter_class=CustomFormatter)

    io_group = parser.add_argument_group("Input/Output Options")
    io_group.add_argument("-i", "--input-file", nargs="+", type=str, help="Input file with dihedral angles")
    io_group.add_argument("-o", "--output-file", type=str, help="Output file")
    io_group.add_argument("--input-format", choices=["auto", "xvg", "csv"], default="auto",
                        help="Format of input file")
    io_group.add_argument("--input-degrees", action="store_true",
                        help="Add this flag if input is provided in degrees (instead of radians)")
    io_group.add_argument("--output-format", choices=["auto", "csv", "json"], default="auto",
                        help="Format of output file")

    kde_group = parser.add_argument_group("KDE options")
    kde_group.add_argument("-k", "--kdes", metavar="<file.pkl>", help="Load KDEs from the given file")
    kde_group.add_argument("-d", "--kde-from-data", metavar="<data.json>", help=(
        "Fir KDEs from the given data. The data is provided in a json file, with "
        "the name of the conformational states as keys and a list of lists [[phi, "
        "psi], [phi, psi], ...] as values"))
    kde_group.add_argument("--kde-from-degrees", action="store_true", help=(
        "Add this flag if the data to fit the KDEs is provided in degrees (instead of radians)"))
    kde_group.add_argument("--dump-kdes", metavar="<file.pkl>", type=str, help=(
        "Dump the fitted KDEs as a pickled file, so that they can be reused "
        "later on using the --kdes flag."))
    kde_group.add_argument("--kde-bandwidth", metavar="<float>", type=float, default=.13)
    #kde_group.add_argument("--use-publication-kdes", type=str,
    #                       help="Load KDEs used in publication. This requires sklearn version x.x.xx")

    misc_group = parser.add_argument_group("Miscellaneous Options")
    misc_group.add_argument("--window", metavar="<int>", type=int, nargs='+', help="Do inference using a moving reading-frame of <int> consecutive samples.")
    misc_group.add_argument("--bootstrap", metavar="<int>", type=int, nargs='+', help="Do inference using <Int> samples obtained through bootstrapping. (By default a run with 3 and 25 is performed.)")
    misc_group.add_argument("--bootstrap-samples", metavar="<int>", type=int, default=500, help="If bootstrap, sample <Int> times from the input data (default: 500)")
    misc_group.add_argument("--seed", metavar="<int>", type=int, default=None, required=False, help="Set random seed for bootstrap sampling (default: None)")
    misc_group.add_argument("--quick", action="store_true", help="Use grid-interpolation instead of KDEs")
    misc_group.add_argument("--precision", type=int, default=4, help='Sets de number of decimals in the output files. By default, 4 decimal.')

    # Do actual parameter parsing
    args = parser.parse_args(arguments)

    # Some validity checks for Input/output
    if args.input_file is not None and args.output_file is None:
        raise argparse.ArgumentError("Missing argument: --output-file")
    if args.input_file is None and args.output_file is not None:
        raise argparse.ArgumentError("Missing argument: --input-file")
    if args.input_file is None and args.output_file is None and args.kde_from_data is None:
        raise argparse.ArgumentError((
            "Missing argument: You need to provide either input and output file "
            "or kde-training-data"))
    
    # Some validity checks for KDEs
    if args.kdes is not None and args.kde_from_data is not None:
        warn("Loading KDEs from provided file. Training data will be ignored")
    if args.kdes is None and args.kde_from_data is None:
        if os.path.isfile(DEFAULT_KDE_PATH):
            args.kdes = DEFAULT_KDE_PATH
        else:
            args.kde_from_data = DEFAULT_TRAINING_DATA_PATH
            args.dump_kdes = DEFAULT_KDE_PATH

    # Some validity checks for Misc
    if args.bootstrap is None and args.window is None:
        args.bootstrap = [3, 25]
    elif args.bootstrap is None and args.window is not None:
        args.bootstrap = []
    if args.window is None:
        args.window = []

    return args


def main():
    args = parse_commandline_arguments(sys.argv[1:])
    #print(args)

    # Read input files
    reader = EnsembleReader(filetype_str = args.input_format, 
                            degrees2radians = args.input_degrees)
    ensemble = reader.readFiles(*args.input_file)

    # Load/Fit KDEs
    PDFEstimator = GridStatePdf if  args.quick else KDEStatePdf
    if args.kdes is not None:
        pdfestimator = PDFEstimator.from_pickle(args.kdes)
    elif args.kde_from_data is not None:
        pdfestimator = PDFEstimator.from_fitting(
            args.kde_from_data,
            bandwidth = args.kde_bandwidth,
            degrees2radians = args.kde_from_degrees)
    else:
        # TODO raise no KDE Exception
        pass
    # IF dump_kdes, THEN save the KDEs
    if args.dump_kdes is not None:
        pdfestimator.dump_pickle(args.dump_kdes)

    # Load the calculation methods
    cscalc = ConfStateCalculator(pdfestimator)
    for window_size in args.window:
        cscalc.add_method(SubsamplingWindow(window_size))
    for sample_size in args.bootstrap:
        cscalc.add_method(SubsamplingBootstrap(sample_size, args.bootstrap_samples, seed=args.seed))

    # Do the inference
    results = cscalc.calculate(ensemble)

    # Write output
    writer = ResultWriter(args.output_format, args.precision)
    writer.writeToFile(results, args.output_file)
    
if __name__ == "__main__":
    sys.exit(main())
