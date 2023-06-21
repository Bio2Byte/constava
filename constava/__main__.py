import sys
import argparse

from constava.calculator import ConfStateCalculator
from constava.ensemblereader import EnsembleReader
from constava.methods import ConstavaBootstrap, ConstavaWindow
from constava.resultwriter import ResultWriter
from constava.pdfestimators import KDEStatePdf


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

    required_group = parser.add_argument_group("required arguments")
    required_group.add_argument("-i", "--infile", nargs="+", type=str, required=True,
                                help="Input file with dihedral angles")
    required_group.add_argument("-o", "--outprefix", type=str, required=True, help="Output file")

    format_group = parser.add_argument_group("file format options")
    format_group.add_argument("--input-format", choices=["auto", "xvg", "csv"], default="auto",
                              help="Format of input file")
    format_group.add_argument("--output-format", choices=["auto", "csv", "json"], default="auto",
                              help="Format of output file")

    kde_group = parser.add_argument_group("KDE options")
    kde_group.add_argument("-k", "--kdes", type=str, help="Load KDEs from the given file")
    kde_group.add_argument("-d", "--training-data", type=str, help=(
        "Train KDEs with given data. This must be provided in a json file, with "
        "the name of the conformational states as keys and a list of lists "
        "[[phi, psi], [phi, psi], ...] as values"))
    kde_group.add_argument("--dump-kde", type=str, help=(
        "If desired to dump trained KDEs, write filename for to save "
        "KDEs. File format must be .XXX"))
    kde_group.add_argument("--kde-bandwidth", type=float, default=.13)
    #kde_group.add_argument("--use-publication-kdes", type=str,
    #                       help="Load KDEs used in publication. This requires sklearn version x.x.xx")

    misc_group = parser.add_argument_group("miscellaneous options")
    misc_group.add_argument("--window", type=int, nargs='+',
                            help="Subsampling using moving reading-frame of size <Int>")
    misc_group.add_argument("--bootstrap", type=int, nargs='+', help="Subsampling using <Int> bootstrapped samples")
    misc_group.add_argument("--bootstrap-samples", type=int, default=500, help="If bootstrap, sample <Int> times")
    misc_group.add_argument("--quick", action="store_true", help="Use grid-interpolation instead of KDEs")
    misc_group.add_argument("--degrees", action="store_true",
                            help="Convert degrees to radians. This is REQUIRED if the data is provided in degrees")
    misc_group.add_argument("--precision", type=int, default=4)

    # Tweak the parameters so if nothing is provided, the bootstraps used in the publication are returned.
    parsed_arguments = parser.parse_args(arguments)
    if parsed_arguments.bootstrap is None and parsed_arguments.window is None:
        parsed_arguments.bootstrap = [3, 25]
    elif parsed_arguments.bootstrap is None and parsed_arguments.window is not None:
        parsed_arguments.bootstrap = []

    if parsed_arguments.window is None:
        parsed_arguments.window = []

    return parsed_arguments


def main():
    args = parse_commandline_arguments(sys.argv[1:])
    #print(args)

    # Read input files
    reader = EnsembleReader(filetype_str = args.input_format, 
                            degrees2radians = args.degrees)
    ensemble = reader.readFiles(*args.infile)

    # Load/Fit KDEs
    if args.kdes is not None:
        pdfestimator = KDEStatePdf.from_pickle(args.kdes)
    elif args.training_data is not None:
        pdfestimator = KDEStatePdf.from_fitting(
            args.training_data,
            bandwidth = args.kde_bandwidth,
            degrees2radians=args.degrees)
    else:
        # TODO raise no KDE Exception
        pass
    # IF dump_kde, THEN save the KDEs
    if args.dump_kde is not None:
        pdfestimator.dump_pickle(args.dump_kde)

    # Load the calculation methods
    cscalc = ConfStateCalculator(pdfestimator)
    for window_size in args.window:
        cscalc.add_method(ConstavaWindow(window_size))
    for sample_size in args.bootstrap:
        cscalc.add_method(ConstavaBootstrap(sample_size, args.bootstrap_samples))

    # Do the inference
    results = cscalc.calculate(ensemble)

    # Write output
    writer = ResultWriter(args.output_format, args.precision)
    for result in results:
        method_name = result.method.getShortName()
        output_file = f"{args.outprefix}_{method_name}.{args.output_format.lower()}"
        writer.writeToFile(result, output_file)
    
if __name__ == "__main__":
    sys.exit(main())