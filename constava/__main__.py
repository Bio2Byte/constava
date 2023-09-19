"""constava.__main__ is the executable for the command line functionality of 
the tool."""
import os, sys
import argparse
import textwrap as tw
from constava import Constava, ConstavaParameters, __version__
from constava.utils.dihedrals import calculate_dihedrals


def parse_parameters(cmdline_arguments):
    """Parse command line arguments and return them as ConstavaParameters object
    """
    parser = argparse.ArgumentParser(description=tw.dedent(
        """\
        Constava analyzes conformational ensembles calculating conformational state 
        propensities and conformational state variability. The conformational state 
        propensities indicate the likelihood of a residue residing in a given 
        conformational state, while the conformational state variability is a measure 
        of the residues ability to transiton between conformational states.
        
        Each conformational state is a statistical model of based on the backbone 
        dihedrals (phi, psi). The default models were derived from an analysis of NMR
        ensembles and chemical shifts. To analyze a conformational ensemble, the phi- 
        and psi-angles for each conformational state in the ensemble need to be 
        provided. 
        
        The `constava dihedrals` submodule provides a simple way to extract backbone 
        dihedral angles from MD simulations or PDB ensembles. For more information
        run: `constava dihedrals -h`. Alternatively, the backbone dihedrals may be
        extracted with GROMACS' `gmx chi` module.

        The `constava analyze` submodule analyzes the provided backbone dihedral angles
        and infers the propensities for each residue to reside in a given 
        conformational state. For more information run: `constava analyze -h`.

        The `constava fit-model` can be used to train a custom probabilistic model of
        confromational states.  The default models were derived from an analysis of NMR
        ensembles and chemical shifts; they cover six conformational states:
            * Core Helix - Exclusively alpha-helical, low backbone dynamics
            * Surrounding Helix - Mostly alpha-helical, high backbone dynamics
            * Core Sheet - Exclusively beta-sheet, low backbone dynamics
            * Surrounding Sheet - Mostly extended conformation, high backbone dynamics
            * Turn - Mostly turn, high backbone dynamics
            * Other - Mostly coil, high backbone dynamics"""),
        formatter_class=argparse.RawTextHelpFormatter, add_help=False)
    
    # General flags for the main software
    genOpt = parser.add_argument_group("Generic options")
    genOpt.add_argument("-h", "--help", action="help", help=tw.dedent(
        """\
        Show this help message and exit. For detailled 
        information on the subcommands, run: 
        `%(prog)s SUBCOMMAND -h`"""))
    genOpt.add_argument("--version", action="version", version=f"%(prog)s {__version__}", 
        help="Show the program's version number and exit")

    subparsers = parser.add_subparsers(title="Subcommands", 
        dest="subcommand", required=True, help=tw.dedent(
        """\
        fit-model
            Fit a custom conformational state model
        analyze
            Analyze a conformational ensemble using a 
            conformational state model
        dihedrals
            Obtain the phi/psi backbone dihedral angles from a MD 
            simulation"""))

    # ======================
    #  Subparser: fit-model
    # ======================
    parser_fit_model = subparsers.add_parser("fit-model", description=tw.dedent(
        """\
        The `constava fit-model` submodule is used to generate the probabilistic
        conformational state models used in the analysis. By default, when running
        `constava analyze` these models are generated on-the-fly. In selected cases 
        generating a model beforehand and loading it can be useful, though.

        We provide two model types. kde-Models are the default. They are fast to fit
        but may be slow in the inference in large conformational ensembles (e.g., 
        long-timescale MD simulations). The idea of grid-Models is, to replace
        the continuous probability density function of the kde-Model by a fixed set
        of grid-points. The PDF for any sample is then estimated by linear 
        interpolation between the nearest grid points. This is slightly less
        accurate then the kde-Model but speeds up inference significantly."""),
        formatter_class=argparse.RawTextHelpFormatter)
    
    fitIO = parser_fit_model.add_argument_group("Input and output options")
    fitIO.add_argument("-i", "--input", type=str, metavar="<file.json>", help=tw.dedent(
        """\
        The data to which the new conformational state models will
        be fitted. It should be provided as a JSON file. The 
        top-most key should indicate the names of the 
        conformational states. On the level below, lists of phi-/
        psi pairs for each stat should be provided. If not provided 
        the default data from the publication will be used."""))
    fitIO.add_argument("-o", "--output", type=str, metavar="<file.pkl>", required=True, help=tw.dedent(
        """\
        Write the generated model to a pickled file, that can be
        loaded gain using `constava analyze --load-model`"""))
    
    fitMdl = parser_fit_model.add_argument_group("Conformational state model options")
    fitMdl.add_argument("--model-type", choices=["kde", "grid"], default="kde", help=tw.dedent(
        """\
        The probabilistic conformational state model used. The 
        default is `kde`. The alternative `grid` runs significantly
        faster while slightly sacrificing accuracy: {'kde', 'grid'}
        (default: 'kde')"""))
    fitMdl.add_argument("--kde-bandwidth", type=float, metavar="<float>", default=.13, help=tw.dedent(
        """\
        This flag controls the bandwidth of the Gaussian kernel 
        density estimator. (default: 0.13)"""))
    fitMdl.add_argument("--grid-points", type=int, metavar="<int>", default=10_000, help=tw.dedent(
        """\
        This flag controls how many grid points are used to 
        describe the probability density function. Only applies if
        `--model-type` is set to `grid`. (default: 10000)"""))
    
    fitMisc = parser_fit_model.add_argument_group("Miscellaneous options")
    fitMisc.add_argument("--degrees", action="store_true", help=tw.dedent(
        """\
        Set this flag, if dihedrals in `model-data` are in degrees 
        instead of radians."""))
    fitMisc.add_argument("-v", "--verbose", action="count", default=0, help=tw.dedent(
        """\
        Set verbosity level of screen output. Flag can be given 
        multiple times (up to 2) to gradually increase output to 
        debugging mode."""))
    
    # ====================
    #  Subparser: analyze
    # ====================
    parser_analyze = subparsers.add_parser("analyze", description=tw.dedent(
        """\
        The `constava analyze` submodule analyzes the provided backbone dihedral angles
        and infers the propensities for each residue to reside in a given 
        conformational state. 

        Each conformational state is a statistical model of based on the backbone 
        dihedrals (phi, psi). The default models were derived from an analysis of NMR
        ensembles and chemical shifts. To analyze a conformational ensemble, the phi- 
        and psi-angles for each conformational state in the ensemble need to be 
        provided. 
        
        As input data the backbone dihedral angles extracted from the conformational 
        ensemble need to be provided. Those can be generated using the 
        `constava dihedrals` submodule (`--input-format csv`) or GROMACS'
        `gmx chi` module (`--input-format xvg`)."""),
        formatter_class=argparse.RawTextHelpFormatter)
    
    anaIO = parser_analyze.add_argument_group("Input & output options")
    anaIO.add_argument("-i", "--input", nargs="+", type=str, metavar="<file.csv>", 
        help="Input file(s) that contain the dihedral angles.")
    anaIO.add_argument("--input-format", choices=["auto", "xvg", "csv"], default="auto", 
        help="Format of the input file: {'auto', 'csv', 'xvg'}")
    anaIO.add_argument("-o", "--output", type=str, metavar="<file.csv>",
        help="The file to write the results to.")
    anaIO.add_argument("--output-format", choices=["auto", "csv", "json", "tsv"], default="auto",
        help="Format of output file: {'csv', 'json', 'tsv'}. (default: 'auto')")

    anaMdl = parser_analyze.add_argument_group("Conformational state model options")
    anaMdl.add_argument("-m", "--load-model", type=str, metavar="<file.pkl>", help=tw.dedent(
        """\
        Load a conformational state model from the given pickled 
        file. If not provided, the default model will be used."""))
    
    anaSmpl = parser_analyze.add_argument_group("Subsampling options")
    anaSmpl.add_argument("--window", metavar="<int>", type=int, nargs='+', help=tw.dedent(
        """\
        Do inference using a moving reading-frame. Each reading 
        frame consists of <int> consecutive samples. Multiple 
        values can be provided."""))
    anaSmpl.add_argument("--window-series", metavar="<int>", type=int, nargs='+', help=tw.dedent(
        """\
        Do inference using a moving reading-frame. Each reading 
        frame consists of <int> consecutive samples. Return the 
        results for every window rather than the average. This can
        result in very large output files. Multiple values can be 
        provided."""))
    anaSmpl.add_argument("--bootstrap", metavar="<int>", type=int, nargs='+',  help=tw.dedent(
        """\
        Do inference using <Int> samples obtained through 
        bootstrapping. Multiple values can be provided."""))
    anaSmpl.add_argument("--bootstrap-samples", metavar="<int>", type=int, default=500, help=tw.dedent(
        """\
        When bootstrapping, sample <Int> times from the input data.
        (default: 500)"""))
    
    anaMisc = parser_analyze.add_argument_group("Miscellaneous options")
    anaMisc.add_argument("--degrees", action="store_true", help=tw.dedent(
        """\
        Set this flag, if dihedrals in the input files are in 
        degrees."""))
    anaMisc.add_argument("--precision", type=int, default=4, 
        help="Sets the number of decimals in the output files.")
    anaMisc.add_argument("--seed", metavar="<int>", type=int, default=None, 
        required=False, help="Set random seed for bootstrap sampling")
    anaMisc.add_argument("-v", "--verbose", action="count", default=0, help=tw.dedent(
        """\
        Set verbosity level of screen output. Flag can be given 
        multiple times (up to 2) to gradually increase output to 
        debugging mode."""))

    # =====================
    # Subparser: dihedrals
    # =====================
    parser_dihedrals = subparsers.add_parser("dihedrals", description=tw.dedent(
        """\
        The `constava dihedrals` submodule is used to extract the backbone dihedrals
        needed for the analysis from confromational ensembles. By default the results
        are written out in radians as this is the preferred format for 
        `constava analyze`.
        
        Note: For the first and last residue in a protein only one backbone dihedral
        can be extracted. Thus, those residues are omitted by default."""),
        formatter_class=argparse.RawTextHelpFormatter)

    dihIO = parser_dihedrals.add_argument_group("Input & output options")
    dihIO.add_argument("-s", "--structure", metavar="<file.pdb>", 
        help="Structure file with atomic information: [pdb, gro, tpr]")
    dihIO.add_argument("-f", "--trajectory", nargs="+", metavar="<file.xtc>", 
        help="Trajectory file with coordinates: [pdb, gro, trr, xtc, crd, nc]")
    dihIO.add_argument("-o", "--output", default=None, required=False,
        help="CSV file to write dihedral information to. (default: dihedrals.csv)")
    
    dihMisc = parser_dihedrals.add_argument_group("Input & output options")
    dihMisc.add_argument("--selection", default="protein",
        help="Selection for the dihedral calculation. (default: 'protein')")
    dihMisc.add_argument("--precision", default=5, type=int,
        help="Defines the number of decimals written for the dihedrals. (default: 5)")
    dihMisc.add_argument("--degrees", action="store_true",
        help="If set results are written in degrees instead of radians.")
    dihMisc.add_argument("-O", "--overwrite", action="store_true",
        help="If set any previously generated output will be overwritten.")

    # Parse command line arguments
    return parser.parse_args(cmdline_arguments)

def run_fit_model(args):
    """Run fit-model subcommand when invoked from command line"""
    # Initialze and run Constava
    cva = Constava(ConstavaParameters(verbose=args.verbose))
    csmodel = cva.fit_csmodel(
        model_type = args.model_type,
        model_data = args.input,
        kde_bandwidth = args.kde_bandwidth,
        grid_points = args.grid_points,
        model_data_degrees = args.degrees)
    # Write the fitted model out as a pickle
    csmodel.dump_pickle(args.output)

def run_analyze(args):
    """Run analyze subcommand when invoked from command line."""
    # Convert command line arguments to ConstavaParameters
    params = ConstavaParameters(verbose=args.verbose)
    params.input_files = args.input
    params.input_format = args.input_format
    params.output_file = args.output
    params.output_format = args.output_format
    params.model_load = args.load_model
    params.window = args.window
    params.window_series = args.window_series
    params.bootstrap = args.bootstrap
    params.bootstrap_samples = args.bootstrap_samples
    params.input_degrees = args.degrees
    params.precision = args.precision
    params.seed = args.seed
    # Initialze and run Constava
    cva = Constava(params)
    cva.run()

def run_dihedrals(args):
    """Run analyze subcommand when invoked from command line."""
    # Set output to default value if needed
    args.output = args.output or "dihedrals.csv"
    if not args.overwrite and os.path.exists(args.output):
        raise FileExistsError(f"Cannot overwrite existing file: {args.output}")
    # Calculate dihedrals
    dihedrals = calculate_dihedrals(args.structure, args.trajectory, 
                                    args.selection, args.degrees)
    # Write results
    float2str = f"%.{args.precision}f" # Definition of float format in output
    dihedrals.to_csv(args.output, header=True, index=False, float_format=float2str)

def main():
    """main function executed when running script in command line mode"""
    # Parse command line parameters
    args = parse_parameters(sys.argv[1:])
    if args.subcommand == "fit-model":
        run_fit_model(args)
    elif args.subcommand == "analyze":
        run_analyze(args)
    elif args.subcommand == "dihedrals":
        run_dihedrals(args)

if __name__ == "__main__":
    sys.exit(main())