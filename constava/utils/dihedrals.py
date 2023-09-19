"""constava.dihedrals is a stand-alone executable script that extracts the
phi/psi backbone dihedral angles from a conformational ensmeble. """

import argparse, os, sys
from typing import List, NamedTuple

import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
import numpy as np
import pandas as pd

def parse_arguments(cmdline_arguments: List[str]) -> NamedTuple:
    """Parses the command line arguments and does some minor sanity checking.
    
    Parameters:
    -----------
        cmdline_arguments : List[str]
            Parses the command line arguments provided as a list of strings, as
            returned by `sys.argv`
    
    Returns:
    --------
        args : NameSpace
            Object containing the parsed parameters.    
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--structure",  
        help="Structure file with atomic information: [pdb, gro, tpr]")
    parser.add_argument("-f", "--trajectory", nargs="+", 
        help="Trajectory file with coordinates: [pdb, gro, trr, xtc, crd, nc]")
    parser.add_argument("-o", "--output", default=None, required=False,
        help="(Optional) CSV file to write dihedral information to. (default: dihedrals.csv)")
    parser.add_argument("--selection", default="protein",
        help="(Optional) Selection for the dihedral calculation. (default: 'protein')")
    parser.add_argument("--precision", default=5, type=int,
        help="(Optional) Defines the number of decimals written for the dihedrals. (default: 5)")
    parser.add_argument("--degrees", action="store_true",
        help="(Optional) Results are written in degrees instead of radians. (default: radians)")
    args = parser.parse_args(cmdline_arguments)

    if args.output is None:
        args.output = "dihedrals.csv"
        if os.path.exists(args.output):
            raise FileExistsError((
                "Output file `{args.outfile}` exists already. Please, specify "
                "--outfile explicitely to overwrite or choose a different file name."
                ))
    return args

def calculate_dihedrals(structure: str, trajectory: List[str], selection: str = "protein", degree: bool = False):
    """Calculates the backbone dihedral angles and returns the result as a 
    DataFrame.
    
    Parameters:
    -----------
        structure : str
            Path to the structure File to use [*.gro, *.pdb, *.tpr, ...]
        trajectory : List[str]
            List of paths to trajectory files corresponding to the structure 
            file [*.trr, *.xtc, ...]
        selection : str
            Atom selection for which the dihedral angles should be calculated. 
            Corresponds to the MDAnalysis selection syntax. [default: `protein`]
        degrees : bool
            If True dihedrals are reported in degrees, else in radians. 
            [default: False]
    
    Returns:
    --------
        dihedrals : DataFrame
            DataFrame containing the backbone dihedral angles calculated along 
            the trajectory
    """
    # Load trajectory files
    u = mda.Universe(structure, trajectory)
    # Calculate dihedral angles.
    rama = Ramachandran(u.select_atoms(selection))
    rama.run()
    if degree:
        phicol, psicol = "Phi[deg]", "Psi[deg]"
        results = rama.results.angles
    else:
        phicol, psicol = "Phi[rad]", "Psi[rad]"
        results = np.radians(rama.results.angles)
    # Convert resulting numpy array into a comprehensive DataFrame
    dihedrals = pd.DataFrame(columns = ["#Frame", "ResIndex", "ResName", phicol, psicol])
    # Iterate over selection (C-alphas) used to define the dihedrals
    for i, atom in enumerate(rama.ag3):
        # Select results for given
        df = pd.DataFrame({
            "#Frame": rama.frames,
            "ResIndex": np.full(rama.frames.shape, atom.resid, dtype=int),
            "ResName": np.full(rama.frames.shape, atom.resname),
            phicol: results[:,i,0], 
            psicol: results[:,i,1]})
        dihedrals = pd.concat([dihedrals, df], ignore_index=True)
    return dihedrals

def main(cmdline_arguments: List[str]):
    """Main function executed when the script is run from the command line"""
    # Parse command line arguments
    args = parse_arguments(cmdline_arguments)
    # Calculate dihedrals
    dihedrals = calculate_dihedrals(
        args.structure, args.trajectory, args.selection, args.degrees)
    # Write results
    float2str = f"%.{args.precision}f" # Definition of float format in output
    dihedrals.to_csv(args.output, header=True, index=False, float_format=float2str)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))