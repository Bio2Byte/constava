import argparse, os, sys
import csv
from typing import List, NamedTuple

import MDAnalysis as mda
from MDAnalysis.analysis.dihedrals import Ramachandran
import numpy as np


def parse_arguments(arguments: List[str]) -> NamedTuple:
    ''' Parse command line arguments '''
    parser = argparse.ArgumentParser()
    parser.add_argument("-s", "--structure",  
        help=("Structure file with atomic information: [pdb, gro, tpr]"))
    parser.add_argument("-f", "--trajectory", nargs="+", 
        help=("Trajectory file with coordinates: [pdb, gro, trr, xtc, crd, nc]"))
    parser.add_argument("-o", "--outfile", default=None, required=False,
        help=("(Optional) CSV file to write dihedral information to. Default: dihedrals.csv"))
    parser.add_argument("--selection", default="protein",
        help=("(Optional) Selection for the dihedral calculation. Default: 'protein'"))
    parser.add_argument("--precision", default=5, type=int,
        help=("(Optional) Defines the number of decimals written for the dihedrals "
              "(in radians). Default: 5"))
    args = parser.parse_args(arguments)

    if args.outfile is None:
        args.outfile = "dihedrals.csv"
        if os.path.exists(args.outfile):
            raise FileExistsError((
                "Output file `{args.outfile}` exists already. Please, specify "
                "--outfile explicitely to overwrite or choose a different file name."
                ))
    return args

def main():
    # Parse command line arguments
    args = parse_arguments(sys.argv[1:])
    # Load conformational ensemble
    u = mda.Universe(args.structure, args.trajectory)
    # Calculate dihedral angles.
    rama = Ramachandran(u.select_atoms(args.selection))
    rama.run()
    dihedrals = np.radians(rama.results.angles)

    # Write the output
    float2str = lambda x: "{0:.{1}f}".format(x, args.precision)
    with open(args.outfile, "w") as ohandle:
        writer = csv.writer(ohandle)
        writer.writerow(["#Frame", "ResIndex", "ResName", "Phi[rad]", "Psi[rad]"])
        # Iterate over C-alphas defining the dihedrals
        for i, atom in enumerate(rama.ag3):
            resname = atom.resname
            resid = atom.resid
            for j, frame in enumerate(rama.frames):
                phi, psi = dihedrals[j,i,:]
                writer.writerow([frame, resid, resname, float2str(phi), float2str(psi)])

if __name__ == "__main__":
    sys.exit(main())