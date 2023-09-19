"""Stand-alone functions that are used in multiple submodules"""

from warnings import warn
import numpy as np

class DihedralRangeError(ValueError):
    """Raised if any dihedral angles are not correctly in radians"""
    pass

class DihedralRangeWarning(UserWarning):
    """Raised on the suspicion that  dihedral angles are not correctly in radians"""
    pass

def check_dihedral_range(arr: np.ndarray):
    """Helper method that checks if the dihedral angles are correctly in 
    radians.
    
    Parameters:
    -----------
        arr : Array[N,2]
            An array of N (phi, psi) pairs.

    Raises:
    -------
        DihedralRangeError
            If any dihedrals fall outside the range [-pi, pi]

        DihedralRangeWarning
            If all dihedrals fall in the range of [-(pi*pi/180), (pi*pi/180)], 
            as this suggests that angles were converted to radians twice.
    """
    vmin, vmax = np.min(arr), np.max(arr)
    if vmin < -np.pi or vmax > np.pi:
        raise DihedralRangeError(f"Dihedrals outside the range [-pi, pi] detected: [{vmin:.3f}, {vmax:.3f}]")
    elif vmin >= np.radians(np.pi) and vmax <= np.radians(np.pi):
        warn(("Provided dihedrals a very small: [{vmin:.3f}, {vmax:.3f}]. "
            "Please check that convertion to radians was only applied once."), 
            DihedralRangeWarning())
