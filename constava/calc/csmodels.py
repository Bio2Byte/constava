""" constava.csmodels contains various definitions how the conformational 
states can be represented as probabilistic models.
"""

import abc
import json
import pickle
import math
from typing import List, Tuple
import numpy as np
from scipy.interpolate import interpn
from sklearn.neighbors import KernelDensity
from ..utils.utils import check_dihedral_range


class ConfStateModelLoadingError(ValueError):
    pass


class ConfStateModelABC(metaclass=abc.ABCMeta):
    """AbstractBaseClass for the probabilistic model describing the 
    conformational states.
    """

    model_type : str = None

    def __init__(self, state_labels: List[str], **kwargs):
        """Initialize probabilistic model of conformational states. The first
        parameter are the labels of the conformational states, while further
        parameters describe the probabilistic model.
        """
        self.state_labels = tuple(state_labels)

    def __repr__(self):
        return f"<{self.model_type}Model states={self.state_labels}>"

    def get_labels(self):
        """Returns a copy of list of labels for the conformational state models"""
        return list(self.state_labels)

    @abc.abstractmethod
    def get_logpdf(self, data: np.ndarray) -> np.ndarray:
        """Inference of log-probability densities based on the probabilistic
        model.
        
        Parameters:
        -----------
            data : Array[N,2]
                N original observations of (phi, psi) angle pairs.
        
        Returns:
        --------
            logpdf : Array[M,N]
                Log-Probability densities for all M probabilistic conformational
                state models across the N original observations.
        """
        pass
    
    @abc.abstractclassmethod
    def from_fitting(cls, training_data_json: str, **kwargs):
        pass

    def dump_pickle(self, output_file: str):
        """Save the probabilistic as a pickle.
        
        Parameters:
        -----------
            output_file : str
                File path to store the model at.
        """
        with open(output_file, "wb") as fhandle:
            pickle.dump(self, fhandle)

    @classmethod
    def from_pickle(cls, pickled_file: str):
        """Load the probabilistic from a pickle.
        
        Parameters:
        -----------
            pickled_file : str
                File path to the pickled model.
        """
        with open(pickled_file, "rb") as fhandle:
            csmodel = pickle.load(fhandle)
        if not isinstance(csmodel, cls):
            raise ConfStateModelLoadingError((
                "Loaded conformational state models of wrong type: `{0}` "
                "(expected: `{1}`)").format(csmodel.model_type, cls.model_type))
        return csmodel


class ConfStateModelKDE(ConfStateModelABC):
    """Probabilistic model of conformational states based on a Gaussian kernel 
    density estimator.

    Attributes:
    -----------
        state_labels : Tuple
            List labels for the conformational state models
        state_kdes: Tuple
            List of the Gaussian kernel density estimators representing the 
            probabilistic models for the conformational states
    
    Methods:
    --------
        get_logpdf(data)
            Infer log-probability densities based on the probabilistic model.
        dump_pickle(output_file)
            Save the conformational state model as a pickle.

    Class methods:
    --------------
        from_pickle(pickled_file)
            Load conformational state models from pickled file
        from_fitting(json_file)
            Fit conformational state models to data provided in json file
    """
    model_type = "kde"

    def __init__(self, state_labels: List[str], state_kdes: List[KernelDensity]):
        """Initialize probabilistic model of conformational states based on 
        Gaussian kernel density estimator.
        
        Parameters:
        -----------
            state_labels : Tuple
                List labels for the conformational state models
            state_kdes: Tuple
                List of the Gaussian kernel density estimators representing the 
                probabilistic models for the conformational states
        """
        self.state_labels = tuple(state_labels)
        self.state_kdes = tuple(state_kdes)

    def get_logpdf(self, data: np.ndarray) -> np.ndarray:
        result = np.stack([
            kde.score_samples(data) for kde in self.state_kdes
        ])
        return result

    @classmethod
    def from_fitting(cls, training_data_json: str, *, in_degrees=False, bandwidth=.13, **_):
        """Generate the probabilistic models at runtime, by fitting the models
        to the provided training data. Training data must be a json.
        
        Parameters:
        -----------
            training_data_json : str
                File path to the training data <json-file>.
            in_degrees : bool
                Set `True` if the training data is in degrees
            bandwidth : float
                Bandwidth of the Gaussian kernel density estimator

        Returns:
        --------
            model : KdeStatePdf
                Probabilistic model of conformational states based on a Gaussian 
                kernel density estimator
        """
        with open(training_data_json, "r") as fhandle:
            training_data = json.load(fhandle)
        # Iterate over conformational states and train pdf estimator
        kde_list, lbl_list = [], []    
        for label, data in training_data.items():
            data = np.radians(data) if in_degrees else np.array(data)
            check_dihedral_range(data)
            # Fit Gaussian kernel density estimator
            kde = KernelDensity(bandwidth=bandwidth)
            kde.fit(data)
            kde_list.append(kde)
            lbl_list.append(label)
        return cls(state_labels=lbl_list, state_kdes=kde_list)


class ConfStateModelGrid(ConfStateModelABC):
    """Probabilistic model of conformational states based on a Gaussian kernel 
    density estimator. The actual inference is done by linear interpolation 
    between fixed grid points. This significantly speeds up the inference, while 
    sacrificing slightly on the accuracy of the inference.

    Attributes:
    -----------
        state_labels : Tuple
            List labels for the conformational state models
        state_grid : Array[M,N,N]
            Grid from which probability densities for conformational states are
            inferred by interpolation between gridpoints.
        grid_crds: Tuple[Array[N], Array[N]]
            Tuple of two arrays that describe the (phi,psi) coordinates of the
            grid points.
    
    Methods:
    --------
        get_logpdf(data)
            Inference of log-probability densities based on the probabilistic model.
        dump_pickle(output_file)
            Save the probabilistic as a pickle.

    Class methods:
    --------------
        from_pickle(pickled_file)
            Load conformational state models from pickled file
        from_fitting(json_file)
            Fit conformational state models to data provided in json file
    """
    model_type = "grid"

    def __init__(self, state_labels: List[str], state_grids: np.ndarray, grid_crds: Tuple):
        self.state_labels = tuple(state_labels)
        self.state_grids = state_grids
        self.grid_crds = grid_crds
    
    def get_logpdf(self, data: np.ndarray) -> np.ndarray:
        result = np.stack([
            interpn(self.grid_crds, grid, data) for grid in self.state_grids
        ])
        return result
    
    @classmethod
    def from_fitting(cls, training_data_json: str, *, in_degrees=False, bandwidth=.13, grid_points=10_000, **_):
        """Generate the probabilistic models at runtime, by fitting the models
        to the provided training data. Training data must be a json.
        
        Parameters:
        -----------
            training_data_json : str
                File path to the training data <json-file>.
            in_degrees : bool
                Set `True` if the training data is in degrees
            bandwidth : float
                Bandwidth of the Gaussian kernel density estimator
            grid_points : int
                The number of gridpoints between which the PDF will be estimated
                by interpolation. Note, that for generating the grid for both 
                axes sqrt(`grid_points`) are used. Thus, if `grid_points` is not 
                a square number, the final grid_points may be less.

        Returns:
        --------
            model : ConfStateModelGrid
                Probabilistic model of conformational states using a grid from
                which the PDF is estimated by interpolation
        """
        # Generate grid points in the (phi,psi)-space
        n = math.isqrt(grid_points)
        _phi, _psi = np.linspace(-np.pi, np.pi, n), np.linspace(-np.pi, np.pi, n)
        gridcrds = np.stack([
            arr.flatten() for arr in np.meshgrid(_phi, _psi, indexing="ij")],
            axis=1)
        # Infer PDF grid from a KDE model
        kde_model = ConfStateModelKDE.from_fitting(training_data_json, bandwidth=.13, in_degrees=in_degrees)
        _labels = kde_model.get_labels()
        _grids = np.reshape(kde_model.get_logpdf(gridcrds), (-1,n,n), order="C")
        return cls(state_labels=_labels, state_grids=_grids, grid_crds=(_phi, _psi))