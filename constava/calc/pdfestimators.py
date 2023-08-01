""" constava.pdfestimators contains various definitions how the conformational 
states can be represented as probabilistic models.
"""

import abc
import json
import pickle
from typing import List, Tuple
import numpy as np
from scipy.interpolate import interpn
from sklearn.neighbors import KernelDensity
from ..io.ensemblereader import check_dihedral_range


class StatePdfABC(metaclass=abc.ABCMeta):
    """AbstractBaseClass for the probabilistic model describing the 
    conformational states.
    """    
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

    @abc.abstractmethod
    def dump_pickle(self, output_file: str):
        """Save the probabilistic as a pickle.
        
        Parameters:
        -----------
            output_file : str
                File path to store the model at.
        """
        pass

    @abc.abstractclassmethod
    def from_pickle(cls, pickled_file: str):
        """Load the probabilistic from a pickle.
        
        Parameters:
        -----------
            pickled_file : str
                File path to the pickled model.
        """
        pass

    @abc.abstractclassmethod
    def from_fitting(cls, training_data_json: str):
        """Generate the probabilistic models at runtime, by fitting the models
        to the provided training data. Training data must be a json.
        
        Parameters:
        -----------
            training_data_json : str
                File path to the training data <json-file>.
        """
        pass


class KDEStatePdf(StatePdfABC):
    """Probabilistic model of conformational states based on a Gaussian kernel 
    density estimator.

    Attributes:
    -----------
        kdes : Tuple
            List of the Gaussian kernel density estimators representing the 
            probabilistic models for the conformational states
        kde_labels: Tuple
            List labels for the conformational state models
    
    Methods:
    --------
        get_logpdf(data)
            Inference of log-probability densities based on the probabilistic model.
        dump_pickle(output_file)
            Save the probabilistic as a pickle.
    """
    def __init__(self, kdes: List[KernelDensity], kde_labels: List[str]):
        """Probabilistic model of conformational states based on a Gaussian kernel 
        density estimator.

        Parameters:
        -----------
            kdes : List or Tuple
                List of the Gaussian kernel density estimators representing the 
                probabilistic models for the conformational states
            kde_labels: List or Tuple
                List labels for the conformational state models
        """
        self.kdes = tuple(kdes)
        self.kde_labels = tuple(kde_labels)
    
    @property
    def labels(self):
        """Returns a copy of list of labels for the conformational state models"""
        return list(self.kde_labels)

    def get_logpdf(self, data: np.ndarray) -> np.ndarray:
        result = np.stack([
            kde.score_samples(data) for kde in self.kdes
        ])
        return result
    
    def dump_pickle(self, output_file: str):
        with open(output_file, "wb") as fhandle:
            pickle.dump({
                "labels": self.kde_labels,
                "kdes": self.kdes
            }, fhandle)
    
    @classmethod
    def from_pickle(cls, pickled_file: str):
        with open(pickled_file, "rb") as fhandle:
            _data = pickle.load(fhandle)
        return cls(kdes=_data["kdes"], kde_labels=_data["labels"])

    @classmethod
    def from_fitting(cls, training_data_json: str, *, bandwidth=.13, degrees2radians=False):
        with open(training_data_json, "r") as fhandle:
            training_data = json.load(fhandle)
        # Iterate over conformational states and train pdf estimator
        kdes, kde_labels = [], []    
        for label, data in training_data.items():
            if degrees2radians:
                data = np.radians(data)
            else:
                data = np.array(data)
            check_dihedral_range(data)
            kde = KernelDensity(bandwidth=bandwidth)
            kde.fit(data)
            kdes.append(kde)
            kde_labels.append(label)
        return cls(kdes=kdes, kde_labels=kde_labels)



class GridStatePdf(StatePdfABC):
    """Probabilistic model of conformational states based on a Gaussian kernel 
    density estimator. The actual inference is done by linear interpolation 
    between fixed grid points. This significantly speeds up the inference, while 
    sacrificing slightly on the accuracy of the inference.

    Attributes:
    -----------
        kdes : Tuple
            List of the Gaussian kernel density estimators representing the 
            probabilistic models for the conformational states
        kde_labels : Tuple
            List labels for the conformational state models
        grids : Array[M,Phi,Psi]
            Grid from which probability densities are inferred by interpolation 
            between gridpoints. If not provided, generated at runtime.
        gridcrds : Tuple(Array[Phi], Array[Psi])
            Tuple of two arrays that indicate the coordinates in the 
            (phi,psi)-space across the grid.

    Methods:
    --------
        get_logpdf(data)
            Inference of log-probability densities based on the probabilistic model.
        dump_pickle(output_file)
            Save the probabilistic as a pickle.
    """
    DEFAULT_COORDINATES = np.linspace(-np.pi, np.pi, 100)

    def __init__(self, kdes: List[KernelDensity], kde_labels: List[str], 
                 grids: np.ndarray = None, gridcrds: Tuple[np.ndarray] = None):
        """Probabilistic model of conformational states. The inference is done 
        by interpolation on a grid either provided or inferred from a Gaussian 
        kernel density estimator.

        Parameters:
        -----------
            kdes : List or Tuple
                List of the Gaussian kernel density estimators representing the 
                probabilistic models for the conformational states
            kde_labels: List or Tuple
                List labels for the conformational state models
            grids : Array[M,Phi,Psi]
                Grid from which probability densities are inferred by interpolation 
                between gridpoints. If not provided, generated at runtime.
            gridcrds : Tuple(Array[Phi], Array[Psi])
                Tuple of two arrays that indicate the coordinates in the 
                (phi,psi)-space across the grid.
        """
        self.kdes = tuple(kdes)
        self.kde_labels = tuple(kde_labels)
        if grids is None or gridcrds is None:
            self.grids, self.gridcrds =  self.construct_grids()
        else:
            self.grids = grids
            self.gridcrds = gridcrds

    @property
    def labels(self):
        """Returns a copy of list of labels for the conformational state models"""
        return list(self.kde_labels)

    def get_logpdf(self, data: np.ndarray) -> np.ndarray:
        result = np.stack([
            interpn(self.gridcrds, grid, data) for grid in self.grids])
        return result

    def construct_grids(self, phicrd: np.ndarray = None, psicrd: np.ndarray = None):
        """Function to generate the grid for the inference at runtime.

        Parameters:
        -----------
            phicrd: Array[N]
            psicrd: Array[M]
                Coordinates of phi/psi in the dihedral space. If not provided,
                `np.linspace(-np.pi, np.pi, 100)`will be used as default.

        Returns:
        --------
            grids
                Grid from which probability densities are inferred by interpolation 
                between gridpoints. If not provided, generated at runtime.
            (phicrd, psicrd)
                Tuple of two arrays that indicate the coordinates in the 
                (phi,psi)-space across the grid.
        """
        if phicrd is None:
            phicrd = self.DEFAULT_COORDINATES
        if psicrd is None:
            psicrd = phicrd
        gridcrd = np.stack(list(
            map(lambda arr: arr.flatten(), np.meshgrid(phicrd, psicrd, indexing="ij"))), 
            axis=1)
        grids = np.stack([
            kde.score_samples(gridcrd) for kde in self.kdes
        ])
        grids = grids.reshape((len(self.kdes),) + phicrd.shape + psicrd.shape, order="C")
        return grids, (phicrd, psicrd)

    def dump_pickle(self, output_file: str):
        with open(output_file, "wb") as fhandle:
            pickle.dump({
                "labels": self.kde_labels,
                "kdes": self.kdes,
                "grids": self.grids,
                "gridcrds": self.gridcrds,
            }, fhandle)

    @classmethod
    def from_pickle(cls, pickled_file: str):
        with open(pickled_file, "rb") as fhandle:
            data = pickle.load(fhandle)
        return cls(kdes=data["kdes"], kde_labels=data["labels"], 
                   grids=data.get("grids", None), gridcrds=data.get("gridcrds", None))
    
    @classmethod
    def from_fitting(cls, training_data_json: str, *, bandwidth=.13, degrees2radians=False):
        """ Fits a new KDE on the given training data """
        with open(training_data_json, "r") as fhandle:
            training_data = json.load(fhandle)
        # Iterate over conformational states and train pdf estimator
        kdes, kde_labels = [], []    
        for label, data in training_data.items():
            if degrees2radians:
                data = np.radians(data)
            else:
                data = np.array(data)
            check_dihedral_range(data)
            kde = KernelDensity(bandwidth=bandwidth)
            kde.fit(data)
            kdes.append(kde)
            kde_labels.append(label)
        return cls(kdes=kdes, kde_labels=kde_labels)