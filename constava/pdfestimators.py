import abc
import json
import pickle
from typing import List, Tuple
import numpy as np
from scipy.interpolate import interpn
from sklearn.neighbors import KernelDensity
from constava.ensemblereader import check_dihedral_range


class StatePdfABC(metaclass=abc.ABCMeta):
    """ 
    The AbstractBaseClass for all functions providing probability density 
    function estimates for various states.
    """
    
    @abc.abstractmethod
    def get_logpdf(self, data: np.ndarray) -> np.ndarray:
        pass

    @abc.abstractmethod
    def dump_pickle(self, output_file: str):
        pass

    @abc.abstractclassmethod
    def from_pickle(cls, pickled_file: str):
        pass

    @abc.abstractclassmethod
    def from_fitting(cls, training_data_json: str):
        pass


class KDEStatePdf(StatePdfABC):

    def __init__(self, kdes: List[KernelDensity], kde_labels: List[str]):
        self.kdes = tuple(kdes)
        self.kde_labels = tuple(kde_labels)
    
    @property
    def labels(self):
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



class GridStatePdf(StatePdfABC):
    """ A subclass of KDEStatePdf designed for faster inference. """

    DEFAULT_COORDINATES = np.linspace(-np.pi, np.pi, 100)

    def __init__(self, kdes: List[KernelDensity], kde_labels: List[str], 
                 grids: np.ndarray = None, gridcrds: Tuple[np.ndarray] = None):
        self.kdes = tuple(kdes)
        self.kde_labels = tuple(kde_labels)
        if grids is None or gridcrds is None:
            self.grids, self.gridcrds =  self.construct_grids()
        else:
            self.grids = grids
            self.gridcrds = gridcrds

    @property
    def labels(self):
        return list(self.kde_labels)

    def get_logpdf(self, data: np.ndarray) -> np.ndarray:
        result = np.stack([
            interpn(self.gridcrds, grid, data) for grid in self.grids])
        return result

    def construct_grids(self, phicrd: np.ndarray = None, psicrd: np.ndarray = None):
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