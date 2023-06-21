import abc
import json
import pickle
from typing import Generator, List, Tuple
import numpy as np
from sklearn.neighbors import KernelDensity


class StatePdfABC(metaclass=abc.ABCMeta):
    
    @abc.abstractmethod
    def get_logpdf(self, data: np.ndarray) -> np.ndarray:
        pass

    @abc.abstractmethod
    def dump_pickle(self, output_file: str):
        pass

    @abc.abstractclassmethod
    def from_pickle(self, pickled_file: str):
        pass

    @abc.abstractclassmethod
    def from_fitting(self, training_data_json: str):
        pass


class KDEStatePdf(StatePdfABC):

    def __init__(self, kdes: List, kde_labels: List):
        self.kdes = tuple(kdes)
        self.kde_labels = tuple(kde_labels)
    
    @property
    def labels(self):
        return list(self.kde_labels)

    # def __iter__(self) -> Generator[Tuple[str, KernelDensity], None, None]:
    #     return zip(self.kde_labels, self.kdes)

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
            kde = KernelDensity(bandwidth=bandwidth)
            kde.fit(data)
            kdes.append(kde)
            kde_labels.append(label)
        return cls(kdes=kdes, kde_labels=kde_labels)


class GridStatePdf(KDEStatePdf):
    pass