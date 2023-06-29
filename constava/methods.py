import abc
import numpy as np
from typing import Optional


class ConstavaABC(metaclass=abc.ABCMeta):

    def calculate(self, state_logpdfs):
        subsampled_pdf = self._subsampling(state_logpdfs)
        state_propensities = self.calculateStatePropensities(subsampled_pdf)
        state_variability  = self.calculateStateVariability(subsampled_pdf)
        return state_propensities, state_variability

    def calculateStatePropensities(self, state_likelihoods):
        """
        Returns the average likelihoods for the residue to fall into the 
        configurational states 
        """
        return np.mean(state_likelihoods, axis=1)

    def calculateStateVariability(self, state_likelihoods):
        """ 
        Calculation of the variability in the conformational states. This is 
        calculated as the RMSF (root mean square fluctuation) of the state-
        likelihoods across all samples 
        """
        mean_likelihoods = np.mean(state_likelihoods, axis=1)
        squard_dev = np.sum((state_likelihoods.T - mean_likelihoods) ** 2, axis=1)
        state_var = np.sqrt(np.mean(squard_dev))
        return state_var

    # @abc.abstractmethod
    # def getLongName(self):
    #     pass

    @abc.abstractmethod
    def getShortName(self):
        pass

    @abc.abstractmethod
    def _subsampling(self, logpdf):
        pass


class ConstavaWindow(ConstavaABC):

    def __init__(self, window_size: int):
        self.window_size = window_size

    # def getLongName(self):
    #     return f"window({self.window_size:d})"

    def getShortName(self):
        return f"window({self.window_size:d})"

    def _subsampling(self, logpdf):
        # Subsampling using consecutive windows of window_size samples
        logpdf = np.stack([
            np.convolve(x, np.ones((self.window_size,)), mode="valid") 
            for x in logpdf])
        # Exponentiate and normalize to obtain likelihoods
        pdf = np.exp(logpdf)
        pdf /= np.sum(pdf, axis=0)
        return pdf


class ConstavaBootstrap(ConstavaABC):

    def __init__(self, sample_size: int, n_samples = 500, seed: Optional[int] = None):
        self.sample_size = sample_size
        self.n_samples = n_samples
        self.seed = seed

    # def getLongName(self):
    #     return f"bootstrapping({self.n_samples:d} samples of size {self.sample_size:d})"

    def getShortName(self):
        return f"bootstrapping({self.sample_size:d},{self.n_samples:d})"

    def _subsampling(self, logpdf):
        # Get dimensions of the input data
        n_states, n_measurements = logpdf.shape

        # Randomly select the <n_samples> samples by bootstrapping, with each
        # sample containing exactly <sample_size> measurements from the original
        # distribution. -> Array[n_states, n_samples, sample_size]
        np.random.seed(self.seed)
        samples = np.random.randint(n_measurements, size=self.sample_size*self.n_samples)
        logpdf = np.reshape(
            logpdf[:,samples], (n_states, self.n_samples, self.sample_size), 
            order="C")
        # Accumulate logpdfs within a sample = logpdf for each of these samples 
        # to be sampled from the same conformational state
        logpdf = np.sum(logpdf, axis=2)
        # Exponentiate and normalize to obtain likelihoods
        pdf = np.exp(logpdf)
        pdf /= np.sum(pdf, axis=0)
        return pdf