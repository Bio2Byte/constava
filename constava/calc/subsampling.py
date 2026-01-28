"""constava.calc.methods contains classes representing sub-sampling schemes fro the
probability state propensity calculations. Multiple of these schemes can be
dynamically defined on run-time."""

import abc
from typing import Optional
import numpy as np


class SubsamplingMethodError(ValueError):
    """Raised when there are no sub-sampling methods passed to the calculator"""


class SubsamplingABC(metaclass=abc.ABCMeta):
    """Base class to subsample the logPDF values obtained from the probabilistic
    conformational state models and calculate the conformational state
    propensities and conformational state variability.

    Methods:
    --------
        calculate(state_logpdfs)
            Calculates the conformational state likelihoods and conformational
            state variability.
        calculateStatePropensities(state_likelihoods)
            Calculates the average conformational state likelihood.
        calculateStateVariability(state_likelihoods)
            Calculates the  conformational state variability.
        getShortName()
            Name of the method for reference in the output.
        _subsampling(state_logpdfs)
            Sub-samples from the distribution of original data points.
    """

    def calculate(self, state_logpdfs):
        """Calculates the conformational state likelihoods and conformational
        state variability from the sampled state logPDFs.

        Parameters:
        -----------
            state_logpdfs : Array[M, N]
                                An array of the logPDFs of N samples for M states

        Returns:
        --------
            state_propensities : Array[M]
                                Average likelihood for samples to fall in any
                                of the M states
            state_variability : float
                                variability fo the state propensities
                                throughout the sampling

        """
        subsampled_pdf = self._subsampling(state_logpdfs)

        state_propensities = self.calculateStatePropensities(subsampled_pdf)
        state_variability = self.calculateStateVariability(subsampled_pdf)

        return state_propensities, state_variability

    def calculateStatePropensities(self, state_likelihoods):
        """Calculates the average conformational state likelihood."""
        return np.mean(state_likelihoods, axis=1)

    def calculateStateVariability(self, state_likelihoods):
        """Calculation of the variability in the conformational states. This is
        calculated as the RMSF (root mean square fluctuation) of the state
        likelihoods across all samples.

        Parameters:
        -----------
            state_likelihoods : Array[M,N]
                Likelihoods for each of the M states along N samples.

        Returns:
        --------
            state_var : float
                Conformational state variability
        """
        mean_likelihoods = np.mean(state_likelihoods, axis=1)
        squard_dev = np.sum((state_likelihoods.T - mean_likelihoods) ** 2, axis=1)
        state_var = np.sqrt(np.mean(squard_dev))

        return state_var

    @abc.abstractmethod
    def getShortName(self) -> str:
        """Name of the method for reference in the output."""

    @abc.abstractmethod
    def _subsampling(self, logpdf):
        """Method used to subsample from the distribution of logPDF values and
        convert them into individual likelihoods for each conformational state
        model.

        Parameters:
        -----------
            logpdf : Array[M,X]
                log-probability densities for M states across X original data points.

        Retruns:
        --------
            pdf : Array[M,N]
                Likelihoods obtained by sub-sampling N times using the described method.
        """


class SubsamplingWindow(SubsamplingABC):
    """Class to subsample the logPDF values obtained from the probabilistic
    conformational state models and calculate the conformational state
    propensities and conformational state variability. Sub-sampling is done
    using a sliding window.

    Attributes:
    -----------
        window_size : int
            Size of the sliding window used in sub-sampling.

    Methods:
    --------
        calculate(state_logpdfs)
            Calculates the conformational state likelihoods and conformational
            state variability.
        calculateStatePropensities(state_likelihoods)
            Calculates the average conformational state likelihood.
        calculateStateVariability(state_likelihoods)
            Calculates the  conformational state variability.
        getShortName()
            Name of the method for reference in the output.
        _subsampling(state_logpdfs)
            Sub-samples from the distribution of original data points.
    """

    def __init__(self, window_size: int):
        """Initialize class to subsample and calculate conformational state
        propensities and conformational state variability. Sub-sampling is done
        using a sliding window.

        Parameters:
        -----------
            window_size : int
                Size of the sliding window used in sub-sampling.
        """
        self.window_size = window_size

    def getShortName(self) -> str:
        """Name of the method for reference in the output."""
        return f"window/{self.window_size:d}/"

    def _subsampling(self, logpdf):
        """Sub-sampling from the distribution of logPDF using a sliding window of
        size `window_size`. With a `window_size == 1`, this effectively uses the
        original data points. Finally, likelihoods for each conformational state
        model are calculated for each sample.

        Parameters:
        -----------
            logpdf : Array[M,X]
                log-Probability densities for M states across X original data points.

        Returns:
        --------
            pdf : Array[M,N]
                Likelihoods obtained by sub-sampling N times using the described method.
        """

        # Sub-sampling using consecutive windows of window_size samples
        convolved_logpdf = np.stack(
            [np.convolve(x, np.ones((self.window_size,)), mode="valid") for x in logpdf]
        )

        # Exponentiation and normalize to obtain likelihoods
        pdf = np.exp(convolved_logpdf)
        pdf /= np.sum(pdf, axis=0)

        return pdf

    def __str__(self):
        return "SubsamplingWindow"


class SubsamplingBootstrap(SubsamplingABC):
    """Class to subsample the logPDF values obtained from the probabilistic
    conformational state models and calculate the conformational state
    propensities and conformational state variability. Sub-sampling is done
    using bootstrapping.

    Attributes:
    -----------
        sample_size : int
            Number of original data points in each bootstrapped sample.
        n_samples : int
            Number of samples to bootstrap.
        seed: int
            Random seed used during bootstrapping

    Methods:
    --------
        calculate(state_logpdfs)
            Calculates the conformational state likelihoods and conformational
            state variability.
        calculateStatePropensities(state_likelihoods)
            Calculates the average conformational state likelihood.
        calculateStateVariability(state_likelihoods)
            Calculates the  conformational state variability.
        getShortName()
            Name of the method for reference in the output.
        _subsampling(state_logpdfs)
            Sub-samples from the distribution of original data points.
    """

    def __init__(self, sample_size: int, n_samples=500, seed: Optional[int] = None):
        """Initialize class to subsample and calcualte conformational state
        propensities and conformational state variability. Sub-sampling is done
        using bootstrapping.

        Parameters:
        -----------
            sample_size : int
                Number of original data points in each bootstrapped sample.
            n_samples : int
                Number of samples to bootstrap.
            seed: int
                Random seed used during bootstrapping
        """
        self.sample_size = sample_size
        self.n_samples = n_samples
        self.seed = seed

    def getShortName(self) -> str:
        """Name of the method for reference in the output."""
        return "bootstrap/{0:d}/{1:d}/{2}/".format(
            self.sample_size, self.n_samples, self.seed or ""
        )

    def _subsampling(self, logpdf):
        """Sub-sampling from the distribution of logPDF using bootstrapping.
        `n_samples` are sub-sampled from the distribution, where each sample
        contains `sample_size`randomly selected original data points. Finally,
        likelihoods for each conformational state model are calculated for each
        sample.

        Parameters:
        -----------
            logpdf : Array[M,X]
                The log-Probability densities for M states across X
                original data points.

        Returns:
        --------
            pdf : Array[M,N]
                Likelihoods obtained by sub-sampling N times using the described method.
        """
        # Get dimensions of the input data
        n_states, n_measurements = logpdf.shape

        # Randomly select the <n_samples> samples by bootstrapping, with each
        # sample containing exactly <sample_size> measurements from the original
        # distribution. -> Array[n_states, n_samples, sample_size]
        rng = np.random.default_rng(self.seed)
        samples = rng.integers(n_measurements, size=self.sample_size * self.n_samples)

        reshaped_logpdf = np.reshape(
            logpdf[:, samples], (n_states, self.n_samples, self.sample_size), order="C"
        )

        # Accumulate logpdfs within a sample = logpdf for each of these samples
        # to be sampled from the same conformational state
        aggregated_logpdf = np.sum(reshaped_logpdf, axis=2)

        # Exponential and normalize to obtain likelihoods
        pdf = np.exp(aggregated_logpdf)
        pdf /= np.sum(pdf, axis=0)

        return pdf

    def __str__(self):
        return "SubsamplingBootstrap"


class SubsamplingWindowSeries(SubsamplingWindow):
    """Class to subsample the logPDF values obtained from the probabilistic
    conformational state models and calculate the conformational state
    propensities and conformational state variability. Sub-sampling is done
    using a sliding window. For each window the results are returned.

    Attributes:
    -----------
        window_size : int
            Size of the sliding window used in sub-sampling.

    Methods:
    --------
        calculate(state_logpdfs)
            Calculates the conformational state likelihoods and conformational
            state variability.
        calculateStatePropensities(state_likelihoods)
            Calculates the samples' conformational state likelihood.
        calculateStateVariability(state_likelihoods)
            Calculates the conformational state variability.
        getShortName()
            Name of the method for reference in the output.
        _subsampling(state_logpdfs)
            Sub-samples from the distribution of original data points.
    """

    def getShortName(self) -> str:
        """Name of the method for reference in the output."""
        return f"window_series/{self.window_size:d}/"

    def calculateStatePropensities(self, state_likelihoods):
        """Calculates the conformational state likelihoods for the given sample."""
        return state_likelihoods

    def calculateStateVariability(self, state_likelihoods):
        """Calculates distance of the conformational states of each sample to
        the average conformational state.

        Parameters:
        -----------
            state_likelihoods : Array[M,N]
                Likelihoods for each of the M states along N samples.

        Returns:
        --------
            state_var : Array[N]
                Conformational state distances from the average
        """
        mean_likelihoods = np.mean(state_likelihoods, axis=1)
        squard_dev = np.sum((state_likelihoods.T - mean_likelihoods) ** 2, axis=1)
        state_var = np.sqrt(squard_dev)
        return state_var


class SubsamplingBootstrapSeries(SubsamplingBootstrap):
    """Class to subsample the logPDF values obtained from the probabilistic
    conformational state models and calculate the conformational state
    propensities and conformational state variability. Sub-sampling is done
    using bootstrapping. For each bootstrapped subsample the results are
    returned.

    Attributes:
    -----------
        sample_size : int
            Number of original data points in each bootstrapped sample.
        n_samples : int
            Number of samples to bootstrap.
        seed: int
            Random seed used during bootstrapping

    Methods:
    --------
        calculate(state_logpdfs)
            Calculates the conformational state likelihoods and conformational
            state variability.
        calculateStatePropensities(state_likelihoods)
            Calculates the samples' conformational state likelihood.
        calculateStateVariability(state_likelihoods)
            Calculates the conformational state variability.
        getShortName()
            Name of the method for reference in the output.
        _subsampling(state_logpdfs)
            Sub-samples from the distribution of original data points.
    """

    def getShortName(self) -> str:
        """Name of the method for reference in the output."""
        return "bootstrap_series/{0:d}/{1:d}/{2}/".format(
            self.sample_size, self.n_samples, self.seed or ""
        )

    def calculateStatePropensities(self, state_likelihoods):
        """Calculates the conformational state likelihoods for the given sample."""
        return state_likelihoods

    def calculateStateVariability(self, state_likelihoods):
        """Calculates distance of the conformational states of each sample to
        the average conformational state.

        Parameters:
        -----------
            state_likelihoods : Array[M,N]
                Likelihoods for each of the M states along N samples.

        Returns:
        --------
            state_var : Array[N]
                Conformational state distances from the average
        """
        mean_likelihoods = np.mean(state_likelihoods, axis=1)

        squard_dev = np.sum((state_likelihoods.T - mean_likelihoods) ** 2, axis=1)
        state_var = np.sqrt(squard_dev)

        return state_var
