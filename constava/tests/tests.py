"""module for unit testing constava"""
import os, glob, filecmp
import tarfile
import tempfile
import unittest
from constava import Constava
from constava.utils.constants import CONSTAVA_DATA_DIR
from constava.calc.csmodels import ConfStateModelLoadingError, ConfStateModelKDE, ConfStateModelGrid
from constava.utils.logging import logging

logger = logging.getLogger("Constava")


class TestWrapper(unittest.TestCase):
    """Class that runs test runs on the wrapper as a whole"""

    TEST_SOURCE = os.path.join(CONSTAVA_DATA_DIR, "constava_testdata.tgz")
    TEST_TEMPORARY_DIRECTORY = None
    TEST_TEMPDIR = None
    TEST_MODELDUMP0 = None
    TEST_MODELDUMP1 = None

    @classmethod
    def setUpClass(cls):
        """Setup method that generates a temporary directory and some file paths
        that are used by all tests"""
        if cls.TEST_TEMPORARY_DIRECTORY is None:
            # Create the temporary directory and files for tests
            cls.TEST_TEMPORARY_DIRECTORY = tempfile.TemporaryDirectory(prefix="ConstavaTest.")
            cls.TEST_TEMPDIR = cls.TEST_TEMPORARY_DIRECTORY.name
            _, cls.TEST_MODELDUMP0 = tempfile.mkstemp(prefix="model.", suffix=".pkl", dir=cls.TEST_TEMPDIR)
            _, cls.TEST_MODELDUMP1 = tempfile.mkstemp(prefix="model.", suffix=".pkl", dir=cls.TEST_TEMPDIR)
            logger.warning(f"TEST PREPARATION: Creating temporary directory: {cls.TEST_TEMPDIR}")

            # Extract the test files into the new temporary directory
            logger.warning(f"TEST PREPARATION: Untarring test data in: {cls.TEST_TEMPDIR}")
            with tarfile.open(cls.TEST_SOURCE, mode="r:gz") as tarchive:
                tarchive.extractall(cls.TEST_TEMPDIR)

    @classmethod
    def get_test_count(cls):
        if not hasattr(cls, "TEST_COUNTER"):
            cls.TEST_COUNTER = 1
        else:
            cls.TEST_COUNTER += 1
            print()
        return cls.TEST_COUNTER
    
    def test_0a_KDEModelFitting(self):
        """Test fitting of a KDE model, dumps model for further tests"""
        logger.warning(f"TEST #{self.get_test_count()}: Fitting of 'kde' model...")
        cva = Constava(verbose=0)
        csmodel = cva.fit_csmodel(kde_bandwidth=.1291)
        self.assertIsInstance(csmodel, ConfStateModelKDE, "Failed to fit conformational state model.")

        logger.warning(f"TEST #{self.get_test_count()}: Caching of 'kde' model...")
        cva.fit_csmodel(kde_bandwidth=.1291)
        self.assertIs(csmodel, cva._csmodel, "Failed to reuse preloaded conformational state model.")

        logger.warning(f"TEST #{self.get_test_count()}: Fitting a new 'kde' model...")
        cva.fit_csmodel(kde_bandwidth=.1)
        self.assertIsNot(csmodel, cva._csmodel, "Failed to update conformational state model after parameter change.")

        logger.warning(f"TEST #{self.get_test_count()}: Storing of pickled 'kde' model...")
        logger.warning(f"... Dumping model in: {self.TEST_MODELDUMP0}")
        cva._csmodel.dump_pickle(self.TEST_MODELDUMP0)
        self.assertTrue(os.path.isfile(self.TEST_MODELDUMP0), "Failed to conformational state dump model.")

    def test_1a_KDEModelLoading(self):
        """Test loading of pre-fitted KDE model"""
        logger.warning(f"TEST #{self.get_test_count()}: Loading of pickled 'kde' models...")
        cva = Constava(verbose=0, model_load=self.TEST_MODELDUMP0)
        cva.load_csmodel(pickled_csmodel=self.TEST_MODELDUMP0)
        self.assertIsInstance(cva._csmodel, ConfStateModelKDE, "Failed to fit conformational state model.")

    def test_2a_KDEModelInference(self):
        """Test inference from pre-fitted KDE model"""
        logger.warning(f"TEST #{self.get_test_count()}: Inference from 'kde' conformational state model...")
        input_files = glob.glob(f"{self.TEST_TEMPDIR}/xvg/ramaPhiPsi*.xvg")
        expected_result = f"{self.TEST_TEMPDIR}/xvg/result_kde.csv"
        output_file  = tempfile.NamedTemporaryFile(prefix="kde.", suffix=".csv", dir=self.TEST_TEMPDIR)
        cva = Constava(
            input_files = input_files,
            output_file = output_file.name,
            model_type = "kde",
            model_load = self.TEST_MODELDUMP0,
            window = [1,3,7,23], 
            bootstrap = [3,7,23],
            window_series = [1,7],
            seed = 42,
            verbose = 0, 
            input_degrees=True)
        cva.run()
        self.assertTrue(filecmp.cmp(output_file.name, expected_result))

    def test_0b_GridModelFitting(self):
        """Test fitting of a grid-interpolation model, dumps model for further tests"""
        logger.warning(f"TEST #{self.get_test_count()}: Fitting of 'grid' models...")
        cva = Constava(verbose=0)
        csmodel = cva.fit_csmodel(model_type="grid", kde_bandwidth=.42, grid_points=145)
        self.assertIsInstance(csmodel, ConfStateModelGrid, "Failed to fit conformational state model.")

        logger.warning(f"TEST #{self.get_test_count()}: Caching of 'grid' models...")
        cva.fit_csmodel(model_type="grid", kde_bandwidth=.42, grid_points=145)
        self.assertIs(csmodel, cva._csmodel, "Failed to reuse preloaded conformational state model.")

        logger.warning(f"TEST #{self.get_test_count()}: Fitting a new 'grid' model...")
        cva.fit_csmodel(model_type="grid", kde_bandwidth=.1, grid_points=3601)
        self.assertIsNot(csmodel, cva._csmodel, "Failed to update conformational state model after parameter change.")

        logger.warning(f"TEST #{self.get_test_count()}: Storing of pickled 'grid' models...")
        logger.warning(f"... Dumping model in: {self.TEST_MODELDUMP1}")
        cva._csmodel.dump_pickle(self.TEST_MODELDUMP1)
        self.assertTrue(os.path.isfile(self.TEST_MODELDUMP1), "Failed to conformational state dump model.")

    def test_1b_GridModelLoading(self):
        """Test loading of pre-fitted grid-inference model"""
        logger.warning(f"TEST #{self.get_test_count()}: Loading of pickled 'grid' models...")
        cva = Constava(verbose=0, model_load=self.TEST_MODELDUMP1)
        cva.load_csmodel(pickled_csmodel=self.TEST_MODELDUMP1)
        self.assertIsInstance(cva._csmodel, ConfStateModelGrid, "Failed to fit conformational state model.")

    def test_2b_GridModelInference(self):
        """Test inference from pre-fitted grid-inference model"""
        logger.warning(f"TEST #{self.get_test_count()}: Inference from 'grid' conformational state model...")
        input_files = f"{self.TEST_TEMPDIR}/csv/dihedrals.csv"
        expected_result = f"{self.TEST_TEMPDIR}/csv/result_grid.csv"
        output_file  = tempfile.NamedTemporaryFile(prefix="grid.", suffix=".csv", dir=self.TEST_TEMPDIR)
        c = Constava(
            input_files = input_files,
            output_file = output_file.name,
            model_type = "grid",
            model_load = self.TEST_MODELDUMP1,
            window = [1,3,7,23], 
            bootstrap = [3,7,23],
            window_series = [1,7],
            seed = 42,
            verbose = 0, 
            input_degrees=False)
        c.run()
        self.assertTrue(filecmp.cmp(output_file.name, expected_result))

def run_unittest():
        """Run the unittests in this module"""
        suite = unittest.TestSuite((
            TestWrapper("test_0a_KDEModelFitting"),
            TestWrapper("test_1a_KDEModelLoading"),
            TestWrapper("test_2a_KDEModelInference"),
            TestWrapper("test_0b_GridModelFitting"),
            TestWrapper("test_1b_GridModelLoading"),
            TestWrapper("test_2b_GridModelInference"),
        ))
        runner = unittest.TextTestRunner()
        runner.run(suite)

if __name__ == "__main__":
    run_unittest()