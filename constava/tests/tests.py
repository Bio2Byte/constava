"""module for unit testing constava"""
import os, glob, filecmp
import tarfile
import tempfile
import unittest
from constava import Constava
from constava.utils.constants import CONSTAVA_DATA_DIR
from constava.calc.csmodels import ConfStateModelLoadingError, ConfStateModelKDE, ConfStateModelGrid
from constava.wrapper.logging import logging

logger = logging.getLogger("Constava")


class TestWrapper(unittest.TestCase):
    """Class that runs test runs on the wrapper as a whole"""

    TEST_SOURCE = os.path.join(CONSTAVA_DATA_DIR, "constava_testdata.tgz")
    
    @classmethod
    def setUpClass(cls):
        if not hasattr(cls, "temporary_directory"):
            # Create the temporary directory to run the tests in
            cls.temporary_directory = tempfile.TemporaryDirectory(prefix="ConstavaTest.")
            cls.tmpdir = cls.temporary_directory.name
            logger.warning(f"TEST: Creating temporary directory: {cls.tmpdir}")
            # Extract the test files into the new temporary directory
            with tarfile.open(cls.TEST_SOURCE, mode="r:gz") as tarchive:
                tarchive.extractall(cls.tmpdir)
        if not hasattr(cls, "kde_model"):
            logger.warning("TEST: Fitting 'kde' conformational state model for testing...")
            cls.kde_model = cls._fit_model("kde", kde_bandwidth=.10)
            logger.warning(f"TEST: ... Fitted model written to: {cls.kde_model}")
        if not hasattr(cls, "grid_model"):
            logger.warning("TEST: Fitting 'grid' conformational state model for testing...")
            cls.grid_model = cls._fit_model("grid", kde_bandwidth=.10, grid_points=3601)
            logger.warning(f"TEST: ... Fitted model written to: {cls.grid_model}")

    @classmethod
    def _fit_model(cls, model_type="kde", **model_kwargs):
        _, dumpfile = tempfile.mkstemp(prefix="model.", suffix=".pkl", dir=cls.tmpdir)
        c = Constava(verbose=1)
        csmodel = c.fit_csmodel(model_type, **model_kwargs)
        csmodel.dump_pickle(dumpfile)
        return dumpfile
    
    @classmethod
    def get_test_count(cls):
        if not hasattr(cls, "TEST_COUNTER"):
            cls.TEST_COUNTER = 1
        else:
            cls.TEST_COUNTER += 1
            print()
        return cls.TEST_COUNTER

    def test_csmodel_loading_kde(self):
        logger.warning(f"TEST #{self.get_test_count()}: Testing loading procedures for 'kde' models...")
        cva = Constava(verbose=1)
        # It should fail to load a 'grid' model as kde
        # with self.assertRaises(ConfStateModelLoadingError):
        #     c.initialize_calculator(model_type = "kde",
        #                             load_model = self.grid_model)
        # Load the KDE model
        csmodel0 = cva.load_csmodel(pickled_csmodel = self.kde_model)
        self.assertIsInstance(csmodel0, ConfStateModelKDE, "Failed to load conformational state model.")
        # Here, a new model should be fittet
        csmodel1 = cva.fit_csmodel(kde_bandwidth=.1291)
        self.assertIsNot(csmodel0, cva._csmodel, "Failed to fit new model after loading a model.")
        # Here, the same model should be reused
        cva.fit_csmodel(kde_bandwidth=.1291)
        self.assertIs(csmodel1, cva._csmodel, "Failed to reuse preloaded conformational state model.")
        # Here, a new model should be fitted
        csmodel2 = cva.fit_csmodel(kde_bandwidth=.13)
        self.assertIsNot(csmodel1, cva._csmodel, "Failed to fit a newmodel after changing parameters")
        # Here, a new model should be loaded
        cva.load_csmodel(pickled_csmodel = self.kde_model)
        self.assertIsNot(csmodel2, cva._csmodel, "Failed to load new model after fitting a model.")

    def test_csmodel_loading_grid(self):
        logger.warning(f"TEST #{self.get_test_count()}: Testing loading procedures for 'grid' models...")
        cva = Constava(verbose=1)
        # # It should fail to load a 'grid' model as kde
        # with self.assertRaises(ConfStateModelLoadingError):
        #     c.initialize_calculator(model_type = "grid",
        #                             load_model = self.kde_model)
        # Load the KDE model
        cva.load_csmodel(pickled_csmodel = self.grid_model)
        csmodel0, cshash0 = cva._csmodel, cva._cshash
        self.assertIsInstance(csmodel0, ConfStateModelGrid, "Failed to load conformational state model.")
        # Here the same csmodel should be used
        cva.load_csmodel(pickled_csmodel = self.grid_model)
        self.assertIs(csmodel0, cva._csmodel, "Failed to reuse preloaded conformational state model.")

    def test_inference_kde(self):
        logger.warning(f"TEST #{self.get_test_count()}: Testing inference for 'kde' models...")
        input_files = glob.glob(f"{self.tmpdir}/xvg/ramaPhiPsi*.xvg")
        expected_result = f"{self.tmpdir}/xvg/result_kde.csv"
        output_file  = tempfile.NamedTemporaryFile(prefix="kde.", suffix=".csv", dir=self.tmpdir)
        c = Constava(
            input_files = input_files,
            output_file = output_file.name,
            model_type = "kde",
            model_load = self.kde_model,
            window = [1,3,7,23], 
            bootstrap = [3,7,23],
            window_series = [1,7],
            seed = 42,
            verbose = 1, 
            input_degrees=True)
        c.run()
        self.assertTrue(filecmp.cmp(output_file.name, expected_result))

    def test_inference_grid(self):
        logger.warning(f"TEST #{self.get_test_count()}: Testing inference for 'grid' models...")
        input_files = [f"{self.tmpdir}/csv/dihedrals.csv"]
        expected_result = f"{self.tmpdir}/csv/result_grid.csv"
        output_file  = tempfile.NamedTemporaryFile(prefix="grid.", suffix=".csv", dir=self.tmpdir)
        c = Constava(
            input_files = input_files,
            output_file = output_file.name,
            model_type = "grid",
            model_load = self.grid_model,
            window = [1,3,7,23], 
            bootstrap = [3,7,23],
            window_series = [1,7],
            seed = 42,
            verbose = 1, 
            input_degrees=False)
        c.run()
        self.assertTrue(filecmp.cmp(output_file.name, expected_result))


if __name__ == "__main__":
    unittest.main()