"""module for unit testing constava"""

import json
import os
import glob
import filecmp
import tempfile
import unittest
from constava import Constava
from constava.utils.constants import CONSTAVA_DATA_DIR
from constava.calc.csmodels import (
    ConfStateModelKDE,
    ConfStateModelGrid,
)
from constava.utils.logging import logging

logger = logging.getLogger("Constava")


class TestWrapper(unittest.TestCase):
    """Class that runs unit tests as a whole"""

    def __init__(self, testName, verbosity):
        super(TestWrapper, self).__init__(testName)
        self.verbosity = verbosity

    TEST_SOURCE = os.path.join(CONSTAVA_DATA_DIR, "constava_testdata.tgz")
    TEST_TEMPORARY_DIRECTORY = None
    TEST_TEMPDIR = None
    TEST_MODELDUMP0 = None
    TEST_MODELDUMP1 = None
    CONSTAVA_TEST_DATA_DIR = None

    @classmethod
    def setUpClass(cls):
        """Setup method that generates a temporary directory and some file paths
        that are used by all tests"""
        if cls.TEST_TEMPORARY_DIRECTORY is None:
            # Create the temporary directory and files for tests
            cls.TEST_TEMPORARY_DIRECTORY = tempfile.TemporaryDirectory(
                prefix="ConstavaTest."
            )
            cls.TEST_TEMPDIR = cls.TEST_TEMPORARY_DIRECTORY.name
            _, cls.TEST_MODELDUMP0 = tempfile.mkstemp(
                prefix="model.", suffix=".pkl", dir=cls.TEST_TEMPDIR
            )
            _, cls.TEST_MODELDUMP1 = tempfile.mkstemp(
                prefix="model.", suffix=".pkl", dir=cls.TEST_TEMPDIR
            )
            cls.CONSTAVA_TEST_DATA_DIR = os.path.join(
                CONSTAVA_DATA_DIR, "constava_testdata"
            )

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
        cva = Constava(verbose=self.verbosity)
        csmodel = cva.fit_csmodel(kde_bandwidth=0.1291)
        self.assertIsInstance(
            csmodel, ConfStateModelKDE, "Failed to fit conformational state model."
        )

        logger.warning(f"TEST #{self.get_test_count()}: Caching of 'kde' model...")
        cva.fit_csmodel(kde_bandwidth=0.1291)
        self.assertIs(
            csmodel,
            cva._csmodel,
            "Failed to reuse preloaded conformational state model.",
        )

        logger.warning(f"TEST #{self.get_test_count()}: Fitting a new 'kde' model...")
        cva.fit_csmodel(kde_bandwidth=0.1)
        self.assertIsNot(
            csmodel,
            cva._csmodel,
            "Failed to update conformational state model after parameter change.",
        )

        logger.warning(
            f"TEST #{self.get_test_count()}: Storing of pickled 'kde' model..."
        )
        logger.warning(f"... Dumping model in: {self.TEST_MODELDUMP0}")
        cva._csmodel.dump_pickle(self.TEST_MODELDUMP0)
        self.assertTrue(
            os.path.isfile(self.TEST_MODELDUMP0),
            "Failed to conformational state dump model.",
        )

    def test_1a_KDEModelLoading(self):
        """Test loading of pre-fitted KDE model"""
        logger.warning(
            f"TEST #{self.get_test_count()}: Loading of pickled 'kde' models..."
        )
        cva = Constava(verbose=self.verbosity, model_load=self.TEST_MODELDUMP0)
        cva.load_csmodel(pickled_csmodel=self.TEST_MODELDUMP0)
        self.assertIsInstance(
            cva._csmodel, ConfStateModelKDE, "Failed to fit conformational state model."
        )

    def test_2a_KDEModelInference(self):
        """Test inference from pre-fitted KDE model"""
        logger.warning(
            f"TEST #{self.get_test_count()}: Inference from 'kde' conformational state model..."
        )
        input_files = glob.glob(f"{self.CONSTAVA_TEST_DATA_DIR}/xvg/ramaPhiPsi*.xvg")
        expected_result = f"{self.CONSTAVA_TEST_DATA_DIR}/xvg/result_kde.csv"
        output_file = tempfile.NamedTemporaryFile(
            prefix="kde.", suffix=".csv", dir=self.TEST_TEMPDIR
        )
        cva = Constava(
            input_files=input_files,
            output_file=output_file.name,
            model_type="kde",
            model_load=self.TEST_MODELDUMP0,
            window=[1, 3, 7, 23],
            bootstrap=[3, 7, 23],
            window_series=[1, 7],
            seed=42,
            verbose=self.verbosity,
            input_degrees=True,
        )
        cva.run()
        filecmp.clear_cache()
        self.assertTrue(filecmp.cmp(output_file.name, expected_result))

    def test_3a_KDEModelInference(self):
        """Test inference from pre-fitted KDE model"""
        logger.warning(
            f"TEST #{self.get_test_count()}: Inference from 'kde' conformational state model in JSON format..."
        )
        input_files = glob.glob(f"{self.CONSTAVA_TEST_DATA_DIR}/xvg/ramaPhiPsi*.xvg")
        expected_result = f"{self.CONSTAVA_TEST_DATA_DIR}/xvg/result_kde.json"
        output_file = tempfile.NamedTemporaryFile(
            prefix="kde.", suffix=".json", dir=self.TEST_TEMPDIR
        )
        cva = Constava(
            input_files=input_files,
            output_file=output_file.name,
            model_type="kde",
            model_load=self.TEST_MODELDUMP0,
            window=[1, 3, 7, 23],
            bootstrap=[3, 7, 23],
            window_series=[1, 7],
            seed=42,
            verbose=self.verbosity,
            input_degrees=True,
        )
        cva.run()

        with open(expected_result, "r", encoding="utf-8") as f_exp:
            expected_obj = json.load(f_exp)

        with open(
            os.path.join(self.TEST_TEMPDIR, "expected_test_3a_KDEModelInference.json"),
            "w",
            encoding="utf-8",
        ) as f:
            json.dump(expected_obj["results"], f)

        with open(output_file.name, "r", encoding="utf-8") as f_act:
            actual_obj = json.load(f_act)

        with open(
            os.path.join(self.TEST_TEMPDIR, "actual_test_3a_KDEModelInference.json"),
            "w",
            encoding="utf-8",
        ) as f:
            json.dump(actual_obj["results"], f)

        filecmp.clear_cache()
        self.assertTrue(
            filecmp.cmp(
                os.path.join(
                    self.TEST_TEMPDIR, "expected_test_3a_KDEModelInference.json"
                ),
                os.path.join(
                    self.TEST_TEMPDIR, "actual_test_3a_KDEModelInference.json"
                ),
            )
        )

    def test_0b_GridModelFitting(self):
        """Test fitting of a grid-interpolation model, dumps model for further tests"""
        logger.warning(f"TEST #{self.get_test_count()}: Fitting of 'grid' models...")
        cva = Constava(verbose=self.verbosity)
        csmodel = cva.fit_csmodel(
            model_type="grid", kde_bandwidth=0.42, grid_points=145
        )
        self.assertIsInstance(
            csmodel, ConfStateModelGrid, "Failed to fit conformational state model."
        )

        logger.warning(f"TEST #{self.get_test_count()}: Caching of 'grid' models...")
        cva.fit_csmodel(model_type="grid", kde_bandwidth=0.42, grid_points=145)
        self.assertIs(
            csmodel,
            cva._csmodel,
            "Failed to reuse preloaded conformational state model.",
        )

        logger.warning(f"TEST #{self.get_test_count()}: Fitting a new 'grid' model...")
        cva.fit_csmodel(model_type="grid", kde_bandwidth=0.1, grid_points=3601)
        self.assertIsNot(
            csmodel,
            cva._csmodel,
            "Failed to update conformational state model after parameter change.",
        )

        logger.warning(
            f"TEST #{self.get_test_count()}: Storing of pickled 'grid' models..."
        )
        logger.warning(f"... Dumping model in: {self.TEST_MODELDUMP1}")
        cva._csmodel.dump_pickle(self.TEST_MODELDUMP1)
        self.assertTrue(
            os.path.isfile(self.TEST_MODELDUMP1),
            "Failed to conformational state dump model.",
        )

    def test_1b_GridModelLoading(self):
        """Test loading of pre-fitted grid-inference model"""
        logger.warning(
            f"TEST #{self.get_test_count()}: Loading of pickled 'grid' models..."
        )
        cva = Constava(verbose=self.verbosity, model_load=self.TEST_MODELDUMP1)
        cva.load_csmodel(pickled_csmodel=self.TEST_MODELDUMP1)
        self.assertIsInstance(
            cva._csmodel,
            ConfStateModelGrid,
            "Failed to fit conformational state model.",
        )

    def test_2b_GridModelInference(self):
        """Test inference from pre-fitted grid-inference model"""
        logger.warning(
            f"TEST #{self.get_test_count()}: Inference from 'grid' conformational state model..."
        )
        input_files = f"{self.CONSTAVA_TEST_DATA_DIR}/csv/dihedrals.csv"
        expected_result = (
            f"{self.CONSTAVA_TEST_DATA_DIR}/csv/result_grid.bandwidth.42.csv"
        )
        output_file = tempfile.NamedTemporaryFile(
            prefix="grid.", suffix=".csv", dir=self.TEST_TEMPDIR
        )
        c = Constava(
            input_files=input_files,
            output_file=output_file.name,
            model_type="grid",
            model_load=self.TEST_MODELDUMP1,
            window=[1, 3, 7, 23],
            bootstrap=[3, 7, 23],
            window_series=[1, 7],
            seed=42,
            verbose=self.verbosity,
            input_degrees=False,
        )
        c.run()

        filecmp.clear_cache()
        self.assertTrue(filecmp.cmp(output_file.name, expected_result))

    def test_3b_GridModelInference(self):
        """Test inference from pre-fitted grid-inference model"""
        logger.warning(
            f"TEST #{self.get_test_count()}: Inference from 'grid' conformational state model in JSON format..."
        )
        input_files = f"{self.CONSTAVA_TEST_DATA_DIR}/csv/dihedrals.csv"
        expected_result = f"{self.CONSTAVA_TEST_DATA_DIR}/csv/result_grid.json"
        output_file = tempfile.NamedTemporaryFile(
            prefix="grid.", suffix=".json", dir=self.TEST_TEMPDIR
        )

        c = Constava(
            input_files=input_files,
            output_file=output_file.name,
            model_type="grid",
            model_load=self.TEST_MODELDUMP1,
            window=[1, 3, 7, 23],
            bootstrap=[3, 7, 23],
            window_series=[1, 7],
            seed=42,
            verbose=self.verbosity,
            input_degrees=False,
        )
        c.run()

        with open(expected_result, "r", encoding="utf-8") as f_exp:
            expected_obj = json.load(f_exp)

        with open(
            os.path.join(self.TEST_TEMPDIR, "expected_test_3b_GridModelInference.json"),
            "w",
            encoding="utf-8",
        ) as f:
            json.dump(expected_obj["results"], f)

        with open(output_file.name, "r", encoding="utf-8") as f_act:
            actual_obj = json.load(f_act)

        with open(
            os.path.join(self.TEST_TEMPDIR, "actual_test_3b_GridModelInference.json"),
            "w",
            encoding="utf-8",
        ) as f:
            json.dump(actual_obj["results"], f)

        filecmp.clear_cache()
        self.assertTrue(
            filecmp.cmp(
                os.path.join(
                    self.TEST_TEMPDIR, "expected_test_3b_GridModelInference.json"
                ),
                os.path.join(
                    self.TEST_TEMPDIR, "actual_test_3b_GridModelInference.json"
                ),
            )
        )


def run_unittest(verbose: int = 0):
    """Run the unit tests in this module"""

    if verbose == 0:
        logger.setLevel(logging.WARNING)
    elif verbose == 1:
        logger.setLevel(logging.INFO)
    else:
        logger.setLevel(logging.DEBUG)

    suite = unittest.TestSuite(
        (
            TestWrapper("test_0a_KDEModelFitting", verbose),
            TestWrapper("test_1a_KDEModelLoading", verbose),
            TestWrapper("test_2a_KDEModelInference", verbose),
            TestWrapper("test_3a_KDEModelInference", verbose),
            TestWrapper("test_0b_GridModelFitting", verbose),
            TestWrapper("test_1b_GridModelLoading", verbose),
            TestWrapper("test_2b_GridModelInference", verbose),
            TestWrapper("test_3b_GridModelInference", verbose),
        )
    )
    runner = unittest.TextTestRunner(verbosity=verbose)
    runner.run(suite)


if __name__ == "__main__":
    run_unittest()
