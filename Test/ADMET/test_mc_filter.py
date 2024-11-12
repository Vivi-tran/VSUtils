import unittest
import pandas as pd
from vsutils.ADMET.mc_filter import (
    MCFilter,
)


class TestMCFilter(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.filter = MCFilter()

        # Test data
        cls.data = pd.DataFrame(
            {
                "smiles": [
                    "CN1C=C(C2=CC=CC=C21)C3=NC(=NC=C3)NC4=C(C=C(C(=C4)NC(=O)C=C)N(C)CCN(C)C)OC",
                    "CCO",
                    (
                        "CC1(COC(C(C1NC)O)OC2C(CC3C(C2O)OC4C(CC=C(O4)C=NC5CC(C(C(C5OC6C(CC=C(O6)C=N3)"
                        + "N)O)OC7C(C(C(CO7)(C)O)NC)O)N)N)O"
                    ),
                ]
            }
        )

    def test_process_dataframe(self):
        results_df = self.filter.process_dataframe(self.data, "smiles")
        self.assertEqual(results_df.shape[0], 3)
        self.assertTrue("ro5_rule" in results_df.columns)
        self.assertTrue(results_df.loc[0, "ro5_rule"])
        self.assertFalse(results_df.loc[0, "pfizer_rule"])

    def test_calculate_des(self):

        series_result = self.filter.calculate_des(
            "CN1C=C(C2=CC=CC=C21)C3=NC(=NC=C3)NC4=C(C=C(C(=C4)NC(=O)C=C)N(C)CCN(C)C)OC"
        )
        self.assertIsInstance(series_result, pd.Series)
        self.assertTrue(series_result["ro5_rule"])
        self.assertFalse(series_result["pfizer_rule"])


if __name__ == "__main__":
    unittest.main()
