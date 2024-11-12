import unittest
from vsutils.ADMET.filters import (
    calculate_ro5_properties,
    calculate_pfizer_rule,
    calculate_gsk_rule,
    calculate_goldentriangle_rule,
    # calculate_qed,
    calculate_sascore,
    calculate_fsp3,
    pains_filter,
    # calculate_mce18,
    calculate_npscore,
    alarm_nmr_filter,
    bms_filter,
    chelator_filter,
)


class TestChemicalProperties(unittest.TestCase):
    def test_calculate_ro5_properties(self):
        self.assertTrue(calculate_ro5_properties("CCO"))
        self.assertFalse(
            calculate_ro5_properties(
                "CC1(COC(C(C1NC)O)"
                + "OC2C(CC3C(C2O)OC4C(CC=C(O4)C=NC5CC(C(C(C5OC6C(CC=C(O6)C=N3)N)O)OC7C(C(C(CO7)(C)O)NC)O)N)N)N)O"
            )
        )

    def test_calculate_pfizer_rule(self):
        # self.assertTrue(calculate_pfizer_rule("CCCCCC(=O)O"))
        self.assertFalse(calculate_pfizer_rule("CCO"))

    def test_calculate_gsk_rule(self):
        self.assertTrue(calculate_gsk_rule("CCO"))
        self.assertFalse(calculate_gsk_rule("C1CCCCC1CCCCCCCCC(=O)O"))

    def test_calculate_goldentriangle_rule(self):
        self.assertFalse(calculate_goldentriangle_rule("CCO"))
        self.assertTrue(calculate_goldentriangle_rule("C1CCCCC1CCCCCCCCC(=O)O"))

    # def test_calculate_qed(self):
    #     qed, is_attractive = calculate_qed('CC(=O)OC1=CC=CC=C1C(=O)O')
    #     self.assertGreater(qed, 0.67)
    #     self.assertTrue(is_attractive)

    def test_calculate_sascore(self):
        sascore, is_easy_to_synthesize = calculate_sascore("CCO")
        self.assertLessEqual(sascore, 6)
        self.assertTrue(is_easy_to_synthesize)

    def test_calculate_fsp3(self):
        fsp3, is_suitable = calculate_fsp3("CCC(C)(C)C")
        self.assertGreaterEqual(fsp3, 0.42)
        self.assertTrue(is_suitable)

    def test_pains_filter(self):
        is_accepted, _, _ = pains_filter("c1ccccc1")
        self.assertTrue(is_accepted)

    # def test_calculate_mce18(self):
    #     mce18, is_interesting = calculate_mce18('c1ccccc1O')
    #     self.assertGreaterEqual(mce18, 45)
    #     self.assertTrue(is_interesting)

    def test_calculate_npscore(self):
        np_score = calculate_npscore("CCO")
        self.assertGreater(np_score, 0)

    def test_alarm_nmr_filter(self):
        is_accepted, _, _ = alarm_nmr_filter("CCO")
        self.assertTrue(is_accepted)

    def test_bms_filter(self):
        is_accepted, _, _ = bms_filter("CCO")
        self.assertTrue(is_accepted)

    def test_chelator_filter(self):
        is_accepted, _, _ = chelator_filter("CCO")
        self.assertTrue(is_accepted)


if __name__ == "__main__":
    unittest.main()
