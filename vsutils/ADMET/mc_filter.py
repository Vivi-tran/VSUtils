from vsutils.ADMET.filters import (
    calculate_ro5_properties,
    calculate_pfizer_rule,
    calculate_gsk_rule,
    calculate_goldentriangle_rule,
    calculate_qed,
    calculate_sascore,
    calculate_fsp3,
    pains_filter,
    calculate_mce18,
    calculate_npscore,
    alarm_nmr_filter,
    bms_filter,
    chelator_filter,
)
from joblib import Parallel, delayed
import pandas as pd


class MCFilter:
    def __init__(self, n_jobs=1, verbose=0) -> None:
        self.n_jobs = n_jobs
        self.verbose = verbose

    @staticmethod
    def calculate_des(smiles: str) -> pd.Series:
        """
        Calculate all rules for a single molecule SMILES string.

        Args:
        - smiles (str): SMILES representation of the molecule.

        Returns:
        - pandas.Series: Series containing results of various filters and descriptors.
        """
        ro5_rule = calculate_ro5_properties(smiles)
        pfizer_rule = calculate_pfizer_rule(smiles)
        gsk_rule = calculate_gsk_rule(smiles)
        goldentriangle_rule = calculate_goldentriangle_rule(smiles)
        qed = calculate_qed(smiles)
        sascore = calculate_sascore(smiles)
        fsp3 = calculate_fsp3(smiles)
        mce18 = calculate_mce18(smiles)
        npscore = calculate_npscore(smiles)
        pains = pains_filter(smiles)
        alarmnmr = alarm_nmr_filter(smiles)
        bms = bms_filter(smiles)
        chelator = chelator_filter(smiles)

        return pd.Series(
            [
                ro5_rule,
                pfizer_rule,
                gsk_rule,
                goldentriangle_rule,
                qed[0],
                qed[1],
                sascore[0],
                sascore[1],
                fsp3[0],
                fsp3[1],
                mce18[0],
                mce18[1],
                npscore,
                pains[0],
                pains[1],
                pains[2],
                alarmnmr[0],
                alarmnmr[1],
                alarmnmr[2],
                bms[0],
                bms[1],
                bms[2],
                chelator[0],
                chelator[1],
                chelator[2],
            ],
            index=[
                "ro5_rule",
                "pfizer_rule",
                "gsk_rule",
                "goldentriangle_rule",
                "qed",
                "qed_excellent",
                "sascore",
                "sascore_excellent",
                "fsp3",
                "fsp3_excellent",
                "mce18",
                "mce18_excellent",
                "npscore",
                "pains_accepted",
                "pains_matched_names",
                "pains_matched_atoms",
                "alarmnmr_accepted",
                "alarmnmr_matched_names",
                "alarmnmr_matched_atoms",
                "bms_accepted",
                "bms_matched_names",
                "bms_matched_atoms",
                "chelator_accepted",
                "chelator_matched_names",
                "chelator_matched_atoms",
            ],
        )

    def process_dataframe(self, df: pd.DataFrame, smiles_column: str) -> pd.DataFrame:
        """
        Processes a DataFrame in parallel to calculate descriptors and filter properties for each SMILES string.

        Args:
        - df (pd.DataFrame): DataFrame containing a column of SMILES strings.
        - smiles_column (str): The name of the column containing the SMILES strings.

        Returns:
        - pd.DataFrame: DataFrame with the results of filter application for each SMILES string.
        """
        smiles_list = df[smiles_column].tolist()
        results = Parallel(n_jobs=self.n_jobs, verbose=self.verbose)(
            delayed(self.calculate_des)(smiles) for smiles in smiles_list
        )

        return pd.DataFrame(results, index=df.index)
