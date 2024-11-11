import os
import sys
from typing import Tuple, List
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import RDConfig
from scopy.ScoFH import fh_filter

sascore_path = os.path.join(RDConfig.RDContribDir, "SA_Score")
if sascore_path not in sys.path:
    sys.path.append(sascore_path)
import sascorer

npscore_path = os.path.join(RDConfig.RDContribDir, "NP_Score")
if npscore_path not in sys.path:
    sys.path.append(npscore_path)
import npscorer


def calculate_ro5_properties(smiles: str, fulfill: int = 4) -> bool:
    """
    Determines if a molecule represented by a SMILES string fulfills the Lipinski's rule of five.

    Args:
    - smiles (str): SMILES representation of the molecule.
    - fulfill (int): The number of rules that must be met to consider the molecule compliant.

    Returns:
    - bool: True if the molecule meets the specified number of rules, False otherwise.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return False
        molecular_weight = Descriptors.ExactMolWt(molecule)
        n_hba = Descriptors.NumHAcceptors(molecule)
        n_hbd = Descriptors.NumHDonors(molecule)
        logp = Descriptors.MolLogP(molecule)

        conditions = [molecular_weight <= 500, n_hba <= 10, n_hbd <= 5, logp <= 5]
        return sum(conditions) >= fulfill
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return False


def calculate_pfizer_rule(smiles: str) -> bool:
    """
    Determines if a molecule represented by a SMILES string fulfills the Pfizer Rule.

    Args:
    - smiles (str): SMILES representation of the molecule.

    Returns:
    - bool: True if the molecule fulfills the Pfizer Rule, False otherwise.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return False
        logp = Descriptors.MolLogP(molecule)
        tpsa = Descriptors.TPSA(molecule)

        conditions = [logp > 3, tpsa < 75]
        return sum(conditions) == len(conditions)
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return False


def calculate_gsk_rule(smiles: str) -> bool:
    """
    Evaluates whether a molecule, defined by its SMILES string, satisfies the GSK Rule,
    which is defined by specific thresholds for molecular weight and LogP.

    Args:
    - smiles (str): SMILES representation of the molecule.

    Returns:
    - bool: True if the molecule fulfills both criteria of the GSK Rule, False otherwise.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return False
        molecular_weight = Descriptors.ExactMolWt(molecule)
        logp = Descriptors.MolLogP(molecule)

        conditions = [molecular_weight <= 400, logp <= 4]
        return sum(conditions) == len(conditions)
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return False


def calculate_goldentriangle_rule(smiles: str) -> bool:
    """
    Determines if a molecule represented by a SMILES string fulfills the GoldenTriangle Rule,
    defined by specific ranges for molecular weight and LogP.

    Args:
    - smiles (str): SMILES representation of the molecule.

    Returns:
    - bool: True if the molecule fulfills both criteria of the GoldenTriangle Rule, False otherwise.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return False
        molecular_weight = Descriptors.ExactMolWt(molecule)
        logp = Descriptors.MolLogP(molecule)

        conditions = [200 <= molecular_weight <= 450, -2 <= logp <= 5]
        return sum(conditions) == len(conditions)
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return False


def calculate_qed(smiles: str) -> Tuple[float, bool]:
    """
    Calculate the Quantitative Estimate of Drug-likeness (QED) and determine if the molecule is
    deemed 'attractive' based on the QED score.

    Args:
    - smiles (str): SMILES string representing the molecule.

    Returns:
    - Tuple[float, bool]: A tuple containing the QED value and a boolean indicating if the QED is above 0.67,
    deemed 'attractive'.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return 0.0, False
        qed = Chem.QED.qed(molecule)
        qed_excellent = qed > 0.67
        return qed, qed_excellent
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return 0.0, False


def calculate_sascore(smiles: str) -> Tuple[float, bool]:
    """
    Calculate the synthetic accessibility score (SAscore) and determine if the molecule is easy
    to synthesize based on this score.

    Args:
    - smiles (str): SMILES string representing the molecule.

    Returns:
    - Tuple[float, bool]: A tuple containing the SAscore and a boolean indicating if the score is 6 or less,
    which means better synthetic accessibility.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return (
                10.0,
                False,
            )  # Assuming 10 as a default high SAscore for inaccessible synthesis
        SAscore = sascorer.calculateScore(molecule)
        SAscore_excellent = SAscore <= 6
        return SAscore, SAscore_excellent
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return 10.0, False


def calculate_fsp3(smiles: str) -> Tuple[float, bool]:
    """
    Calculate the fraction of sp3 hybridized carbons (Fsp3) and determine if the molecule has
    a suitable value for Fsp3.

    Args:
    - smiles (str): SMILES string representing the molecule.

    Returns:
    - Tuple[float, bool]: A tuple containing the Fsp3 ratio and a boolean indicating if this ratio
    is 0.42 or higher, considered suitable.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return 0.0, False
        fsp3 = Chem.rdMolDescriptors.CalcFractionCSP3(molecule)
        fsp3_excellent = fsp3 >= 0.42
        return fsp3, fsp3_excellent
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return 0.0, False


def pains_filter(smiles: str) -> Tuple[bool, List[str], List[int]]:
    """
    Applies the PAINS filter to a molecule given its SMILES string and returns the filter result
    and details.

    Args:
    - smiles (str): SMILES representation of the molecule.

    Returns:
    - Tuple[bool, List[str], List[int]]: A tuple containing a boolean indicating if
    the molecule passed the PAINS filter, a list of matched PAINS names, and a list
    of matched atom indices.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return False, [], []
        pains = fh_filter.Check_PAINS(molecule, detail=True)
        return (
            pains["Disposed"] == "Accepted",
            pains["MatchedNames"],
            pains["MatchedAtoms"],
        )
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return False, [], []


def calculate_mce18(smiles: str) -> Tuple[float, bool]:
    """
    Calculates the MCE-18 score for a molecule to assess its complexity and structural features.

    Args:
    - smiles (str): SMILES representation of the molecule.

    Returns:
    - Tuple[float, bool]: A tuple containing the MCE-18 score and a boolean indicating if the score
    is excellent (>= 45).
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return 0.0, False
        from rdkit.Chem import rdMolDescriptors

        # Extract properties
        properties = {
            "AR": rdMolDescriptors.CalcNumAromaticRings(molecule) > 0,
            "NAR": rdMolDescriptors.CalcNumAliphaticRings(molecule) > 0,
            "CHIRAL": len(
                Chem.FindMolChiralCenters(molecule, force=True, includeUnassigned=True)
            )
            > 0,
            "SPIRO": rdMolDescriptors.CalcNumSpiroAtoms(molecule),
            "SP3": rdMolDescriptors.CalcFractionCSP3(molecule),
        }
        # Additional calculations for MCE-18
        cyclic = [
            1
            for atom in molecule.GetAtoms()
            if atom.GetAtomicNum() == 6
            and atom.IsInRing()
            and atom.GetHybridization() == Chem.HybridizationType.SP3
        ]
        acyclic = [
            1
            for atom in molecule.GetAtoms()
            if atom.GetAtomicNum() == 6
            and not atom.IsInRing()
            and atom.GetHybridization() == Chem.HybridizationType.SP3
        ]
        CYC, ACYC = sum(cyclic), sum(acyclic)
        total_C = sum([1 for atom in molecule.GetAtoms() if atom.GetAtomicNum() == 6])
        CYC_ratio = CYC / total_C if total_C else 0
        ACYC_ratio = ACYC / total_C if total_C else 0
        Q1 = (
            3
            - 2 * molecule.GetNumAtoms()
            + sum(atom.GetDegree() ** 2 for atom in molecule.GetAtoms()) / 2.0
        )
        mce18 = (
            properties["AR"]
            + properties["NAR"]
            + properties["CHIRAL"]
            + properties["SPIRO"]
            + (properties["SP3"] + CYC_ratio - ACYC_ratio) / (1 + properties["SP3"])
        ) * Q1
        return mce18, mce18 >= 45
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return 0.0, False


def calculate_npscore(smiles: str) -> float:
    """
    Calculates the Natural Product likeness score (NPscore) for a molecule based on its SMILES string.

    Args:
    - smiles (str): SMILES representation of the molecule.

    Returns:
    - float: NPscore of the molecule.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return 0.0
        fscore = npscorer.readNPModel()
        return npscorer.scoreMol(molecule, fscore)
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return 0.0


def alarm_nmr_filter(smiles: str) -> Tuple[bool, List[str], List[int]]:
    """
    Applies the ALARM NMR filter to assess if a molecule violates known NMR problem structures.

    Args:
    - smiles (str): SMILES representation of the molecule.

    Returns:
    - Tuple[bool, List[str], List[int]]: A tuple containing a boolean indicating if the molecule
    passed the ALARM NMR filter, a list of matched problematic structures, and a list of
    matched atom indices.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return False, [], []
        alarm_nmr = fh_filter.Check_Alarm_NMR(molecule, detail=True)
        return (
            alarm_nmr["Disposed"] == "Accepted",
            alarm_nmr["MatchedNames"],
            alarm_nmr["MatchedAtoms"],
        )
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return False, [], []


def bms_filter(smiles: str) -> Tuple[bool, List[str], List[int]]:
    """
    Applies the BMS (Bristol-Myers Squibb) filter to a molecule represented by a SMILES string to identify
    problematic substructures.

    Args:
    - smiles (str): SMILES string representing the molecule.

    Returns:
    - Tuple[bool, List[str], List[int]]: A tuple containing whether the molecule passes the BMS filter,
      a list of substructure names causing violations if any, and their corresponding atom indices.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return False, [], []
        bms = fh_filter.Check_BMS(molecule, detail=True)
        bms_accepted = bms["Disposed"] == "Accepted"
        return bms_accepted, bms["MatchedNames"], bms["MatchedAtoms"]
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return False, [], []


def chelator_filter(smiles: str) -> Tuple[bool, List[str], List[int]]:
    """
    Applies the Chelator filter to a molecule given its SMILES string to assess whether it has properties that
    make it act as a chelating agent.

    Args:
    - smiles (str): SMILES representation of the molecule.

    Returns:
    - Tuple[bool, List[str], List[int]]: A tuple indicating if the molecule is accepted by the Chelator filter,
      along with lists of names and atom indices of any matched chelating substructures.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        if molecule is None:
            return False, [], []
        chelator = fh_filter.Check_Chelating(molecule, detail=True)
        chelator_accepted = chelator["Disposed"] == "Accepted"
        return chelator_accepted, chelator["MatchedNames"], chelator["MatchedAtoms"]
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {str(e)}")
        return False, [], []
