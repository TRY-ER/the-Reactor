import pytest
from benchmarking.metrics.SMILES import SMILESMetric
from benchmarking.metrics.base import CompareObj, BaseReference

# Each entry: (list of equivalent SMILES, list of alternative representations)
smiles_variations = [
    # Ethanol
    (["CCO"], ["C(C)O", "OCC"]),
    # Benzene
    (["c1ccccc1"], ["C1=CC=CC=C1"]),
    # Cyclohexane
    (["C1CCCCC1"], ["C1CCCCC1"]),
    # Acetic acid
    (["CC(=O)O"], ["OC(=O)C"]),
    # L-alanine (using canonical form in reference)
    (["C[C@H](N)C(=O)O"], ["N[C@@H](C)C(=O)O"]),
    # Ammonium ion
    (["[NH4+]"], ["[NH4+]"]),
    # Sodium chloride (using canonical form in reference)
    (["[Cl-].[Na+]"], ["[Na+].[Cl-]"]),
    # Acetylene
    (["C#C"], ["C#C"]),
    # Pyridine
    (["c1ccncc1"], ["n1ccccc1"]),
]

def test_smiles_metric_equivalence():
    metric = SMILESMetric()
    for base, variations in smiles_variations:
        base_ref = BaseReference(value_list=[base])
        for var in variations:
            compare_obj = CompareObj(value_list=[[var]])
            result = metric.compute(compare_obj, base_ref)
            assert result.correct == 1.0, f"Failed for base: {base[0]}, variation: {var}, got correct={result.correct}"
            assert result.valid == 1.0, f"Failed for base: {base[0]}, variation: {var}, got valid={result.valid}"
            assert result.found == 1.0, f"Failed for base: {base[0]}, variation: {var}, got found={result.found}"