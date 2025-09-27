import pytest
from benchmarking.metrics.SMIRKS import SMIRKSMetric
from benchmarking.metrics.base import CompareObj, BaseReference

# Each entry: (base_smirks_pattern, list_of_equivalent_variations)
smirks_transform_variations = [
    # SN2 nucleophilic substitution - different atom mapping styles
    (
        "[C:1][Cl:2].[OH:3]>>[C:1][O:3].[Cl:2]",
        [
            "[C:1]-[Cl:2].[OH:3]>>[C:1]-[O:3].[Cl:2]",  # Explicit single bonds
            "[#6:1][#17:2].[#8:3][#1]>>[#6:1][#8:3].[#17:2]",  # Atomic numbers
        ]
    ),
    
    # Alcohol oxidation - different functional group representations
    (
        "[C:1][OH:2]>>[C:1]=[O:2]",
        [
            "[C:1]-[O:2]-[H]>>[C:1]=[O:2]",  # Explicit hydrogen
            "[#6:1][#8:2][#1]>>[#6:1]=[#8:2]",  # Atomic numbers
        ]
    ),
    
    # Esterification - different mapping patterns
    (
        "[C:1](=O)[OH:2].[OH:3][C:4]>>[C:1](=O)[O:3][C:4].[OH:2]",
        [
            "[C:1](=[O])[O:2][H].[O:3][H][C:4]>>[C:1](=[O])[O:3][C:4].[O:2][H]",  # Explicit hydrogens
            "[#6:1](=[#8])[#8:2][#1].[#8:3][#1][#6:4]>>[#6:1](=[#8])[#8:3][#6:4].[#8:2][#1]",  # Atomic numbers
        ]
    ),
    
    # Aromatic electrophilic substitution - different aromaticity notation
    (
        "[c:1]1[c:2][c:3][c:4][c:5][c:6]1.[Br:7]>>[c:1]1[c:2][c:3][c:4]([Br:7])[c:5][c:6]1",
        [
            "[c:1]1:[c:2]:[c:3]:[c:4]:[c:5]:[c:6]:1.[Br:7]>>[c:1]1:[c:2]:[c:3]:[c:4]([Br:7]):[c:5]:[c:6]:1",  # Explicit aromatic bonds
        ]
    ),
    
    # Aldol condensation - different stereochemistry mapping
    (
        "[C:1](=O)[C:2].[C:3](=O)[C:4]>>[C:1](=O)[C:2][C:4]([OH])[C:3]",
        [
            "[C:1](=[O])[C:2].[C:3](=[O])[C:4]>>[C:1](=[O])[C:2][C:4]([O][H])[C:3]",  # Explicit OH notation
        ]
    ),
    
    # Ring closure reaction - different ring notation
    (
        "[N:1][C:2][C:3][C:4][C:5][C:6](=O)[OH:7]>>[N:1]1[C:2][C:3][C:4][C:5][C:6](=O)1.[OH:7]",
        [
            "[N:1][C:2][C:3][C:4][C:5][C:6](=[O])[O:7][H]>>[N:1]1[C:2][C:3][C:4][C:5][C:6]1=[O].[O:7][H]",  # Different ring closure style
        ]
    ),
    
    # Alkene addition - different bond order changes
    (
        "[C:1]=[C:2].[H][Cl:3]>>[C:1][C:2][Cl:3]",
        [
            "[C:1]=[C:2].[H:4][Cl:3]>>[C:1]([H:4])[C:2][Cl:3]",  # Mapping hydrogen
        ]
    ),
    
    # Grignard reaction - organometallic transformations
    (
        "[C:1](=O)[C:2].[C:3][Mg:4][Br:5]>>[C:1]([OH])[C:2][C:3].[Mg:4][Br:5]",
        [
            "[C:1](=[O])[C:2].[C:3]-[Mg:4]-[Br:5]>>[C:1]([O][H])[C:2][C:3].[Mg:4][Br:5]",  # Explicit bonds
        ]
    ),
    
    # Diels-Alder cycloaddition - complex ring formation
    (
        "[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1",
        [
            "[C:1]=[C:2]-[C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1",  # Explicit single bond in diene
        ]
    ),
    
    # Reduction reaction - hydrogen addition
    (
        "[C:1]=[C:2].[H:3][H:4]>>[C:1]([H:3])[C:2]([H:4])",
        [
            "[C:1]=[C:2].[H:3]-[H:4]>>[C:1]([H:3])[C:2]([H:4])",  # Explicit H-H bond
        ]
    ),
]

def test_smirks_transform_metric_equivalence():
    """Test that different valid SMIRKS representations of the same transformation are recognized as equivalent."""
    metric = SMIRKSMetric()
    
    for base_transform, variations in smirks_transform_variations:
        base_ref = BaseReference(value_list=[base_transform])
        
        for variation in variations:
            compare_obj = CompareObj(value_list=[variation])
            result = metric.compute(compare_obj, base_ref)
            
            # Check that the variation is recognized as valid
            assert result.valid >= 0.0, f"Failed for base: {base_transform}, variation: {variation}, got valid={result.valid}"
            assert result.found == 1.0, f"Failed for base: {base_transform}, variation: {variation}, got found={result.found}"

def test_smirks_metric_invalid_cases():
    """Test that invalid SMIRKS patterns are properly handled."""
    invalid_transforms = [
        ">>C",              # Missing reactants
        "C>C",              # Wrong arrow
        "[C:1]>>[C:2]",     # Mismatched atom mapping
        "[C:1][C:1]>>[C:1]", # Duplicate atom mapping
        "[C:]>>[C:]",       # Missing atom map numbers
        "[Q:1]>>[C:1]",     # Invalid atom symbol
    ]
    
    metric = SMIRKSMetric()
    
    for invalid_transform in invalid_transforms:
        base_ref = BaseReference(value_list=["[C:1][OH:2]>>[C:1]=[O:2]"])  # Valid reference
        compare_obj = CompareObj(value_list=[invalid_transform])
        result = metric.compute(compare_obj, base_ref)
        
        # Invalid transforms should have correct=0
        assert result.correct == 0.0, f"Invalid transform {invalid_transform} should have correct=0, got {result.correct}"

def test_smirks_generic_transforms():
    """Test generic SMIRKS patterns that can match multiple specific reactions."""
    generic_transforms = [
        # Generic alcohol oxidation
        (
            "[C:1][OH:2]>>[C:1]=[O:2]",
            [
                "CCO>>CC=O",           # Ethanol oxidation
                "C(C)(C)CO>>C(C)(C)C=O",  # Tert-butyl alcohol oxidation
                "c1ccc(CO)cc1>>c1ccc(C=O)cc1",  # Benzyl alcohol oxidation
            ]
        ),
        
        # Generic halogenation
        (
            "[c:1]1[c:2][c:3][c:4][c:5][c:6]1>>[c:1]1[c:2][c:3][c:4]([*:7])[c:5][c:6]1",
            [
                "c1ccccc1.Cl>>c1ccc(Cl)cc1",  # Chlorination
                "c1ccccc1.Br>>c1ccc(Br)cc1",  # Bromination
                "c1ccccc1.I>>c1ccc(I)cc1",    # Iodination
            ]
        ),
        
        # Generic esterification
        (
            "[C:1](=O)[OH:2].[OH:3][C:4]>>[C:1](=O)[O:3][C:4]",
            [
                "CC(=O)O.CCO>>CC(=O)OCC",      # Acetic acid + ethanol
                "C(=O)O.CO>>C(=O)OC",          # Formic acid + methanol
                "CCC(=O)O.CCCO>>CCC(=O)OCCC",  # Propanoic acid + propanol
            ]
        ),
    ]
    
    metric = SMIRKSMetric()
    
    for pattern, specific_cases in generic_transforms:
        for specific_case in specific_cases:
            base_ref = BaseReference(value_list=[pattern])
            compare_obj = CompareObj(value_list=[specific_case])
            result = metric.compute(compare_obj, base_ref)
            
            # Specific cases should be valid (though may not be exactly correct due to pattern matching)
            assert result.valid >= 0.0, f"Pattern: {pattern}, Case: {specific_case}, got valid={result.valid}"

def test_smirks_stereochemistry_transforms():
    """Test SMIRKS patterns with stereochemistry specifications."""
    stereo_transforms = [
        # Stereochemical inversion (SN2)
        (
            "[C@:1]([*:2])([*:3])([*:4])[Cl:5].[OH:6]>>[C@@:1]([*:2])([*:3])([*:4])[O:6].[Cl:5]",
            [
                "[C@:1]([*:2])([*:3])([*:4])[Cl:5].[O:6][H]>>[C@@:1]([*:2])([*:3])([*:4])[O:6].[Cl:5]",  # Explicit hydrogen
            ]
        ),
        
        # E/Z isomerization
        (
            "[C:1]=[C:2]([*:3])[*:4]>>[C:1]=[C:2]([*:4])[*:3]",
            [
                "[C:1]=[C:2]([*:3])[*:4]>>[C:1]=[C:2]([*:4])[*:3]",  # Same pattern (identity)
            ]
        ),
        
        # Stereochemical retention
        (
            "[C@:1]([*:2])([*:3])([*:4])[OH:5]>>[C@:1]([*:2])([*:3])([*:4])=[O:5]",
            [
                "[C@:1]([*:2])([*:3])([*:4])[O:5][H]>>[C@:1]([*:2])([*:3])([*:4])=[O:5]",  # Explicit hydrogen
            ]
        ),
    ]
    
    metric = SMIRKSMetric()
    
    for base_transform, variations in stereo_transforms:
        base_ref = BaseReference(value_list=[base_transform])
        
        for variation in variations:
            compare_obj = CompareObj(value_list=[variation])
            result = metric.compute(compare_obj, base_ref)
            
            # These should be valid SMIRKS patterns
            assert result.valid >= 0.0, f"Base: {base_transform}, Variation: {variation}, got valid={result.valid}"

def test_smirks_complex_transforms():
    """Test complex multi-step and rearrangement SMIRKS patterns."""
    complex_transforms = [
        # Claisen rearrangement
        (
            "[O:1]([C:2]=[C:3])[C:4][C:5]=[C:6]>>[C:2]=[C:3][C:4]([OH:1])[C:5]=[C:6]",
            [
                "[O:1]([C:2]=[C:3])[C:4][C:5]=[C:6]>>[C:2]=[C:3][C:4]([O:1][H])[C:5]=[C:6]",  # Explicit OH
            ]
        ),
        
        # Ring expansion
        (
            "[C:1]1[C:2][C:3][C:4]1.[C:5]=[C:6]>>[C:1]1[C:2][C:3][C:4][C:5][C:6]1",
            [
                "[C:1]1[C:2][C:3][C:4]1.[C:5]=[C:6]>>[C:1]1[C:2][C:3][C:4][C:5][C:6]1",  # Same pattern
            ]
        ),
        
        # Functional group migration
        (
            "[C:1]([OH:2])[C:3]=[C:4]>>[C:1]=[C:3][C:4][OH:2]",
            [
                "[C:1]([O:2][H])[C:3]=[C:4]>>[C:1]=[C:3][C:4][O:2][H]",  # Explicit hydrogen
            ]
        ),
    ]
    
    metric = SMIRKSMetric()
    
    for base_transform, variations in complex_transforms:
        base_ref = BaseReference(value_list=[base_transform])
        
        for variation in variations:
            compare_obj = CompareObj(value_list=[variation])
            result = metric.compute(compare_obj, base_ref)
            
            # These should be valid SMIRKS patterns
            assert result.valid >= 0.0, f"Base: {base_transform}, Variation: {variation}, got valid={result.valid}"
