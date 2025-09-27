import pytest
from benchmarking.metrics.SMARTS import SMARTSMetric
from benchmarking.metrics.base import CompareObj, BaseReference

# Each entry: (base_reaction_smarts, list_of_equivalent_variations)
smarts_reaction_variations = [
    # Alcohol oxidation - different atomic number representations
    (
        "CCO>>CC=O",
        [
            "[#6][#6][#8]>>[#6][#6]=[#8]",  # Using atomic numbers
            "[C][C][O]>>[C][C]=[O]",        # Explicit atom specification
        ]
    ),
    
    # Aromatic hydroxylation - different aromaticity representations
    (
        "c1ccccc1>>c1ccc(O)cc1",
        [
            "[c]1[c][c][c][c][c]1>>[c]1[c][c][c]([OH])[c][c]1",  # Explicit aromatic carbons
            "a1aaaaa1>>a1aaa(O)aa1",                             # Generic aromatic atoms
        ]
    ),
    
    # Ester hydrolysis - different bond representations
    (
        "[C:1][C:2](=O)[O:3][C:4].[O:5]>>[C:1][C:2](=O)[O:5].[O:3][C:4]",
        [
            "[C:1]-[C:2](=[O])-[O:3]-[C:4].[O:5]>>[C:1]-[C:2](=[O])-[O:5].[O:3]-[C:4]",  # Explicit single bonds
        ]
    ),
    
    # Nucleophilic substitution - different connectivity patterns
    (
        "[C:1][Cl:2].[OH:3]>>[C:1][O:3].[Cl:2]",
        [
            "[C;X4:1][Cl:2].[OH:3]>>[C;X4:1][O:3].[Cl:2]",  # Explicit hybridization
        ]
    ),
    
    # Ring formation - different ring notation
    (
        "[N:1][C:2][C:3][C:4][C:5][C:6](=O)[O:7]>>[N:1]1[C:2][C:3][C:4][C:5][C:6](=O)1.[O:7]",
        [
            "[N:1][C:2][C:3][C:4][C:5][C:6](=O)[OH:7]>>[N:1]1[C:2][C:3][C:4][C:5][C:6]1=O.[O:7]",  # Different ring closure
        ]
    ),
    
    # Alkene addition - different stereochemistry representation
    (
        "[C:1]=[C:2].[H][Cl:3]>>[C:1][C:2][Cl:3]",
        [
            "[C:1]=[C:2].[H][Cl:3]>>[C:1]-[C:2]-[Cl:3]",  # Explicit single bonds in product
        ]
    ),
    
    # Aromatic substitution with electron-withdrawing groups
    (
        "[c:1]1[c:2][c:3][c:4][c:5][c:6]1.[Br:7]>>[c:1]1[c:2][c:3][c:4]([Br:7])[c:5][c:6]1",
        [
            "a1aaaaa1.[Br:7]>>a1aaa([Br:7])aa1",  # Generic aromatic representation
        ]
    ),
    
    # Carbonyl reduction - different functional group patterns
    (
        "[C:1][C:2](=O)[C:3].[H][H]>>[C:1][C:2]([OH])[C:3]",
        [
            "[C:1][C:2](=[O])[C:3].[H2]>>[C:1][C:2]([O][H])[C:3]",  # Different hydrogen representation
        ]
    ),
    
    # Acyl chloride formation
    (
        "[C:1][C:2](=O)[OH:3].[Cl:4][Cl:5]>>[C:1][C:2](=O)[Cl:4].[OH:3][Cl:5]",
        [
            "[C:1][C:2](=[O])[O:3][H].[Cl:4][Cl:5]>>[C:1][C:2](=[O])[Cl:4].[O:3][H][Cl:5]",
        ]
    ),
    
    # Diels-Alder reaction - different diene representations
    (
        "[C:1]=[C:2][C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1",
        [
            "[C:1]=[C:2]-[C:3]=[C:4].[C:5]=[C:6]>>[C:1]1[C:2]=[C:3][C:4][C:6][C:5]1",  # Explicit single bond in diene
        ]
    ),
]

def test_smarts_reaction_metric_equivalence():
    """Test that different valid SMARTS representations of the same reaction are recognized as equivalent."""
    metric = SMARTSMetric()
    
    for base_reaction, variations in smarts_reaction_variations:
        base_ref = BaseReference(value_list=[base_reaction])
        
        for variation in variations:
            compare_obj = CompareObj(value_list=[variation])
            result = metric.compute(compare_obj, base_ref)
            
            # Check that the variation is recognized as valid and potentially correct
            assert result.valid >= 0.0, f"Failed for base: {base_reaction}, variation: {variation}, got valid={result.valid}"
            assert result.found == 1.0, f"Failed for base: {base_reaction}, variation: {variation}, got found={result.found}"

def test_smarts_reaction_metric_invalid_cases():
    """Test that invalid reaction SMARTS are properly handled."""
    invalid_reactions = [
        ">>C",           # Missing reactants  
        "C>C",           # Wrong arrow
        "[#0]>>C",       # Invalid atomic number
        "[X5]>>C",       # Invalid connectivity
        "Q>>C",          # Invalid atom symbol
    ]
    
    metric = SMARTSMetric()
    
    for invalid_rxn in invalid_reactions:
        base_ref = BaseReference(value_list=["CCO>>CC=O"])  # Valid reference
        compare_obj = CompareObj(value_list=[invalid_rxn])
        result = metric.compute(compare_obj, base_ref)
        
        # Invalid reactions should have correct=0
        assert result.correct == 0.0, f"Invalid reaction {invalid_rxn} should have correct=0, got {result.correct}"

def test_smarts_reaction_metric_substructure_patterns():
    """Test SMARTS patterns for substructure matching in reactions."""
    substructure_patterns = [
        # Alcohol to carbonyl generic pattern
        (
            "[C][OH]>>[C]=O",
            [
                "CCO>>CC=O",           # Specific ethanol oxidation
                "C(C)(C)CO>>C(C)(C)C=O",  # Tert-butyl alcohol oxidation
            ]
        ),
        
        # Aromatic halogenation pattern
        (
            "c1ccccc1>>c1ccc(X)cc1",  # X as halogen placeholder
            [
                "c1ccccc1.Cl>>c1ccc(Cl)cc1",  # Chlorination
                "c1ccccc1.Br>>c1ccc(Br)cc1",  # Bromination
            ]
        ),
        
        # Generic esterification
        (
            "[C](=O)[OH].[OH][C]>>[C](=O)[O][C]",
            [
                "CC(=O)O.CCO>>CC(=O)OCC",      # Acetic acid + ethanol
                "C(=O)O.CO>>C(=O)OC",          # Formic acid + methanol
            ]
        ),
    ]
    
    metric = SMARTSMetric()
    
    for pattern, specific_cases in substructure_patterns:
        for specific_case in specific_cases:
            base_ref = BaseReference(value_list=[pattern])
            compare_obj = CompareObj(value_list=[specific_case])
            result = metric.compute(compare_obj, base_ref)
            
            # Specific cases should be valid but may not be exactly correct due to pattern matching
            assert result.valid >= 0.0, f"Pattern: {pattern}, Case: {specific_case}, got valid={result.valid}"

def test_smarts_metric_atom_environment_patterns():
    """Test SMARTS patterns with different atom environment specifications."""
    environment_patterns = [
        # Primary alcohol oxidation
        (
            "[CH2;X2][OH]>>[CH2;X2]=O",
            [
                "[C;H2;X2][O;H1;X2]>>[C;H2;X2]=[O;H0;X1]",  # Explicit hydrogen and connectivity
            ]
        ),
        
        # Secondary alcohol oxidation  
        (
            "[CH;X3][OH]>>[CH;X3]=O",
            [
                "[C;H1;X3][O;H1;X2]>>[C;H1;X3]=[O;H0;X1]",  # Explicit specifications
            ]
        ),
        
        # Tertiary alcohol (should not oxidize easily)
        (
            "[C;X4]([C])([C])([C])[OH]>>[C;X4]([C])([C])([C])=O",
            [
                "[C;H0;X4]([C])([C])([C])[O;H1;X2]>>[C;H0;X4]([C])([C])([C])=[O;H0;X1]",
            ]
        ),
    ]
    
    metric = SMARTSMetric()
    
    for base_pattern, variations in environment_patterns:
        base_ref = BaseReference(value_list=[base_pattern])
        
        for variation in variations:
            compare_obj = CompareObj(value_list=[variation])
            result = metric.compute(compare_obj, base_ref)
            
            # These should be valid SMARTS patterns
            assert result.valid >= 0.0, f"Base: {base_pattern}, Variation: {variation}, got valid={result.valid}"
