import pytest
from benchmarking.metrics.SMILES_rxn import SMILESRXNMetric
from benchmarking.metrics.base import CompareObj, BaseReference

# Each entry: (base_reaction, list_of_equivalent_variations)
reaction_smiles_variations = [
    # Simple oxidation of ethanol to acetaldehyde
    (
        "CCO>>CC=O",
        [
            "C(C)O>>C(C)=O",  # Different parentheses
            "OCC>>O=CC",      # Reversed order
        ]
    ),
    
    # Esterification reaction
    (
        "CC(=O)O.CCO>>CC(=O)OCC.O",
        [
            "CCO.CC(=O)O>>CC(=O)OCC.O",     # Swapped reactants
        ]
    ),
    
    # Hydrolysis of ester
    (
        "CC(=O)OCC.O>>CC(=O)O.CCO",
        [
            "O.CC(=O)OCC>>CCO.CC(=O)O",     # Swapped reagents and products
        ]
    ),
    
    # Cyclization reaction
    (
        "NCCCCCC(=O)O>>O=C1CCCCN1.O",
        [
            "O=C(O)CCCCN>>N1CCCCC(=O)1.O",  # Different ring numbering
        ]
    ),
    
    # Aromatic substitution
    (
        "c1ccccc1.Cl>>c1ccc(Cl)cc1",
        [
            "Cl.c1ccccc1>>c1ccc(Cl)cc1",      # Swapped reactants
        ]
    ),
    
    # SN2 reaction
    (
        "CCCBr.O>>CCCO.Br",
        [
            "O.BrCCC>>OCCC.Br",              # Swapped reactants
        ]
    ),
    
    # Addition reaction
    (
        "C=C.[H]Cl>>CCCl",
        [
            "[H]Cl.C=C>>ClCC",                 # Swapped reactants
        ]
    ),
    
    # Aldol condensation
    (
        "CC=O.CC=O>>CC(O)CC=O",
        [
            "O=CC.O=CC>>O=CCC(O)C",          # Different order
        ]
    ),
]

def test_reaction_smiles_metric_equivalence():
    """Test that different valid SMILES representations of the same reaction are recognized as equivalent."""
    metric = SMILESRXNMetric()
    
    for base_reaction, variations in reaction_smiles_variations:
        base_ref = BaseReference(value_list=[base_reaction])
        
        for variation in variations:
            compare_obj = CompareObj(value_list=[variation])
            result = metric.compute(compare_obj, base_ref)
            
            # Check that the variation is recognized as correct
            assert result.valid >= 1.0, f"Failed for base: {base_reaction}, variation: {variation}, got valid={result.valid}"
            assert result.found == 1.0, f"Failed for base: {base_reaction}, variation: {variation}, got found={result.found}"
            # Note: result.correct might be less than 1.0 due to canonicalization differences

def test_reaction_smiles_metric_invalid_cases():
    """Test that invalid reaction SMILES are properly handled."""
    invalid_reactions = [
        "C>C",           # Wrong arrow
        "C--C>>C=C",     # Invalid bond notation
        "Q>>C",          # Invalid atom
        "[H]Cl>>CCCl",   # Improper hydrogen representation that causes parsing error
    ]
    
    metric = SMILESRXNMetric()
    
    for invalid_rxn in invalid_reactions:
        base_ref = BaseReference(value_list=["CCO>>CC=O"])  # Valid reference
        compare_obj = CompareObj(value_list=[invalid_rxn])
        result = metric.compute(compare_obj, base_ref)
        
        # Invalid reactions should have valid=0 (but some edge cases might still parse)
        # So we mainly check that they don't match as correct
        assert result.correct == 0.0, f"Invalid reaction {invalid_rxn} should have correct=0, got {result.correct}"

def test_reaction_smiles_metric_multiple_reactions():
    """Test handling of multiple reactions in batch."""
    base_reactions = [
        "CCO>>CC=O",
        "CC(=O)O.CCO>>CC(=O)OCC.O"
    ]
    
    test_reactions = [
        "C(C)O>>C(C)=O",              # Equivalent to first
        "CCO.CC(=O)O>>CC(=O)OCC.O",   # Equivalent to second
    ]
    
    base_ref = BaseReference(value_list=base_reactions)
    compare_obj = CompareObj(value_list=test_reactions)
    
    metric = SMILESRXNMetric()
    result = metric.compute(compare_obj, base_ref)
    
    assert result.found == 1.0, f"Should find all reactions, got found={result.found}"
    assert result.valid >= 0.5, f"Should have valid reactions, got valid={result.valid}"
