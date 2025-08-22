"""
Comprehensive test cases for SMILES, SMARTS, and SMIRKS validators
Based on Daylight documentation examples and RDKit behavior
"""

import pytest
from app.validators import (
    validate_smiles, 
    validate_smarts, 
    validate_smirks,
    validate_smiles_detailed
)


class TestSMILESValidation:
    """Test cases for SMILES validation based on Daylight documentation."""
    
    def test_basic_smiles_from_documentation(self):
        """Test basic SMILES examples from Daylight documentation."""
        # Direct examples from https://www.daylight.com/dayhtml/doc/theory/theory.smiles.html
        valid_smiles = [
            # Basic molecules - from docs
            "CC",                   # ethane
            "C=C",                  # ethene  
            "C#C",                  # ethyne
            "CCO",                  # ethanol
            "CC(=O)O",              # acetic acid
            "c1ccccc1",             # benzene
            "C1CCCCC1",             # cyclohexane
            
            # Branching - from docs
            "CC(C)C",               # isobutane
            "CC(C)CC",              # isopentane  
            "CC(C)(C)C",            # neopentane
            
            # Ring examples - from docs
            "C1CCC1",               # cyclobutane
            "C1CCCC1",              # cyclopentane
            "C1CCCCC1",             # cyclohexane
            "C1CCCCCCC1",           # cycloheptane
            
            # Aromatic examples - from docs
            "c1ccccc1",             # benzene
            "c1ccncc1",             # pyridine
            "c1c[nH]cc1",           # pyrrole
            "c1ccc2ccccc2c1",       # naphthalene
            
            # Heteroatoms - from docs
            "CCO",                  # ethanol
            "CCN",                  # ethylamine
            "CCS",                  # ethanethiol
            "CCCF",                 # 1-fluoropropane
            "CCCl",                 # 1-chloropropane
            "CCBr",                 # 1-bromoethane
            "CCI",                  # 1-iodoethane
            
            # Charged species - from docs
            "[Na+]",                # sodium cation
            "[Cl-]",                # chloride anion
            "[NH4+]",               # ammonium
            "[OH-]",                # hydroxide
            "C[N+](C)(C)C",         # tetramethylammonium
            
            # Stereochemistry - from docs
            "C[C@H](N)C(=O)O",      # L-alanine
            "C[C@@H](N)C(=O)O",     # D-alanine  
            "C/C=C/C",              # trans-2-butene
            "C/C=C\\C",             # cis-2-butene
            
            # Isotopes - from docs
            "[13C]",                # carbon-13
            "[2H]",                 # deuterium
            "[3H]",                 # tritium
            
            # Disconnected structures - from docs
            "CCO.CCN",              # ethanol and ethylamine
            "[Na+].[Cl-]",          # sodium chloride
            "c1ccccc1.c1ccccc1",    # two benzene rings
        ]
        
        for smiles in valid_smiles:
            assert validate_smiles(smiles), f"Should be valid SMILES: {smiles}"
    
    def test_invalid_smiles_based_on_rdkit(self):
        """Test invalid SMILES that RDKit actually rejects."""
        invalid_smiles = [
            "",                     # empty string
            "   ",                  # whitespace only
            "C(",                   # unmatched opening parenthesis
            "C)",                   # unmatched closing parenthesis  
            "C[",                   # unmatched opening bracket
            "C]",                   # unmatched closing bracket
            "[",                    # lone opening bracket
            "]",                    # lone closing bracket
            "(",                    # lone opening parenthesis
            ")",                    # lone closing parenthesis
            "C1CC",                 # unclosed ring
            "CX",                   # invalid atom (X without brackets)
            "C=",                   # incomplete double bond
            "C#",                   # incomplete triple bond
            "C==C",                 # invalid double equals
            "C##C",                 # invalid double hash
            "C>>",                  # SMIRKS arrow in SMILES
            ">>C",                  # SMIRKS arrow in SMILES
            "C>C",                  # partial SMIRKS arrow in SMILES
        ]
        
        for smiles in invalid_smiles:
            assert not validate_smiles(smiles), f"Should be invalid SMILES: {smiles}"


class TestSMARTSValidation:
    """Test cases for SMARTS validation based on Daylight documentation."""
    
    def test_basic_smarts_from_documentation(self):
        """Test basic SMARTS examples from Daylight documentation."""
        # Examples from https://www.daylight.com/dayhtml/doc/theory/theory.smarts.html
        valid_smarts = [
            # Phenol example from docs
            "[OH]c1ccccc1",         # phenol-containing structures (key example)
            
            # All valid SMILES are valid SMARTS (from docs)
            "C",                    # carbon
            "CC",                   # ethane
            "c1ccccc1",             # benzene
            
            # Atomic number queries - from docs
            "[#6]",                 # carbon by atomic number
            "[#7]",                 # nitrogen by atomic number
            "[#8]",                 # oxygen by atomic number
            
            # Logical operators - from docs
            "[C,N]",                # carbon OR nitrogen
            "[!C]",                 # NOT carbon
            "[C&H2]",               # carbon AND exactly 2 hydrogens
            "[C;H2]",               # carbon AND exactly 2 hydrogens (high precedence)
            
            # Connection queries - from docs
            "[CX4]",                # carbon with 4 connections
            "[CX3]",                # carbon with 3 connections
            "[CX2]",                # carbon with 2 connections
            "[NX3]",                # nitrogen with 3 connections
            "[OX2]",                # oxygen with 2 connections
            
            # Hydrogen count - from docs
            "[CH3]",                # carbon with 3 hydrogens
            "[CH2]",                # carbon with 2 hydrogens
            "[CH1]",                # carbon with 1 hydrogen
            "[CH0]",                # carbon with 0 hydrogens
            "[NH2]",                # nitrogen with 2 hydrogens
            
            # Ring queries - from docs
            "[R]",                  # any ring atom
            "[!R]",                 # not in ring
            "[R1]",                 # in exactly 1 ring
            "[R2]",                 # in exactly 2 rings
            "[r5]",                 # in 5-membered ring
            "[r6]",                 # in 6-membered ring
            
            # Bond queries - from docs
            "C-C",                  # single bond
            "C=C",                  # double bond
            "C#C",                  # triple bond
            "C:C",                  # aromatic bond
            "C~C",                  # any bond
            "C@C",                  # ring bond
            "C!@C",                 # non-ring bond
            
            # Recursive SMARTS - from docs
            "[$(C=O)]",             # carbon in carbonyl
            "[$(CC)]",              # carbon connected to carbon
            "[$(c1ccccc1)]",        # atom in benzene ring
            "[!$(CC)]",             # NOT connected to carbon
            
            # Charge queries - from docs
            "[+]",                  # positive charge
            "[-]",                  # negative charge
            "[+1]",                 # charge +1
            "[-1]",                 # charge -1
        ]
        
        for smarts in valid_smarts:
            assert validate_smarts(smarts), f"Should be valid SMARTS: {smarts}"
    
    def test_invalid_smarts_based_on_rdkit(self):
        """Test invalid SMARTS that RDKit actually rejects."""
        invalid_smarts = [
            "",                     # empty string
            "   ",                  # whitespace only
            "(",                    # lone opening parenthesis
            ")",                    # lone closing parenthesis
            "[",                    # lone opening bracket
            "]",                    # lone closing bracket
            "[]",                   # empty brackets
            "[)",                   # mismatched bracket/parenthesis
            "C[",                   # unmatched opening bracket
            "C]",                   # unmatched closing bracket
            "C(",                   # unmatched opening parenthesis
            "C)",                   # unmatched closing parenthesis
        ]
        
        for smarts in invalid_smarts:
            assert not validate_smarts(smarts), f"Should be invalid SMARTS: {smarts}"


class TestSMIRKSValidation:
    """Test cases for SMIRKS validation based on Daylight documentation."""
    
    def test_basic_smirks_from_documentation(self):
        """Test basic SMIRKS examples from Daylight documentation."""
        # Examples from https://www.daylight.com/dayhtml/doc/theory/theory.smirks.html
        valid_smirks = [
            # Key SN2 example from docs
            "[C:1][Cl:2].[OH:3]>>[C:1][OH:3].[Cl:2]",
            
            # Basic transformations
            "[C:1]>>[C:1]",                         # identity transformation
            "[C:1][H:2]>>[C:1][OH:2]",              # H to OH substitution
            "[C:1]=[C:2]>>[C:1]-[C:2]",             # double to single bond
            "[C:1]#[C:2]>>[C:1]=[C:2]",             # triple to double bond
            
            # Functional group transformations
            "[C:1](=O)[OH:2]>>[C:1](=O)[O-:2]",     # acid deprotonation
            "[C:1][OH:2]>>[C:1]=O",                 # alcohol oxidation
            "[C:1]=O>>[C:1][OH:2]",                 # carbonyl reduction
            
            # Addition reactions
            "[C:1]=[C:2].[H:3][Br:4]>>[C:1]([H:3])[C:2][Br:4]",
            
            # Elimination reactions
            "[C:1][C:2][Br:3]>>[C:1]=[C:2].[H][Br:3]",
            
            # Multi-component reactions
            "[C:1]=[C:2].[C:3]=[C:4]>>[C:1][C:2][C:3][C:4]",
            "[C:1](=O)[OH:2].[N:3]>>[C:1](=O)[N:3].[OH:2]",
            
            # Stereochemical changes
            "[C@:1]>>[C@@:1]",                      # inversion
            "[C:1]/[C:2]>>[C:1]\\[C:2]",            # E/Z isomerization
            
            # RDKit accepts these patterns (empty reactants/products)
            "C>>",                                  # reactant only
            ">>C",                                  # product only
            
            # Valid but different atom mappings (RDKit allows this)
            "[C:1]>>[C:2]",                         # different mapping numbers
        ]
        
        for smirks in valid_smirks:
            assert validate_smirks(smirks), f"Should be valid SMIRKS: {smirks}"
    
    def test_invalid_smirks_based_on_rdkit(self):
        """Test invalid SMIRKS patterns that RDKit actually rejects."""
        invalid_smirks = [
            "",                             # empty string
            "   ",                          # whitespace only
            "C",                            # no reaction arrow
            "C>C",                          # single arrow (not >>)
            "C>>>C",                        # triple arrow
            "C>",                           # incomplete arrow
            ">C",                           # incomplete arrow
            "[C:1>>C:1]",                   # arrow inside brackets
            "C[",                           # unmatched bracket
            "C(",                           # unmatched parenthesis
            "[",                            # lone bracket
            "(",                            # lone parenthesis
        ]
        
        for smirks in invalid_smirks:
            assert not validate_smirks(smirks), f"Should be invalid SMIRKS: {smirks}"


class TestDetailedValidation:
    """Test detailed validation functions."""
    
    def test_detailed_smiles_validation(self):
        """Test the detailed SMILES validation function."""
        # Valid SMILES
        result = validate_smiles_detailed("CCO")
        assert result['valid'] is True
        assert result['error'] is None
        assert result['molecule'] is not None
        assert result['num_atoms'] == 3
        assert result['num_bonds'] == 2
        assert result['molecular_weight'] > 0
        
        # Invalid SMILES
        result = validate_smiles_detailed("C(")
        assert result['valid'] is False
        assert result['error'] is not None
        assert result['molecule'] is None
        assert result['num_atoms'] == 0
        
        # Empty string
        result = validate_smiles_detailed("")
        assert result['valid'] is False
        assert 'Empty' in result['error']
        
        # Non-string input
        result = validate_smiles_detailed(None)  # type: ignore
        assert result['valid'] is False
        assert 'string' in result['error']


class TestDocumentationCompliance:
    """Test compliance with specific documentation examples."""
    
    def test_smiles_documentation_examples(self):
        """Test exact SMILES examples from documentation."""
        # Key examples from Daylight SMILES documentation
        doc_examples = [
            ("CC", "ethane"),
            ("C=C", "ethene"),
            ("C#C", "ethyne"),
            ("CCO", "ethanol"),
            ("CC(=O)O", "acetic acid"),
            ("c1ccccc1", "benzene"),
            ("C1CCCCC1", "cyclohexane"),
            ("CC(C)C", "isobutane"),
            ("C[C@H](N)C(=O)O", "alanine with stereochemistry"),
            ("CCO.CCN", "disconnected components"),
            ("[Na+].[Cl-]", "ionic compound"),
        ]
        
        for smiles, name in doc_examples:
            assert validate_smiles(smiles), f"Documentation example {name} should be valid: {smiles}"
    
    def test_smarts_documentation_examples(self):
        """Test exact SMARTS examples from documentation."""
        # Key example from Daylight SMARTS documentation
        assert validate_smarts("[OH]c1ccccc1"), "Key phenol example from docs should be valid"
        
        # Other documented patterns
        doc_examples = [
            ("[#6]", "carbon by atomic number"),
            ("[C,N]", "carbon or nitrogen"),
            ("[!C]", "not carbon"),
            ("[CX4]", "sp3 carbon"),
            ("[$(C=O)]", "recursive carbonyl pattern"),
            ("[R]", "ring atom"),
            ("[r6]", "6-membered ring"),
        ]
        
        for smarts, name in doc_examples:
            assert validate_smarts(smarts), f"Documentation example {name} should be valid: {smarts}"
    
    def test_smirks_documentation_examples(self):
        """Test exact SMIRKS examples from documentation."""
        # Key SN2 example from Daylight SMIRKS documentation
        sn2_example = "[C:1][Cl:2].[OH:3]>>[C:1][OH:3].[Cl:2]"
        assert validate_smirks(sn2_example), "Key SN2 example from docs should be valid"
        
        # Other documented patterns
        doc_examples = [
            ("[C:1]=[C:2]>>[C:1]-[C:2]", "alkene reduction"),
            ("[C:1][OH:2]>>[C:1]=O", "alcohol oxidation"),
            ("[C:1]>>[C:1]", "identity transformation"),
        ]
        
        for smirks, name in doc_examples:
            assert validate_smirks(smirks), f"Documentation example {name} should be valid: {smirks}"


class TestEdgeCases:
    """Test edge cases and boundary conditions."""
    
    def test_type_safety(self):
        """Test type safety of validation functions."""
        invalid_types = [None, 123, [], {}, True, object()]
        
        for invalid_input in invalid_types:
            assert not validate_smiles(invalid_input), f"Should reject type: {type(invalid_input)}"
            assert not validate_smarts(invalid_input), f"Should reject type: {type(invalid_input)}"
            assert not validate_smirks(invalid_input), f"Should reject type: {type(invalid_input)}"
    
    def test_whitespace_handling(self):
        """Test whitespace handling."""
        # Leading/trailing whitespace should be handled
        assert validate_smiles("CCO"), "Basic SMILES should be valid"
        
        # Empty and whitespace-only strings should be invalid
        assert not validate_smiles(""), "Empty string should be invalid"
        assert not validate_smiles("   "), "Whitespace-only should be invalid"
        assert not validate_smiles("\n"), "Newline should be invalid"
        assert not validate_smiles("\t"), "Tab should be invalid"


if __name__ == "__main__":
    pytest.main([__file__])
