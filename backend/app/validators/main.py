import os
import sys
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../..')))

from rdkit import Chem
from rdkit.Chem import rdChemReactions

from app.validators.logger_wrapper import capture_rdkit_errors
from app.validators.smirk_text_val import validate_smirk_text

class StringValidator:
    """
    Example -> [5*]NCN[5*].*Nc1ccc(NCCC[4*])cc1|[4*][*:1].[5*][*:2]>>[*:1][*:2] |[4*]-[*:1].[5*]-[*:2]>>[$([C&!D1&!$(C=*)]-&!@[#6]):1]-&!@[$([N&!D1&!$(N=*)&!$(N-[!#6&!#16&!#0&!#1])&!$([N&R]@[C&R]=O)]):2] 

    Part 1 -> SMILES ->  [5*]NCN[5*].*Nc1ccc(NCCC[4*])cc1
    Part 2 -> Reaction SMARTS ->  [4*]-[*:1].[5*]-[*:2]>>[$([C&!D1&!$(C=*)]-&!@[#6]):1]-&!@[$([N&!D1&!$(N=*)&!$(N-[!#6&!#16&!#0&!#1])&!$([N&R]@[C&R]=O)]):2]

    Validation Layer:

    1. Grammar Validation
        - Check for valid SMILES, SMARTS syntax
        - Ensure correct use of atom and bond symbols
    2. Syntax Validation
        - Ensure right indexing from the reaction components
        - Check for correct use of brackets and parentheses and transformations
        - Check for placement of "|" in the final representation
    """

    def __init__(self):
        pass 
    
    @capture_rdkit_errors
    def validate_smiles_grammar_candidate(self, str_part: str) -> tuple[bool, str]:
        try:
            mol = Chem.MolFromSmiles(str_part)
            
            if mol is None:
                return False, "Invalid SMILES received"
            
            return True, "Valid SMILES"
        except Exception as e:
            return False, str(e) 

    @capture_rdkit_errors
    def validate_smiles_reaction_grammar_candidate(self, str_part: str) -> tuple[bool, str]:
        try:
            assert ">>" in str_part, "Invalid reaction SMILES format [ must contain `>>` ]"
            reactant = str_part.split(">>")[0].strip()
            product = str_part.split(">>")[1].strip()

            re_mol = Chem.MolFromSmiles(reactant)

            if re_mol is None:
                return False, "Invalid SMILES received"
            
            pro_mol = Chem.MolFromSmiles(product)

            if pro_mol is None:
                return False, "Invalid SMILES received"

            return True, "Valid SMILES"
        except Exception as e:
            return False, str(e) 

    @capture_rdkit_errors
    def validate_smarts_reaction(self, str_part: str) -> tuple[bool, str]:
        try:
            rxn = rdChemReactions.ReactionFromSmarts(str_part) 
            if rxn is None:
                return False, "Invalid Reaction SMARTS"
            return True, "Valid Reaction SMARTS"
        except Exception as e:
            return False, str(e)


    def validate_smirks_reaction(self, str_part: str) -> tuple[bool, str]:
        try:
            validate_smirk_text(str_part)
            return True, "Valid SMIRKS"
        except Exception as e:
            return False, str(e)


if __name__ == "__main__":
    validator = StringValidator()
    
    # Test invalid SMILES
    print("Testing invalid SMILES:")
    smiles_result = validator.validate_smiles_grammar_candidate("[5*]NC(N[5*].[*]Nc1ccc(NCCC[4*])cc1")
    print("SMILES result >>", smiles_result)

    # Test valid reaction SMILES
    print('Testing invalid SMILES:')
    smiles_result = validator.validate_smiles_grammar_candidate("[H][H].BrBr>>[H]Br.[H]Br")
    print("SMILES reaction result >>", smiles_result)

    # Test invalid SMARTS
    print("\nTesting invalid SMARTS:")
    smarts_result = validator.validate_smarts_reaction("invalid>>smarts")
    print("SMARTS result >>", smarts_result)
    
    # Test valid SMILES
    print("\nTesting valid SMILES:")
    valid_smiles_result = validator.validate_smiles_grammar_candidate("CCO")
    print("Valid SMILES result >>", valid_smiles_result)

    # Test valid SMIRKS
    candidate = "[*:1][C@:2]([*:3])([*:4])[*:5]>>[*:1][C@@:2]([*:3])([*:4])[*:5]" # chiral changes
    valid_smirks_result = validator.validate_smirks_reaction(candidate)
    print("Valid SMIRKS result >>", valid_smirks_result)


    # Test invalid SMIRKS
    candidate = "[C:1](=[O:2])[Cl:3].[H:99][NH:4][C:0]>> [C:1](=[O:2])[NH:4][C:0].[Cl:3]"
    invalid_smirks_result = validator.validate_smirks_reaction(candidate)
    print("Invalid SMIRKS result >>", invalid_smirks_result)