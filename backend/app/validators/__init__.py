"""
This file contains validators that validate mostly chemical string annotations
like SMILES, SMARTS, and SMIRKS using RDKit for parsing and validation
"""
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors
from typing import Dict, Any


def validate_smiles(smiles: str) -> bool:
    """
    Validate SMILES string using RDKit parsing.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        bool: True if valid SMILES, False otherwise
    """
    if not isinstance(smiles, str) or not smiles.strip():
        return False
        
    try:
        mol = Chem.MolFromSmiles(smiles)
        return mol is not None
    except Exception:
        return False


def validate_smarts(smarts: str) -> bool:
    """
    Validate SMARTS string using RDKit parsing.
    
    Args:
        smarts: SMARTS string to validate
        
    Returns:
        bool: True if valid SMARTS, False otherwise
    """
    if not isinstance(smarts, str) or not smarts.strip():
        return False
        
    try:
        mol = Chem.MolFromSmarts(smarts)
        return mol is not None
    except Exception:
        return False


def validate_smirks(smirks: str) -> bool:
    """
    Validate SMIRKS string using RDKit parsing.
    
    Args:
        smirks: SMIRKS string to validate
        
    Returns:
        bool: True if valid SMIRKS, False otherwise
    """
    if not isinstance(smirks, str) or not smirks.strip():
        return False
        
    try:
        rxn = AllChem.ReactionFromSmarts(smirks)
        return rxn is not None
    except Exception:
        return False


def validate_smiles_detailed(smiles: str) -> Dict[str, Any]:
    """
    Detailed SMILES validation with error information.
    
    Args:
        smiles: SMILES string to validate
        
    Returns:
        Dict containing validation result and additional information
    """
    result = {
        'valid': False,
        'error': None,
        'molecule': None,
        'num_atoms': 0,
        'num_bonds': 0,
        'molecular_weight': 0.0
    }
    
    if not isinstance(smiles, str):
        result['error'] = 'Input must be a string'
        return result
        
    if not smiles.strip():
        result['error'] = 'Empty SMILES string'
        return result
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            result['error'] = 'Invalid SMILES: RDKit could not parse'
            return result
            
        result['valid'] = True
        result['molecule'] = mol
        result['num_atoms'] = mol.GetNumAtoms()
        result['num_bonds'] = mol.GetNumBonds()
        result['molecular_weight'] = Chem.Descriptors.MolWt(mol)
        
    except Exception as e:
        result['error'] = f'Parsing error: {str(e)}'
        
    return result