from pydantic import BaseModel

class OutputSchema(BaseModel):
    reactions_type: str
    reactions_text : list[str]
    reactions_composition_SMILES: list[list[str]] 
    reactions_SMILES: list[str]
    reactions_SMARTS: list[str]
    reactions_SMIRKS: list[str]