from benchmarking.metrics.base import BaseMetric, CompareObj, BaseReference, MetricReturn
from rdkit.Chem.rdChemReactions import ReactionFromSmarts
from rdkit import Chem
    
class SMIRKSMetric(BaseMetric):
    def __init__(self):
        super().__init__(name="SMIRKS Metric")

    def compute(self, compare_obj: CompareObj, base_ref: BaseReference) -> MetricReturn:
        num_found = len(compare_obj.value_list)
        num_valid = 0
        num_correct = 0
        for value, ref_value in zip(compare_obj.value_list, base_ref.value_list):
            if type(value) is not str:
                continue
            try:
                rxn = ReactionFromSmarts(value) 
                ref_rxn = ReactionFromSmarts(ref_value) 
                if rxn is None or ref_rxn is None:
                    continue
                num_valid += 1
                correct_reactants = self.compare_reactions(list(rxn.GetReactants()), list(ref_rxn.GetReactants()))
                correct_products = self.compare_reactions(list(rxn.GetProducts()), list(ref_rxn.GetProducts()))
                if correct_reactants and correct_products:
                    num_correct += 1

            except Exception as e:
                print("Error in parsing SMIRKS:", e)
                continue

        return MetricReturn(
            found=num_found/len(base_ref.value_list) if len(base_ref.value_list) > 0 else 0,
            valid=num_valid/len(base_ref.value_list) if len(base_ref.value_list) > 0 else 0,
            correct=num_correct/len(base_ref.value_list) if len(base_ref.value_list) > 0 else 0
        )

    def remove_square_brackets(self, smi):
        if smi.startswith('[') and smi.endswith(']'):
            return smi[1:-1]
        return smi

    def compare_reactions(self, reactant_mols, ref_reactant_mols):
        reactant_smiles = [Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True) for mol in reactant_mols]
        ref_reactant_smiles = [Chem.MolToSmiles(mol, isomericSmiles=True,  canonical=True) for mol in ref_reactant_mols]
        reactant_smiles = sorted(reactant_smiles)
        ref_reactant_smiles = sorted(ref_reactant_smiles)
        if len(reactant_smiles) != len(ref_reactant_smiles):
            return False
        for smi,ref_smi in zip(reactant_smiles, ref_reactant_smiles):
            smi = self.remove_square_brackets(smi)
            ref_smi = self.remove_square_brackets(ref_smi)
            if smi in ref_smi or ref_smi in smi or "*" in smi:
                continue
            else:
                return False
        return True

if __name__ == "__main__":
    compare_obj = CompareObj(
        value_list=[
            "[235#92:1].[1#0]>>[92#36:2].[141#56:3].[1#0].[1#0].[1#0]"
        ]
    )
    base_ref = BaseReference(
        value_list=[
            "[235#92].[1#0]>>[92#36].[141#56].[1#0].[1#0].[1#0]"
        ]
    )


    metric = SMIRKSMetric()
    result = metric.compute(compare_obj, base_ref)
    print(f"Num Found: {result.found}")
    print(f"Num Valid: {result.valid}")
    print(f"Num Correct: {result.correct}")