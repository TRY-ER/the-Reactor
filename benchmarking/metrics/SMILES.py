from benchmarking.metrics.base import BaseMetric, CompareObj, BaseReference, MetricReturn
from rdkit import Chem
    
class SMILESMetric(BaseMetric):
    def __init__(self):
        super().__init__(name="SMILES Metric")

    def compute(self, compare_obj: CompareObj, base_ref: BaseReference) -> MetricReturn:
        num_found = len(compare_obj.value_list)
        num_valid = 0
        num_correct = 0
        # completely extent the list in lists to get number of maximum molecules in the base_ref.value_list
        extended_list = [smiles for sublist in base_ref.value_list for smiles in sublist] 
        max_mols = len(extended_list) 
        for values in compare_obj.value_list:
            if type(values) is not list:
                continue
            for smiles in values:
                if not isinstance(smiles, str):
                    continue 
                can_smiles = self.get_cannonicalized(smiles)
                if can_smiles == "<INVALID SMILES>":
                    continue
                num_valid += 1
                if can_smiles in extended_list:
                    num_correct += 1

        return MetricReturn(
            found=num_found/len(base_ref.value_list) if len(base_ref.value_list) > 0 else 0,
            valid=num_valid/max_mols if max_mols > 0 else 0,
            correct=num_correct/max_mols if max_mols > 0 else 0
        )

    def get_cannonicalized(self, smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol, canonical=True)
        return "<INVALID SMILES>"


if __name__ == "__main__":
    compare_obj = CompareObj(
        value_list=[
            ["CCO", "N[C@@H](C)C(=O)O", "invalid_smiles"],
            ["C1CCCCC1", "C1=CC=CC=C1"]
        ]
    )
    base_ref = BaseReference(
        value_list=[
            ["CCO", "CCN", "CCCN"],
            ["C1CCCCC1", "c1ccccc1"]
        ]
    )
    metric = SMILESMetric()
    result = metric.compute(compare_obj, base_ref)
    print(f"Num Found: {result.found}")
    print(f"Num Valid: {result.valid}")
    print(f"Num Correct: {result.correct}")