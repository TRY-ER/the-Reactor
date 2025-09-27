# import all metrics
from benchmarking.metrics.SMILES import SMILESMetric
from benchmarking.metrics.SMARTS import SMARTSMetric
from benchmarking.metrics.SMIRKS import SMIRKSMetric
from benchmarking.metrics.SMILES_rxn import SMILESRXNMetric
from benchmarking.metrics.base import BaseMetric, CompareObj, BaseReference, MetricReturn
from benchmarking.metrics.regress import mean_absolute_error, root_mean_squared_error, r_squared, mean_absolute_percentage_error

class Marker:
    def __init__(self, column_map: dict):
        self.column_map = column_map
        self.allowed_types = ["smiles", "smarts", "smirks", "smiles_rxn"]
        self.df = None

    def map_columns(self):
        pass

    def get_object(self, index: int, colname: str):
        pass

    def get_reference(self, index: int, colname: str):
        pass


    def call_metric(self, metric_type: str, compare_obj: CompareObj, base_ref: BaseReference) -> MetricReturn:
        if metric_type not in self.allowed_types:
            raise ValueError(f"Type {metric_type} not supported. Allowed types: {self.allowed_types}")

        if metric_type == "smiles":
            metric = SMILESMetric()
        elif metric_type == "smarts":
            metric = SMARTSMetric()
        elif metric_type == "smirks":
            metric = SMIRKSMetric()
        elif metric_type == "smiles_rxn":
            metric = SMILESRXNMetric()
        else:
            raise ValueError(f"Type {metric_type} not supported. Allowed types: {self.allowed_types}")

        result = metric.compute(compare_obj, base_ref)
        return result 

    def mark_one(self, index: int):
        if self.df is None:
            raise ValueError("DataFrame not set. Please set the DataFrame before marking.")
        if type not in self.allowed_types:
            raise ValueError(f"Type {type} not supported. Allowed types: {self.allowed_types}")

        types = ["smiles", "smarts", "smirks", "smiles_rxn"]
        results = {}
        
        for col_type in types:
            compare_obj = CompareObj(value_list=self.get_object(index, colname=col_type))
            base_ref = BaseReference(value_list=self.get_reference(index, colname=col_type)) 
            result = self.call_metric(type=type, compare_obj=compare_obj, base_ref=base_ref)
            results[col_type] = result

        return results
    

    def compute_result(self):
        pass

    def mark(self):
        if self.df is None:
            raise ValueError("DataFrame not set. Please set the DataFrame before marking.")
        all_results = []
        for idx in range(len(self.df)):
            result = self.mark_one(idx)
            all_results.append(result)
        
        final_result = self.compute_result(all_results)
        return final_result 