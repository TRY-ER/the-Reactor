from pydantic import BaseModel

class MetricReturn(BaseModel):
    found: float 
    valid: float 
    correct : float 

class BaseReference(BaseModel):
    value_list: list 

class CompareObj(BaseModel):
    value_list: list

class BaseMetric:
    def __init__(self, name: str):
        self.name = name

    def compute(self, compare_obj: CompareObj, base_ref: BaseReference) -> MetricReturn:
        raise NotImplementedError("Subclasses should implement this!")