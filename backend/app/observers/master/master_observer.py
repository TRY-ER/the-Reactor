import os
import sys
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../../..')))

from app.observers.master.verbose import verbose_response
from app.observers.master.openai_runner import run_agent as master_run_agent
from app.observers.reaction_iter import provide_reaction_details
from app.observers.logger import PersistentDataLogger
from app.observers.model_params import MODEL_PARAMS, MODEL_PROGRESS_PARAMS
import uuid
import pandas as pd
import asyncio
import time

global INIT_TIME 
INIT_TIME = time.time()

global CURRENT_ITERATION
CURRENT_ITERATION = 0

async def agent_run(observation_params: dict) -> dict:
    response = await master_run_agent(**observation_params)
    assert type(response) == list, "the openai agent response should be of type list [list of dictionaries]"
    return response

async def run_observation(experiment_name: str,
                          data_range: list[int] | None = None) -> None:
    init_agent_name = "master-oss-"
    max_data = 67
    if data_range is not None:
        max_data = len(data_range)
    # get current file path + outputs + experiment_name
    logger_file_path = os.path.join(
        os.path.dirname(__file__), 'outputs', experiment_name)
    headers = [
        "sl_no",
        "model_type",
        "model_name",
        "reaction_type",
        "reactions_text",
        "reactions_composition_SMILES",
        "invalid_SMILE",
        "reactions_SMILES",
        "invalid_reaction_SMILES",
        "reactions_SMARTS",
        "invalid_reaction_SMARTS",
        "reactions_SMIRKS",
        "invalid_reaction_SMIRKS"
    ]
    logger = PersistentDataLogger(
        log_dir=logger_file_path,
        filename="findings.csv",
        headers=headers
    )
    dummy_df = pd.DataFrame(columns=headers)
    model_type_count = 0
    for model_type, model_values in MODEL_PARAMS.items():
        model_type_count += 1
        API_KEY = model_values.get("api_key", None)
        MODEL_NAMES = model_values.get("models", None)
        assert API_KEY is not None
        assert MODEL_NAMES is not None
        assert type(MODEL_NAMES) is list
        init_time = time.time()
        for i, model_name in enumerate(MODEL_NAMES):
            reaction_progress = 1 
            for values in provide_reaction_details(
                "/home/kalki/src/valency/reactor/analysis/data/chemical_reaction_types_completed_urls_modified.csv",
                    column_to_access=["Serial No.", "Reaction Type", "Generated Prompt"],
                    custom_range=data_range):
                reaction_type = values.get("Reaction Type", "")
                generated_prompt = values.get("Generated Prompt", "")
                sl_no = values.get("Serial No.", "")

                assert type(model_name) is str
                observation_params = {
                    "observation_id": model_name + "|" + reaction_type + "|" + str(uuid.uuid4()),
                    "model_type": model_type,
                    "model_name": model_name,
                    "model_key": API_KEY,
                    "agent_name": init_agent_name + model_name,
                    "query": generated_prompt,
                }

                try:
                    start_time = time.time()
                    datas = await agent_run(observation_params)
                    end_time = time.time()
                    execution_time = end_time - start_time
                    reaction_type = reaction_type


                    current_data = []
                    for data in datas:  
                        assert type(data) is dict, f"the data returned from the agent run is not Dict it's {type(data)}"
                        reactions_text = data.get("text")
                        smiles_part= data.get("SMILES")
                        smiles_rxn_part = data.get("SMILES_reaction")
                        smarts_rxn_part = data.get("SMARTS_reaction")
                        smirks_rxn_part = data.get("SMIRKS_reaction")

                        current_data.append({ 
                            "sl_no": sl_no,
                            "model_type": model_type,
                            "model_name": model_name,
                            "reaction_type": reaction_type,
                            "reactions_text": reactions_text,
                            "reactions_composition_SMILES": smiles_part["values"],
                            "invalid_SMILE": smiles_part["invalid_count"],
                            "reactions_SMILES": smiles_rxn_part["values"],
                            "invalid_reaction_SMILES": smiles_rxn_part["invalid_count"],
                            "reactions_SMARTS": smarts_rxn_part["values"],
                            "invalid_reaction_SMARTS": smarts_rxn_part["invalid_count"],
                            "reactions_SMIRKS": smirks_rxn_part["values"],
                            "invalid_reaction_SMIRKS": smirks_rxn_part["invalid_count"]
                        })

                    if len(current_data) == 0:
                        current_data = [{
                            "sl_no": sl_no,
                            "model_type": model_type,
                            "model_name": model_name,
                            "reaction_type": reaction_type,
                            "reactions_text": generated_prompt,
                            "reactions_composition_SMILES": [],
                            "invalid_SMILE": -1, 
                            "reactions_SMILES": [],
                            "invalid_reaction_SMILES": -1,
                            "reactions_SMARTS": [],
                            "invalid_reaction_SMARTS": -1,
                            "reactions_SMIRKS": [],
                            "invalid_reaction_SMIRKS": -1 
                        }]

                except Exception as e:
                    print("Error occured for model type >>", model_name)
                    print("Error occurred while processing response:", e)

                    current_data = [{
                        "sl_no": sl_no,
                        "model_type": model_type,
                        "model_name": model_name,
                        "reaction_type": reaction_type,
                        "reactions_text": [f"error >>{str(e)}"],
                        "reactions_composition_SMILES": ["error"],
                        "reactions_SMILES": ['error'],
                        "reactions_SMARTS": ['error'],
                        "reactions_SMIRKS": ['error'],
                    }]
                    execution_time =  "Error Execute" 
                
                for data in current_data:
                    dummy_df.loc[len(dummy_df)] = data 

                progress_map = {
                    "model_type_progress": model_type_count,
                    "model_name_progress": i + 1,
                    "reaction_type_progress": reaction_progress 
                }

                params_map = {
                    "max_model_types":  MODEL_PROGRESS_PARAMS['max_model_type'],
                    "max_num_models": len(MODEL_NAMES),
                    "max_reaction_types": max_data
                }

                try:
                    verbose_response(current_data, params_map,
                                 progress_map, execution_time)
                except Exception as e:
                    print("Error in verbose response:", e)
                reaction_progress += 1
    logger.log_df(dummy_df)


if __name__ == "__main__":
    import dotenv
    dotenv.load_dotenv()
    import asyncio
    asyncio.run(run_observation(
        experiment_name='experiment-29-08-2025-<02>',
        # data_range=[4]
    ))