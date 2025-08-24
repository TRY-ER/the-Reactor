import os
import sys
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../../..')))

from app.observers.verbose import verbose_response
from app.observers.one_shot.runner import run_agent
from app.observers.reaction_iter import provide_reaction_details
from app.observers.one_shot.run_params import QUERY, INSTRUCTIONS
from app.observers.logger import PersistentDataLogger
from app.observers.model_params import MODEL_PARAMS, MODEL_PROGRESS_PARAMS
import uuid
import pandas as pd
import asyncio
import time


async def run_observation(experiment_name: str,
                          data_range: list[int] | None = None) -> None:
    init_agent_name = "one-shot-text"
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
        "reactions_SMILES",
        "reactions_SMARTS",
        "reactions_SMIRKS",
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
        for i, model_name in enumerate(MODEL_NAMES):
            reaction_progress = 1 
            for values in provide_reaction_details(
                "/home/kalki/src/valency/reactor/analysis/data/chemical_reaction_types_completed_urls_modified.csv",
                    column_to_access=["Serial No.", "Reaction Type", "Generated Prompt"],
                    custom_range=data_range):
                reaction_type = values.get("Reaction Type", "")
                generated_prompt = values.get("Generated Prompt", "")
                sl_no = values.get("Serial No.", "")

                query = QUERY.format(
                    reaction_type=reaction_type, generated_prompt=generated_prompt)
                assert type(model_name) is str
                observation_params = {
                    "observation_id": model_name + "|" + reaction_type + "|" + str(uuid.uuid4()),
                    "model_type": model_type,
                    "model_name": model_name,
                    "model_key": API_KEY,
                    "agent_name": init_agent_name + model_name,
                    "query": query,
                    "instruction": INSTRUCTIONS
                }

                try:
                    start_time = time.time()
                    response = await run_agent(**observation_params)
                    end_time = time.time()
                    execution_time = end_time - start_time

                    data = response.final_output.model_dump()
                    reaction_type = data.get("reactions_type", "unknown")
                    reactions_text = data.get("reactions_text", [])
                    reactions_composition_SMILES = data.get(
                        "reactions_composition_SMILES", [])
                    reactions_SMILES = data.get("reactions_SMILES", [])
                    reactions_SMARTS = data.get("reactions_SMARTS", [])
                    reactions_SMIRKS = data.get("reactions_SMIRKS", [])
                    current_data = {
                        "sl_no": sl_no,
                        "model_type": model_type,
                        "model_name": model_name,
                        "reaction_type": reaction_type,
                        "reactions_text": reactions_text,
                        "reactions_composition_SMILES": reactions_composition_SMILES,
                        "reactions_SMILES": reactions_SMILES,
                        "reactions_SMARTS": reactions_SMARTS,
                        "reactions_SMIRKS": reactions_SMIRKS
                    }

                except Exception as e:
                    # print("Error occured for model type >>", model_name)
                    # print("Error occurred while processing response:", e)

                    current_data = {
                        "sl_no": sl_no,
                        "model_type": model_type,
                        "model_name": model_name,
                        "reaction_type": reaction_type,
                        "reactions_text": [f"error >>{str(e)}"],
                        "reactions_composition_SMILES": ["error"],
                        "reactions_SMILES": ['error'],
                        "reactions_SMARTS": ['error'],
                        "reactions_SMIRKS": ['error'],
                    }
                    execution_time =  "Error Execute" 

                dummy_df.loc[len(dummy_df)] = current_data
                # print a formatted view of the findings
                # print(
                #     f"""
                #         Model Type: {model_type}
                #         Model Name: {model_name}
                #         Reaction Type: {reaction_type}
                #         Reactions Text: {reactions_text}
                #         Reactions Composition SMILES: {reactions_composition_SMILES}
                #         Reactions SMILES: {reactions_SMILES}
                #         Reactions SMARTS: {reactions_SMARTS}
                #         Reactions SMIRKS: {reactions_SMIRKS}
                #     """
                # )

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

                verbose_response(current_data, params_map,
                                 progress_map, execution_time)
                reaction_progress += 1
    logger.log_df(dummy_df)


if __name__ == "__main__":
    import dotenv
    dotenv.load_dotenv()
    import asyncio
    asyncio.run(run_observation(
        experiment_name='experiment-22-08-2025-<18-20>',
    ))
