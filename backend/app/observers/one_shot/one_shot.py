import os
import sys
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../../..')))

from app.observers.one_shot.runner import run_agent
from app.observers.reaction_iter import provide_reaction_details
from app.observers.one_shot.run_params import QUERY, INSTRUCTIONS
from app.observers.logger import PersistentDataLogger
from app.observers.model_params import MODEL_PARAMS
import asyncio
import pandas as pd
import uuid

async def run_observation(experiment_name: str) -> None:
    init_agent_name = "one-shot-text"
    # get current file path + outputs + experiment_name
    logger_file_path = os.path.join(
        os.path.dirname(__file__), 'outputs', experiment_name)
    headers = [
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
    for model_type, model_values in MODEL_PARAMS.items():
        API_KEY = model_values.get("api_key", None)
        MODEL_NAMES = model_values.get("models", None)
        assert API_KEY is not None
        assert MODEL_NAMES is not None
        assert type(MODEL_NAMES) is list
        for model_name in MODEL_NAMES:
            for values in provide_reaction_details(
                "/home/kalki/src/valency/reactor/analysis/data/chemical_reaction_types_completed_urls_modified.csv",
                    column_to_access=["Serial No.", "Reaction Type", "Generated Prompt"]):
                reaction_type = values.get("Reaction Type", "")
                generated_prompt = values.get("Generated Prompt", "")

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
                    response = await run_agent(**observation_params)
                    data = response.final_output.model_dump()
                    reaction_type = data.get("reactions_type", "unknown")
                    reactions_text = data.get("reactions_text", [])
                    reactions_composition_SMILES = data.get(
                        "reactions_composition_SMILES", [])
                    reactions_SMILES = data.get("reactions_SMILES", [])
                    reactions_SMARTS = data.get("reactions_SMARTS", [])
                    reactions_SMIRKS = data.get("reactions_SMIRKS", [])
                    dummy_df.loc[len(dummy_df)] = {
                        "model_type": model_type,
                        "model_name": model_name,
                        "reaction_type": reaction_type,
                        "reactions_text": reactions_text,
                        "reactions_composition_SMILES": reactions_composition_SMILES,
                        "reactions_SMILES": reactions_SMILES,
                        "reactions_SMARTS": reactions_SMARTS,
                        "reactions_SMIRKS": reactions_SMIRKS
                    }
                    # print a formatted view of the findings
                    print(
                        f"""
                            Model Type: {model_type}
                            Model Name: {model_name}
                            Reaction Type: {reaction_type}
                            Reactions Text: {reactions_text}
                            Reactions Composition SMILES: {reactions_composition_SMILES}
                            Reactions SMILES: {reactions_SMILES}
                            Reactions SMARTS: {reactions_SMARTS}
                            Reactions SMIRKS: {reactions_SMIRKS}
                        """
                    )
                except Exception as e:
                    print("Error occured for model type >>", model_type)
                    print("Error occurred while processing response:", e)
                break
            break
        break
    logger.log_df(dummy_df)



if __name__ == "__main__":
    import dotenv
    dotenv.load_dotenv()
    import asyncio
    asyncio.run(run_observation(experiment_name='test-experiment'))
