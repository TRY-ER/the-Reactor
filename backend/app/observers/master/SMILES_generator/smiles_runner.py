import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../../..')))

import json
from app.harmony_adapter.openai_adapter import HarmonyOpenAIAdapter
from app.observers.master.SMILES_generator.prompt import (
    QUERY,
    INSTRUCTIONS,
    RETRY_QUERY
)
from app.validators.main import StringValidator
from app.observers.master.base_schema_runner import run_agent as schema_run_agent
from pydantic import BaseModel

ALLOWED_MODELS = [
    "groq/moonshotai/kimi-k2-instruct",
    "groq/openai/gpt-oss-120b",
    'groq/openai/gpt-oss-20b',
]


class SMILESOutputSchema(BaseModel):
    SMILES_response: list[str]  # this is the list of the individual reaction texts


async def run_agent(
    model_type: str,
    model_name: str,
    model_key: str,
    query: str,
    agent_name: str = "SMILES_converter_agent",
    *args,
    **kwargs
):
    model_identity="YOU ARE AN EXPERT AGENT THAT GENERATES SMILES NOTATION FOR CHEMICAL REACTION TEXTS",

    params = {
        "model_type": model_type,
        "model_name": model_name,
        "model_key": model_key,
        "agent_name": agent_name,
        "query": query,
        "instruction": INSTRUCTIONS,
        "master_query": QUERY,
        "retry_query": RETRY_QUERY,
        "output_schema": SMILESOutputSchema,
        "model_identity": model_identity,
        "validator": StringValidator(),
        "validator_key": "SMILES_response"
    }

    result = await schema_run_agent(**params)

    return result
    
                            

if __name__ == "__main__":
    import dotenv
    import asyncio
    sys.path.insert(0, os.path.abspath(
        os.path.join(os.path.dirname(__file__), '../../..')))
    dotenv.load_dotenv()
    api_key = os.getenv("GROQ_API_KEY")

    params = {
        "model_type": "groq",
        "model_name": "groq/openai/gpt-oss-120b",
        "model_key": api_key,
        "agent_name": "openai_multistep_oss_agent",
        "query": "",
    }

    params["query"] = QUERY.format(reaction_text="Overall (net) reaction for the H2/Br2 chain system: H2 + Br2 -> 2 HBr.")

    result = asyncio.run(run_agent(**params))

    print("result >>", result)