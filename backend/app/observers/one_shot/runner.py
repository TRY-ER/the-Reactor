import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from agents import Agent, Runner
from app.observers.one_shot.run_params import QUERY, INSTRUCTIONS
from app.observers.output_schema import OutputSchema
from agents.extensions.models.litellm_model import LitellmModel


ALLOWED_MODELS = ["groq", "xai"]

async def run_agent(
              model_type: str,
              model_name: str,
              model_key: str,
              agent_name: str = "one_shot_text_agent",
              query: str = QUERY,
              output_schema = None,
              instruction = None,
              *args,
              **kwargs):
    # model param initialization
    if model_type in ALLOWED_MODELS:
        model = LitellmModel(model=model_name, api_key=model_key)
    else:
        model = model_name


    # Agent initialization
    agent = Agent(
        name=agent_name,
        model=model,
        instructions=INSTRUCTIONS if instruction is None else instruction,
        tools=[],
        output_type=OutputSchema if output_schema is None else output_schema,
    )

    response = await Runner.run(
        agent,
        query
    )

    return response


# if __name__ == "__main__":
#     import dotenv
#     dotenv.load_dotenv()
#     import asyncio
#     api_key =  os.getenv("OPENAI_API_KEY")


#     from pydantic import BaseModel

#     class TestOutputSchema(BaseModel):
#         response: str

#     params = {
#         "observation_id": "test_observation_id",
#         "model_type": "openai",
#         "model_name": "gpt-5-2025-08-07",
#         "model_key": api_key,
#         "agent_name": "test_agent_name",
#         "query": "how are you ?",
#         "output_schema": TestOutputSchema,
#     }

#     asyncio.run(run_agent(**params))
