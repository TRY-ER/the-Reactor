import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))


from app.observers.one_shot.run_params import QUERY, INSTRUCTIONS
from app.harmony_adapter.openai_adapter import HarmonyOpenAIAdapter

ALLOWED_MODELS = [
                  "groq/moonshotai/kimi-k2-instruct",
                  "groq/deepseek-r1-distill-llama-70b",
                  "groq/openai/gpt-oss-120b", 
                  'groq/openai/gpt-oss-20b',
                  ]

async def run_agent(
    model_type: str,
    model_name: str,
    model_key: str,
    agent_name: str = "one_shot_openai_oss_agent",
    query: str = QUERY,
    output_schema = None,
    instruction = INSTRUCTIONS,
    *args,
    **kwargs
    ):

    if model_type != "groq":
        raise ValueError(f"Model type {model_type} is not supported. Only 'groq' is supported.")    
    if model_name not in ALLOWED_MODELS:
        raise ValueError(f"Model {model_name} is not allowed. Choose one between {ALLOWED_MODELS}")

    if model_name.startswith("groq/"):
        model_name = model_name.replace("groq/", "")

    agent = HarmonyOpenAIAdapter(
        model_name=model_name,
        api_key=model_key
    )

    agent.init_conversation(
        model_identity="",
        instructions=instruction,
        reasoning_effort="High"
    )


    response = agent.invoke(query)


    if output_schema is not None:
        return agent.invoke_parse(output_schema)

    return response


if __name__ == "__main__":
    import dotenv
    import asyncio
    sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))
    dotenv.load_dotenv()
    api_key =  os.getenv("GROQ_API_KEY")

    from pydantic import BaseModel

    class TestOutputSchema(BaseModel):
        response: str

    params = {
        "model_type": "groq",
        "model_name": "openai/gpt-oss-20b",
        "model_key": api_key,
        "agent_name": "one_shot_openai_oss_agent",
        "query": "What is the capital of France?",
        "instruction": "Provide a detailed response."
    }

    result = asyncio.run(run_agent(**params, output_schema=TestOutputSchema))

    print("result >>", result)