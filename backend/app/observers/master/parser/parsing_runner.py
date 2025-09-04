import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../../..')))

import json
from app.harmony_adapter.openai_adapter import HarmonyOpenAIAdapter
from app.observers.master.parser.prompt import (
    QUERY,
    INSTRUCTIONS,
    VERIFIER_QUERY,
    VERIFIER_INSTRUCTION
)
from pydantic import BaseModel

ALLOWED_MODELS = [
    "groq/moonshotai/kimi-k2-instruct",
    "groq/openai/gpt-oss-120b",
    'groq/openai/gpt-oss-20b',
]


class ParsingOutputSchema(BaseModel):
    response: list[str]  # this is the list of the individual reaction texts

class VerifyingOutputSchem(BaseModel):
    response:  bool # True for valid reaction text, False for invalid reaction text


async def run_agent(
    model_type: str,
    model_name: str,
    model_key: str,
    query: str,
    agent_name: str = "parser_agent",
    instruction=INSTRUCTIONS,
    *args,
    **kwargs
):

    if model_type != "groq":
        raise ValueError(
            f"Model type {model_type} is not supported. Only 'groq' is supported.")
    if model_name not in ALLOWED_MODELS:
        raise ValueError(
            f"Model {model_name} is not allowed. Choose one between {ALLOWED_MODELS}")

    if model_name.startswith("groq/"):
        model_name = model_name.replace("groq/", "")

    agent = HarmonyOpenAIAdapter(
        model_name=model_name,
        api_key=model_key
    )

    agent.init_conversation(
        model_identity="You are a reaction parser agent",
        instructions=instruction,
        reasoning_effort="High"
    )

    query = QUERY.format(reaction_text=query)

    response = agent.invoke_parse(ParsingOutputSchema, query)
    # dummy_vals = ["A + B -> AB"] 
                #   "A precipitation reaction involves forming a solid in a solution or within another solid during a chemical reaction, occurring when ion concentrations surpass the solubility limit, resulting in an insoluble salt. This reaction can be facilitated by adding a precipitating agent or removing the solvent. Rapid processes yield amorphous or microcrystalline residues, while slower ones produce single crystals, which can also form via recrystallization from microcrystalline salts."]
    if response:
        values = response["content"][0]["text"]
        values = json.loads(values)

        new_values = []
        if values.get("response", None) is None:
            raise ValueError(
                "Invalid response from agent. 'response' key not found in the output.")
        else:
            for v in values.get("response"):
            # for v in dummy_vals: 
                print("model name in run_agent >>", model_name)
                verifier_agent = HarmonyOpenAIAdapter(
                    model_name=model_name,
                    api_key=model_key
                )
                verifier_agent.init_conversation(
                    model_identity="You are a reaction verifier agent",
                    instructions=VERIFIER_INSTRUCTION,
                    reasoning_effort="High"
                )

                if not isinstance(v, str):
                    raise ValueError(
                        "Invalid response from agent. 'response' key should contain a list of strings.")

                verifier_query = VERIFIER_QUERY.format(reaction_text=v)
                verify_response = verifier_agent.invoke_parse(VerifyingOutputSchem, verifier_query)
                if verify_response:
                    verify_values = verify_response["content"][0]["text"]
                    verify_values = json.loads(verify_values)
                    print("verify_query>>", verifier_query)
                    print("verify values >>", verify_values)
                    if verify_values.get("response", None) is None:
                        raise ValueError(
                            "Invalid response from verifier agent. 'response' key not found in the output.")
                    else:
                        if verify_values.get("response") == True:
                            new_values.append(v)
                        else:
                            print(
                                f"Invalid reaction text: {v}")
                            continue

        return new_values
    return [] 
                            



if __name__ == "__main__":
    import dotenv
    import asyncio
    from app.observers.reaction_iter import provide_reaction_details
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
        "instruction": INSTRUCTIONS
    }

    iter_path = "/home/kalki/src/valency/reactor/analysis/data/chemical_reaction_types_completed_urls_modified.csv"

    for i, reaction in enumerate(provide_reaction_details(
        iter_path,
        column_to_access=["Serial No.", "Reaction Type", "Generated Prompt"],
        custom_range=None
    )):
        generated_prompt = reaction["Generated Prompt"]
        params["query"] = generated_prompt
        print("index >>", i)
        # print("prompt >>", generated_prompt)
        result = asyncio.run(run_agent(**params))
        print("result >>", result)
        break
        # if i > 0 and i % 10 == 0:
        #     break
