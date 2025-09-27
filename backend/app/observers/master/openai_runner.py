import os
import sys
import time
import json

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

from app.observers.one_shot.run_params import QUERY, INSTRUCTIONS
from app.harmony_adapter.openai_adapter import HarmonyOpenAIAdapter

# RUN Agent setup for all the subagents
from app.observers.master.parser.parsing_runner import run_agent as parser_run_agent 
from app.observers.master.SMILES_generator.smiles_runner import run_agent as smiles_run_agent
from app.observers.master.SMILES_rxn_generator.smiles_rxn_runner import run_agent as smiles_rxn_run_agent
from app.observers.master.SMARTS_rxn_generator.smarts_rxn_runner import run_agent as smarts_rxn_run_agent
from app.observers.master.SMIRKS_rxn_generator.smirks_rxn_runner import run_agent as smirks_rxn_run_agent



ALLOWED_MODELS = [
                  "groq/moonshotai/kimi-k2-instruct",
                  "groq/openai/gpt-oss-120b", 
                  'groq/openai/gpt-oss-20b',
                  ]



async def run_agent(
    model_type: str,
    model_name: str,
    model_key: str,
    agent_name: str = "master_openai_oss_agent",
    query: str = QUERY,
    output_schema = None,
    *args,
    **kwargs
    ):

    params = {
        "model_type": model_type,
        "model_name": model_name,
        "model_key": model_key,
        "agent_name": agent_name,
        "query": query,
    }

    parser_response = await parser_run_agent(**params)     

    assert type(parser_response) == list, "The parser is expected to return a list !"

    sub_agent_params = {
        x: y for x, y in params.items() if x not in ["agent_name", "query"]
    }

    master_repo = []

    for text in parser_response:
        assignable = {
            "text": text
        }
        assert type(text)  == str , "The response from the parser has to be in string form !"

        # execute extract SMILES

        smiles_response = await smiles_run_agent(**sub_agent_params, query=text) 

        # print("SMILES PART")
        # print("-----------")

        if smiles_response:
            content = smiles_response.get("content", [])
            content_text = content[0]["text"]
            content_obj = json.loads(content_text) 
            invalid_count = content[0].get("invalid", 0)

            assert "SMILES_response" in content_obj, "SMILES_response not found from the smiles agent"

            values = content_obj["SMILES_response"]
            assert type(values) == list, "SMILES_response is not a list"

            # print('values received >>', values)
            # print('invalid count >>', invalid_count)

            assignable["SMILES"] = {
                "values": values,
                "invalid_count": invalid_count
            }

        # execute extract SMILES reaction form

        print("[--] Sleeping for 10 seconds to ensure no rate limit is triggered")
        time.sleep(10)

        smiles_rxn_response = await smiles_rxn_run_agent(**sub_agent_params, query=text)

        # print("SMILES REACTION PART")
        # print("-----------")

        if smiles_rxn_response:
            content = smiles_rxn_response.get("content", [])
            content_text = content[0]["text"]
            content_obj = json.loads(content_text) 
            invalid_count = content[0].get("invalid", 0)

            assert "SMILES_reaction_response" in content_obj, "SMILES_reaction_response not found from the smiles agent"

            values = content_obj["SMILES_reaction_response"]
            assert type(values) == list, "SMILES_reaction_response is not a list"

            # print('values received >>', values)
            # print('invalid count >>', invalid_count)

            assignable["SMILES_reaction"] = {
                "values": values,
                "invalid_count": invalid_count
            }

        # execute extract SMARTS reaction form

        print("[--] Sleeping for 10 seconds to ensure no rate limit is triggered")
        time.sleep(10)


        smarts_rxn_response = await smarts_rxn_run_agent(**sub_agent_params, query=text)

        # print("SMARTS REACTION PART")
        # print("-----------")

        if smarts_rxn_response:
            content = smarts_rxn_response.get("content", [])
            content_text = content[0]["text"]
            content_obj = json.loads(content_text) 
            invalid_count = content[0].get("invalid", 0)

            assert "SMARTS_reaction_response" in content_obj, "SMARTS_reaction_response not found from the smarts agent"

            values = content_obj["SMARTS_reaction_response"]
            assert type(values) == list, "SMARTS_reaction_response is not a list"

            # print('values received >>', values)
            # print('invalid count >>', invalid_count)

            assignable["SMARTS_reaction"] = {
                "values": values,
                "invalid_count": invalid_count
            }

        # execute extract SMIRKS reaction form

        print("[--] Sleeping for 10 seconds to ensure no rate limit is triggered")
        time.sleep(10)

        smirks_rxn_response = await smirks_rxn_run_agent(**sub_agent_params, query=text)

        # print("SMIRKS REACTION PART")
        # print("-----------")


        if smirks_rxn_response:
            content = smirks_rxn_response.get("content", [])
            content_text = content[0]["text"]
            content_obj = json.loads(content_text) 
            invalid_count = content[0].get("invalid", 0)

            assert "SMIRKS_reaction_response" in content_obj, "SMIRKS_reaction_response not found from the smirks agent"

            values = content_obj["SMIRKS_reaction_response"]
            assert type(values) == list, "SMIRKS_reaction_response is not a list"

            # print('values received >>', values)
            # print('invalid count >>', invalid_count)

            assignable["SMIRKS_reaction"] = {
                "values": values,
                "invalid_count": invalid_count 
            }

        master_repo.append(assignable)
    return master_repo

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
        "model_name": "groq/openai/gpt-oss-120b",
        "model_key": api_key,
        "agent_name": "one_shot_openai_oss_agent",
        "query": "",
    }

    query = """
    A synthesis (or combination) reaction involves two or more simple substances combining to form a more complex compound, typically represented as \( A + B \rightarrow AB \). This type of reaction is characterized by multiple reactants yielding a single product. An example is the formation of iron(II) sulfide from iron and sulfur:

    \[ 8 \text{Fe} + \text{S}_8 \rightarrow 8 \text{FeS} \]

    Another example includes the reaction of hydrogen and oxygen gases to produce water. Synthesis reactions are vital for forming complex molecules from simpler ones and are represented visually among four basic chemical reaction types: synthesis, decomposition, single replacement, and double replacement.
    """
    params["query"] = query

    result = asyncio.run(run_agent(**params, output_schema=TestOutputSchema))
    print("result >>", result)