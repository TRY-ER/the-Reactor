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
from app.observers.master.service_formatter import format_for_stream



ALLOWED_MODELS = [
                  "groq/moonshotai/kimi-k2-instruct",
                  "groq/openai/gpt-oss-120b", 
                  'groq/openai/gpt-oss-20b',
                  ]


async def run_service_agent(
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

    if len(parser_response) == 0:
        yield format_for_stream("info", "No Reactions Found in the Text")
    
    else:
        for i,text in enumerate(parser_response):
            yield format_for_stream("data", f"{text}", aux={"info_type": "reaction_text", "reaction_index": i})

    sub_agent_params = {
        x: y for x, y in params.items() if x not in ["agent_name", "query"]
    }

    master_repo = []

    for i,text in enumerate(parser_response):
        yield format_for_stream("info", f"Starting for reaction {i}", aux={"info_type": "reaction_init", "reaction_index": i})

        assignable = {
            "text": text
        }

        assert type(text)  == str , "The response from the parser has to be in string form !"

        # execute extract SMILES

        yield format_for_stream("info", f"Starting SMILES Agent for reaction {i}", aux={"info_type": "smiles_init", "reaction_index": i})

        smiles_response = await smiles_run_agent(**sub_agent_params, query=text) 

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
            
            returnable = {
                "values": values,
                "invalid_count": invalid_count
            }

            assignable["SMILES"] = returnable 

            yield format_for_stream("data", f"Got SMILES values for reaction {i}", aux={"info_type": "smiles_update", "reaction_index": i, "returnable": returnable})

        yield format_for_stream("info", f"Starting Reaction SMILES Agent for reaction {i}", aux={"info_type": "rxn_smiles_init", "reaction_index": i})

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

            returnable = {
                "values": values,
                "invalid_count": invalid_count
            }

            assignable["SMILES_reaction"] = returnable 

            yield format_for_stream("data", f"Got Reaction SMILES values for reaction {i}", aux={"info_type": "smiles_rxn_update", "reaction_index": i, "returnable": returnable})

        yield format_for_stream("info", f"Starting Reaction SMARTS Agent for reaction {i}", aux={"info_type": "rxn_smarts_init", "reaction_index": i})

        smarts_rxn_response = await smarts_rxn_run_agent(**sub_agent_params, query=text)

        if smarts_rxn_response:
            content = smarts_rxn_response.get("content", [])
            content_text = content[0]["text"]
            content_obj = json.loads(content_text) 
            invalid_count = content[0].get("invalid", 0)

            assert "SMARTS_reaction_response" in content_obj, "SMARTS_reaction_response not found from the smarts agent"

            values = content_obj["SMARTS_reaction_response"]
            assert type(values) == list, "SMARTS_reaction_response is not a list"

            returnable =  {
                "values": values,
                "invalid_count": invalid_count
            }

            assignable["SMARTS_reaction"] =  returnable

            yield format_for_stream("data", f"Got Reaction SMARTS values for reaction {i}", aux={"info_type": "smarts_rxn_update", "reaction_index": i, "returnable": returnable})


        yield format_for_stream("info", f"Starting Reaction SMIRKS Agent for reaction {i}", aux={"info_type": "rxn_smirks_init", "reaction_index": i})

        smirks_rxn_response = await smirks_rxn_run_agent(**sub_agent_params, query=text)

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

            returnable = {
                "values": values,
                "invalid_count": invalid_count 
            }

            assignable["SMIRKS_reaction"] = returnable

            yield format_for_stream("data", f"Got Reaction SMIRKS values for reaction {i}", aux={"info_type": "smirks_rxn_update", "reaction_index": i, "returnable": returnable})

        master_repo.append(assignable)
        yield format_for_stream("info", f"Completed reaction {i}", aux={"info_type": "reaction_complete", "reaction_index": i, "returnable": returnable})

    yield format_for_stream("info", f"Completed the run", aux={"info_type": "run_complete"})
    yield {"final_csv": master_repo}


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

    query = r"""
    A synthesis (or combination) reaction involves two or more simple substances combining to form a more complex compound, typically represented as \( A + B \rightarrow AB \). This type of reaction is characterized by multiple reactants yielding a single product. An example is the formation of iron(II) sulfide from iron and sulfur:

    \[ 8 \text{Fe} + \text{S}_8 \rightarrow 8 \text{FeS} \]

    Another example includes the reaction of hydrogen and oxygen gases to produce water. Synthesis reactions are vital for forming complex molecules from simpler ones and are represented visually among four basic chemical reaction types: synthesis, decomposition, single replacement, and double replacement.
    """
    params["query"] = query

    async def run_and_collect():
        results = []
        async for result in run_service_agent(**params, output_schema=TestOutputSchema):
            results.append(result)
            print("r >>", result)
        return results

    results = asyncio.run(run_and_collect())