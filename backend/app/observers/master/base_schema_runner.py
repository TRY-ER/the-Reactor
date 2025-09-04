from pydantic import BaseModel
import os
import sys
import json

sys.path.insert(0, os.path.abspath(os.path.join(
    os.path.dirname(__file__), '../../..')))

from app.harmony_adapter.openai_adapter import HarmonyOpenAIAdapter
from app.validators.main import StringValidator

ALLOWED_MODELS = [
    "groq/moonshotai/kimi-k2-instruct",
    "groq/openai/gpt-oss-120b",
    'groq/openai/gpt-oss-20b',
]

VALIDATOR_KEY_REPO = [
    "SMILES_response",
    "SMILES_reaction_response",
    "SMARTS_reaction_response",
    "SMIRKS_reaction_response"
]

async def run_agent(
    model_type: str,
    model_name: str,
    model_key: str,
    query: str,
    master_query: str,
    retry_query: str,
    output_schema:  BaseModel ,
    validator: StringValidator,
    validator_key: str,
    agent_name: str = "parser_agent",
    instruction="",
    model_identity: str = "",
    max_retries: int = 3,
    *args,
    **kwargs
):

    assert validator_key in VALIDATOR_KEY_REPO, f"Validator key '{validator_key}' is not recognized. Must be one of {VALIDATOR_KEY_REPO}"

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
        model_identity=model_identity,
        instructions=instruction,
        reasoning_effort="High"
    )

    query = master_query.format(reaction_text=query)

    response = agent.invoke_parse(output_schema, query)

    if response:
        _responses = response["content"][0]["text"]
        _responses = json.loads(_responses)[validator_key]
        assert type(_responses) == list, f"Expected list but got {type(_responses)}"
        invalid_count = 0
        for smiles in _responses:
            if validator_key == "SMILES_response":
                is_valid, validation_message = validator.validate_smiles_grammar_candidate(smiles)
                if is_valid == False:
                    invalid_count += 1
                    agent.add_message_from_dict({
                        "role": "assistant",
                        "content": f"SMILES: {smiles} Validation State: {is_valid}, Message: {validation_message}"
                    })
                # print(f"Validation result for SMILES '{smiles}': {is_valid}, Message: {validation_message}")

            elif validator_key == "SMILES_reaction_response":
                is_valid, validation_message = validator.validate_smiles_reaction_grammar_candidate(smiles)
                if is_valid == False:
                    invalid_count += 1
                    agent.add_message_from_dict({
                        "role": "assistant",
                        "content": f"SMILES Reactions: {smiles} Validation State: {is_valid}, Message: {validation_message}"
                    })
                # print(f"Validation result for SMILES Reactions'{smiles}': {is_valid}, Message: {validation_message}")

            elif validator_key in ["SMARTS_reaction_response"]:
                is_valid, validation_message = validator.validate_smarts_reaction(smiles)
                if is_valid == False:
                    invalid_count += 1
                    agent.add_message_from_dict({
                        "role": "assistant",
                        "content": f"SMARTS: {smiles} Validation State: {is_valid}, Message: {validation_message}"
                    })
                # print(f"Validation result for SMARTS '{smiles}': {is_valid}, Message: {validation_message}")

            elif validator_key in ["SMIRKS_reaction_response"]:
                is_valid, validation_message = validator.validate_smirks_reaction(smiles)
                if is_valid == False:
                    invalid_count += 1
                    agent.add_message_from_dict({
                        "role": "assistant",
                        "content": f"SMIRKS: {smiles} Validation State: {is_valid}, Message: {validation_message}"
                    })
                # print(f"Validation result for SMIRKS '{smiles}': {is_valid}, Message: {validation_message}")
        
        if invalid_count > 0:
            while True:
                if max_retries == 0 :
                    if response is not None:
                        response["invalid"] = invalid_count
                        return response
                    else:
                        break
                retry_query = retry_query.format(reaction_text=query, invalid_count=invalid_count)
                response = agent.invoke_parse(output_schema, retry_query)
                if response:
                    _responses = response["content"][0]["text"]
                    _responses = json.loads(_responses)[validator_key]
                    assert type(_responses) == list, f"Expected list but got {type(_responses)}"
                    invalid_count = 0
                    for candidate in _responses:
                        if validator_key == "SMILES_response":
                            is_valid, validation_message = validator.validate_smiles_grammar_candidate(candidate)
                            if is_valid == False:
                                invalid_count += 1
                                agent.add_message_from_dict({
                                    "role": "assistant",
                                    "content": f"SMILES: {candidate} Validation State: {is_valid}, Message: {validation_message}"
                                })
                            # print(f"Validation result for SMILES '{candidate}': {is_valid}, Message: {validation_message}")

                        elif validator_key == "SMILES_reaction_response":
                            is_valid, validation_message = validator.validate_smiles_reaction_grammar_candidate(candidate)
                            if is_valid == False:
                                invalid_count += 1
                                agent.add_message_from_dict({
                                    "role": "assistant",
                                    "content": f"SMILES Reactions: {candidate} Validation State: {is_valid}, Message: {validation_message}"
                                })
                            # print(f"Validation result for SMILES Reactions'{candidate}': {is_valid}, Message: {validation_message}")

                        elif validator_key in ["SMARTS_reaction_response"]:
                            is_valid, validation_message = validator.validate_smarts_reaction(candidate)
                            if is_valid == False:
                                invalid_count += 1
                                agent.add_message_from_dict({
                                    "role": "assistant",
                                    "content": f"SMARTS: {candidate} Validation State: {is_valid}, Message: {validation_message}"
                                })
                            # print(f"Validation result for SMARTS '{candidate}': {is_valid}, Message: {validation_message}")

                        elif validator_key in ["SMIRKS_reaction_response"]:
                            is_valid, validation_message = validator.validate_smirks_reaction(candidate)
                            if is_valid == False:
                                invalid_count += 1
                                agent.add_message_from_dict({
                                    "role": "assistant",
                                    "content": f"SMIRKS: {candidate} Validation State: {is_valid}, Message: {validation_message}"
                                })
                            # print(f"Validation result for SMIRKS '{candidate}': {is_valid}, Message: {validation_message}")
                    if invalid_count == 0:
                        response["invalid"] = 0
                        return response

                max_retries -= 1
    if response and type(response) == dict:
        response["invalid"] = 0
    return response

# test script
if __name__ == "__main__":
    import dotenv
    import asyncio

    sys.path.insert(0, os.path.abspath(
        os.path.join(os.path.dirname(__file__), '../../..')))

    from app.observers.reaction_iter import provide_reaction_details
    from app.observers.master.SMIRKS_rxn_generator.prompt import (
        QUERY,
        INSTRUCTIONS,
        RETRY_QUERY
    )
    from app.validators.main import StringValidator
    dotenv.load_dotenv()
    api_key = os.getenv("GROQ_API_KEY")

    class SMILESOutputSchema(BaseModel):
        SMILES_response: list[str]  # this is the list of the individual reaction texts

    class RXNSMILESOutputSchema(BaseModel):
        SMILES_reaction_response: list[str]  # this is the list of the individual reaction SMILES texts

    class RXNSMARTSOutputSchema(BaseModel):
        SMARTS_reaction_response: list[str]  # this is the list of the individual reaction SMARTS texts
    
    class RXNSMIRKSOutputSchema(BaseModel):
        SMIRKS_reaction_response: list[str]  # this is the list of the individual reaction SMIRKS texts

    params = {
        "model_type": "groq",
        "model_name": "groq/openai/gpt-oss-120b",
        "model_key": api_key,
        "agent_name": "openai_multistep_oss_agent",
        "query": "Overall (net) reaction for the H2/Br2 chain system: H2 + Br2 -> 2 HBr.",
        "instruction": INSTRUCTIONS,
        "master_query": QUERY,
        "retry_query": RETRY_QUERY,
        "output_schema": RXNSMIRKSOutputSchema,
        "validator": StringValidator(),
        "validator_key": "SMIRKS_reaction_response"
    }

    result = asyncio.run(run_agent(**params))

    print("result >>", result)