import openai
from openai.types.responses import (
   ResponseReasoningItem,
   ResponseOutputMessage,
   ResponseOutputText,
   ResponseFunctionToolCall,
)
import os
import sys
import json 
from pydantic import BaseModel

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../..')))

from app.harmony_adapter.base import HarmonyAdapterBase
from app.harmony_adapter.groq_api_request import make_openai_groq_request

ALLOWED_MODELS = ["openai/gpt-oss-120b", 'openai/gpt-oss-20b']


class HarmonyOpenAIAdapter(HarmonyAdapterBase):
    def __init__(self, 
                 model_name: str,
                 api_key: str, 
                 function_tools: list = [],
                 function_map: dict = {}):
        super().__init__(model_name)
        self.api_key = api_key
        self.function_tools = function_tools
        self.function_map = function_map
        self.client = openai.OpenAI(
            base_url="https://api.groq.com/openai/v1",
            api_key=self.api_key
        )


    def set_tools(self, tools: list):
        self.function_tools = tools

    def set_tool_map(self, tool_map: dict):
        self.function_map = tool_map

    def convert_to_openai_messages(self):
        messages = []
        assert self.conversations is not None, "Conversations are not initialized."
        for c in self.conversations.messages:
            role = str(c.author.role.value)
            contents = []
            for i in c.content:
                content = i.to_dict()
                if role == "system":
                    reasoning_effort = content.get("reasoning_effort", "Low")
                    if type(reasoning_effort) is not str:
                        content["reasoning_effort"] = reasoning_effort.value
                    else:
                        content["reasoning_effort"] = reasoning_effort
                contents.append(content)
            messages.append({
                "role": role,
                "content": contents
            })
        return messages

    def parse_responses(self, response):
        for output in response.output:
            if isinstance(output, ResponseReasoningItem):
                if output.content and len(output.content) > 0:
                    self.add_message_from_dict({
                        "role": "assistant",
                        "content": [{
                            "type": "text",
                            "text": output.content[0].text
                        }]
                    })
            elif isinstance(output, ResponseOutputMessage):
                if output.content and len(output.content) > 0:
                    if isinstance(output.content[0], ResponseOutputText): 
                        self.add_message_from_dict({
                            "role": "assistant",
                            "content": [{
                                "type": "text",
                                "text": output.content[0].text
                            }]
                        })
            elif isinstance(output, ResponseFunctionToolCall):
                self.add_message_from_dict({
                    "role": "assistant",
                    "content": [{
                        "type": "text",
                        "text": output.arguments
                    }],
                    "channel": "commentary",
                    "recipient": f"functions.{output.name}",
                    "content_type": '<|constrain|> json'
                })
            else:
                print("Unknown output type:", type(output))
                print("output >>", output)

        return self.conversations.messages[-1].to_dict() if self.conversations else None

    def call_function(self, function_name: str, arguments: str):
        if function_name not in self.function_map:
            print(f"Function {function_name} is not available.")
            return None

        tool = self.function_map[function_name]
        args = json.loads(arguments)
        tool_response = tool(**args)
        return tool_response

    def parse_tool_params(self, response: dict):
        recipient = response.get("recipient", None)
        if recipient:
            return {
                "name": recipient.split(".")[-1],
                "arguments": response.get("content", [{}])[0].get("text", "")
            }
        return None

    def invoke(self, query: str, *args, **kwargs):
        response = self.client.responses.create(
            model=self.model_name,
            input=[
                *self.convert_to_openai_messages(),
                {
                    "role": "user",
                    "content": query 
                }
            ],
            tools=self.function_tools,
            tool_choice="auto"
        )
        return self.parse_responses(response)

    def invoke_tool_call(self, tool_res: dict):
        tool_name = tool_res["name"]
        tool_args = tool_res["arguments"]
        tool_response = self.call_function(tool_name, tool_args)
        self.add_message_from_dict({
            "role": "assistant",
            "name": f"functions.{tool_name}",
            "content": [
                {
                    "type": "text",
                    "text": str(tool_response)
                }
            ],
            "channel": "commentary",
        })

        response = self.client.responses.create(
            model = self.model_name,
            input = [
                *self.convert_to_openai_messages()
            ],
            tools=self.function_tools,
            tool_choice="auto"
        )

        return self.parse_responses(response)
    
    def invoke_parse(self, output_schema: BaseModel, query: str | None = None,  *args, **kwargs):
        if query:
            self.add_message_from_dict({
                "role": "user",
                "content": [
                    {
                        "type": "text",
                        "text": query
                    }
                ]
            })

        response = self.client.responses.parse(
            model=self.model_name,
            input=[
                *self.convert_to_openai_messages(),
            ],
            text_format=output_schema 
        )
        return self.parse_responses(response)

    def invoke_stream(self, *args, **kwargs):
        raise NotImplementedError(
            "Streaming is not supported in this adapter.")


# if __name__ == "__main__":
#     import dotenv
#     dotenv.load_dotenv()

#     import os
#     api_key = os.getenv('GROQ_API_KEY', "<key>")
#     adapter = HarmonyOpenAIAdapter(model_name="openai/gpt-oss-120b",
#                                    api_key=api_key)
#     adapter.init_conversation(
#         model_identity="You are a helpful assistant.",
#         instructions="Please assist the user with their queries.",
#         reasoning_effort="High",
#     )

#     adapter2 = HarmonyOpenAIAdapter(model_name="openai/gpt-oss-120b",
#                                    api_key=api_key)
#     adapter2.init_conversation(
#         model_identity="You are a helpful assistant.",
#         instructions="Please assist the user with their queries.",
#         reasoning_effort="High",
#     )
    

#     function_tools =[
#         {
#                 'type': 'function',
#                 "name": "get_current_weather",
#                 "description": "Get the current weather in a given location",
#                 "parameters": {
#                     "type": "object",
#                     "properties": {
#                         "location": {
#                             "type": "string",
#                             "description": "The city and state, e.g. San Francisco, CA"
#                         },
#                         "unit": {
#                             "type": "string",
#                             "enum": ["celsius", "fahrenheit"]
#                         }
#                     },
#                     "required": ["location"]
#                 }
#             }
#     ]
    
#     def get_current_weather(location: str, unit: str = "celsius"):
#         # Simulate a weather API call
#         weather_data = {
#             "location": location,
#             "temperature": 25,
#             "unit": unit
#         }
#         return weather_data

    # adapter.set_tools(function_tools)
    # adapter.set_tool_map({
    #     "get_current_weather": get_current_weather
    # })

    # response = adapter.invoke("What is the current wheather in New Delhi")

    # tool_res = adapter.parse_tool_params(response)

    # class ResponseSchema(BaseModel):
    #     location: str
    #     temperature: float
    #     unit: str

    # if tool_res:
    #     tool_response = adapter.invoke_tool_call(tool_res)
    #     print("tool response >>", tool_response)
    #     parsed_response = adapter.invoke_parse(ResponseSchema, None)
    #     print("parsed_response >>", parsed_response)
    # else:
    #     print("response >>", response)

    # query = "How are you ?"

    # parsed_response_2 = adapter2.invoke_parse(ResponseSchema, query)

    # conv = adapter2.decode_conversations()

