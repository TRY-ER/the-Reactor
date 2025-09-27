import os
import sys
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../..')))

from app.harmony_adapter.groq_api_request import make_groq_request, make_openai_groq_request 
from app.harmony_adapter.base import HarmonyAdapterBase


class HarmonyAdapter(HarmonyAdapterBase):
    """
    The adapter follows harmony response format. 
    The major problem with the adapter is that the API endpoints does not allow for direct 
    output export for tools rather check for tool call formatting. 
    """
    def __init__(self, model_name: str, api_key: str):
        super().__init__(model_name)
        self.api_key = api_key

    def invoke(self, query: str, *args, **kwargs):
        if self.conversations is None:
            raise ValueError("Conversation is not initialized. Please call init_conversation first.")
        prompt = self.request_query(query)
        response = make_groq_request(model_name=self.model_name,
                                           api_key=self.api_key,
                                           query=prompt)
        return response

    def invoke_stream(self, *args, **kwargs):
        raise NotImplementedError("Streaming is not supported in this adapter.")


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
#         function_tools=[
#             {
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
#         ]
#     )

#     response = adapter.invoke(query="What's the weather like in San Francisco, CA?")
#     # response = adapter.invoke(query="how are you ?")
#     print("response >>", response)

