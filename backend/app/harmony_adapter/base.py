from openai_harmony import (
    SystemContent,
    DeveloperContent,
    Conversation,
    ReasoningEffort,
    ToolDescription,
    Message,
    Role,
    load_harmony_encoding,
    HarmonyEncodingName,
    Author
)
import time


class HarmonyAdapterBase:
    def __init__(self, 
                 model_name: str):
        self.model_name = model_name
        self.conversations = None

    def format_system_message(self, 
                              model_identity: str | None,
                              reasoning_effort: str = "High", 
                              ):
        if reasoning_effort not in ["High", "Medium", "Low"]:
            raise ValueError("reasoning_effort must be 'High', 'Medium', or 'Low'")
        reasoning_effort = ReasoningEffort(reasoning_effort)
        content = SystemContent.new().with_reasoning_effort(reasoning_effort).with_conversation_start_date(f"{time.strftime('%Y-%m-%d')}") # the conversation data has to be in "YYYY-MM-DD" format
        if model_identity:
            content = content.with_model_identity(model_identity)
        return (
            content
        )

    def _format_functions(self, function_tools: list[dict]):
        """
        The function has to be formatted in following format.
        "name": "function_name",
        "description": "function_description",
        "parameters": {
            "type": "object",
            "properties": {
                "param1": {
                    "type": "string",
                    "description": "Description of param1"
                },
                "param2": {
                    "type": "integer",
                    "description": "Description of param2"
                }
            },
            "required": ["param1", "param2"]
        }
        """
        tools = []
        for func in function_tools:
            tool = ToolDescription.new(
                name=func["name"],
                description=func.get("description", ""),
                parameters=func.get("parameters", {})
            )
            tools.append(tool)
        return tools


    def format_developer_message(self, instructions: str, function_tools: list[dict] | None = None):
        content = DeveloperContent.new().with_instructions(instructions)
        if function_tools:
            content = content.with_function_tools(self._format_functions(function_tools))
        return content

    def format_conversation_message(self, messages: list[Message]):
        content = Conversation.from_messages(
           messages=messages
        )
        return content

    def init_conversation(self,
                          #system params
                          model_identity: str | None,
                          instructions: str,
                          reasoning_effort: str = "High",
                          function_tools: list[dict] | None = None
                          ):
        system_message = self.format_system_message(model_identity=model_identity, reasoning_effort=reasoning_effort)    
        developer_message = self.format_developer_message(instructions=instructions, function_tools=function_tools)
        message_list = [
            Message.from_role_and_content(Role.SYSTEM, system_message),
            Message.from_role_and_content(Role.DEVELOPER, developer_message),
        ]
        self.conversations = self.format_conversation_message(message_list)
        return self.conversations

    def decode_conversations(self):
        if self.conversations is None:
            raise ValueError("Conversation is not initialized. Please call init_conversation first.")
        encoding = load_harmony_encoding(HarmonyEncodingName.HARMONY_GPT_OSS)
        prompt_token = encoding.render_conversation(
            conversation=self.conversations)
        prompt = encoding.decode(prompt_token)       
        return prompt

    def add_message_from_dict(self, message_dict: dict):
        """
        Adds a message to the conversation from a dictionary representation.
        The message_dict will be in the following format.
        {
            "role": "user",
            "name": "<name>", # optional
            "content": [
                {
                    "type": "text", # type could be of text, system_content, developer_content
                    "text": "Hello, how can I help you?"
                }
            "chennel": "default", # optional
            "recipient": "all", # optional
            "content_type": "text", # optional 
            ]
        }
        """
        message = Message.from_dict(message_dict)
        self.add_message(message)

    def add_message(self, message: Message):
        if self.conversations is None:
            raise ValueError("Conversation is not initialized. Please call init_conversation first.")
        self.conversations.messages.append(message)
        return self.conversations

    def request_query(self, query: str):
        if self.conversations is None:
            raise ValueError("Conversation is not initialized. Please call init_conversation first.")
        user_message = Message.from_role_and_content(Role.USER, query)
        self.conversations.messages.append(user_message)
        encoding = load_harmony_encoding(HarmonyEncodingName.HARMONY_GPT_OSS)
        prompt_token = encoding.render_conversation_for_completion(self.conversations, Role.ASSISTANT)
        prompt = encoding.decode(prompt_token)       
        return prompt
    

    def invoke(self, *args, **kwargs):
        raise NotImplementedError("Subclasses must implement the invoke method.")

    def invoke_stream(self, *arts, **kwargs):
        raise NotImplementedError("Subclasses must implement the invoke_stream method.")


# testing script
# if __name__ == "__main__":
#     adapter = HarmonyAdapterBase(model_name="some_model_name")
#     conversation = adapter.init_conversation(
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

#     # user query addition
#     adapter.add_message_from_dict({
#         "role": "user",
#         "content": [
#             {
#                 "type": "text",
#                 "text": "What's the weather like in Boston?"
#             }
#         ]
#     })

#     # tool call addition
#     adapter.add_message_from_dict({
#         "role": "assistant",
#         "content": [
#             {
#                 "type": "text",
#                 "text": '{"location": "Boston", "unit": "celsius"}'
#             }
#         ],
#         "channel": "commentary",
#         "recipient": "functions.get_current_weather",
#         "content_type": '<constrain> json'
#     })

#     test_ext_message = Message.from_role_and_content(Role.ASSISTANT, '{"location": "Tokyo"}').with_channel("commentary").with_recipient("functions.get_current_weather").with_content_type("<|constrain|> json")

#     # tool response addition
#     adapter.add_message_from_dict({
#         "role": "tool",
#         "name": 'functions.get_current_weather',
#         "content": [
#             {
#                 "type": "text",
#                 "text": '{"temperature": "20", "sunny": true }'
#             }
#         ],
#         "channel": "commentary",
#     })

#     test_ext_message = Message.from_author_and_content(
#             Author.new(Role.TOOL, "functions.get_current_weather"),
#             '{ "temperature": 20, "sunny": true }',
#         ).with_channel("commentary")
    
#     adapter.conversations.messages.append(test_ext_message)

#     # prompt = adapter.request_query("How you doing ?")
#     # print("prompt >>", prompt)
#     print("conversations >>", adapter.decode_conversations())
#     for conv in adapter.conversations.messages:
#         print(conv.to_dict())
#     # print("raw conversations >>", adapter.conversations)