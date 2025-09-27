import os
import sys
# sys.path.append(os.path.join(os.path.dirname(os.path.dirname(__file__))))
sys.path.append(os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))))
print("path >>",os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(__file__))))))

from app.agents.scrapper_agent.config import MODEL_NAME, ENDPOINT_TYPE
import asyncio
from agents import Agent, Runner, FunctionTool, RunContextWrapper, function_tool, WebSearchTool
from agents.extensions.models.litellm_model import LitellmModel
from agents._run_impl import QueueCompleteSentinel
from openai.types.responses import (
    ResponseTextDeltaEvent,
)
from agents.stream_events import (
    RawResponsesStreamEvent,
)

from app.agents.scrapper.master_scraper import scrape_website

@function_tool
def scrape_website_tool(url: str):
    scrape_result = scrape_website(url)
    file_names = scrape_result.get("files_created", [])
    # get file name and extract txt file and image file separately
    # use the file path to return those as different response type available
    

class OneshotScrapperAgent:
    """
    This class implements one shot agent with single pass over an LLM. 
    """

    def __init__(self,  model_name: str, endpoint_type: str, agent_name="scrapper_agent",instructions: str = "You are an agent that can browse and scrapes websites"):
        self.model_name = model_name
        self.endpoint_type = endpoint_type
        self.agent = Agent(
            name=agent_name,
            instructions=instructions,
            model=self.get_model(),
            tools= [
                WebSearchTool()
            ]
        )

    def get_model(self):
        try:
            if self.endpoint_type == "groq":
                API_KEY = os.getenv("GROQ_API_KEY")
                assert API_KEY is not None, "GROQ_API_KEY environment variable is not set."
                return LitellmModel(model=self.model_name, api_key=API_KEY)
            else:
                return self.model_name 
        except Exception as e:
            print(f"Error occurred while getting model: {e}")
            return None

    async def run(self, query: str, *args, **kwargs):
        result = await Runner.run(starting_agent=self.agent, input=query, *args, **kwargs)
        return result.final_output

    async def run_streamed(self, query, *args, **kwargs):
        streaming_runner = Runner.run_streamed(
            starting_agent=self.agent, input=query, *args, **kwargs)
        queue = streaming_runner._event_queue
        while not streaming_runner.is_complete:
            try:
                event = await asyncio.wait_for(queue.get(), timeout=1.0)
                if isinstance(event, QueueCompleteSentinel):
                    break
                if isinstance(event, RawResponsesStreamEvent):
                    if isinstance(event.data, ResponseTextDeltaEvent):
                        yield event.data
            except asyncio.QueueEmpty:
                continue


def run_agent(query, model=MODEL_NAME, endpoint_type=ENDPOINT_TYPE, *args, **kwargs):
    """Run the one-shot agent with the given query."""
    agent_instance = OneshotScrapperAgent(
        model_name=model, endpoint_type=endpoint_type)
    response = asyncio.run(agent_instance.run(query, *args, **kwargs))
    return response


if __name__ == "__main__":
    query = "What can you find about websites on organization named xonext ?"
    response = run_agent(query)
    print(response)
