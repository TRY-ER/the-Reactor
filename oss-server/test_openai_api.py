# initiate the environment variables
import os
from dotenv import load_dotenv

load_dotenv()

from agents import Agent, Runner
from agents.extensions.models.litellm_model import LitellmModel
from agents._run_impl import QueueCompleteSentinel
from openai.types.responses import (
    ResponseTextDeltaEvent,
)
from agents.stream_events import (
    RawResponsesStreamEvent,
)

MODEL_NAME = "groq/openai/gpt-oss-120b"

agent = Agent(
    name="test_openai_api",
    instructions="Test OpenAI API with streaming and non-streaming responses.",
    model=LitellmModel(model=MODEL_NAME, api_key=os.getenv("GROQ_API_KEY")),
)

async def run_tests():
    # Test non-streaming response
    print("\n🔄 Testing NON-STREAMING response...")
    result = await Runner.run(agent, "Explain Mars in one sentence.")
    print(f"Non-streaming response: {result}")

    # Test streaming response
    print("\n🌊 Testing STREAMING response...")
    streaming_result = Runner.run_streamed(agent, "write me a poem about the sea")
    queue = streaming_result._event_queue
    while not streaming_result.is_complete:
        try:
            # Get events from queue with timeout
            event = await asyncio.wait_for(queue.get(), timeout=1.0)
            
            if isinstance(event, QueueCompleteSentinel):
                break
                
            # Process the event
            if isinstance(event, RawResponsesStreamEvent):
                if isinstance(event.data, ResponseTextDeltaEvent):
                    print(event.data.delta, end='', flush=True)
            # print(f"Event: {event}")
            
        except asyncio.TimeoutError:
            continue
    
    # print(f"Streaming complete: {streaming_result.final_output}")
    # for chunk in result:
    #     print(chunk, end='', flush=True)

# Run the tests
if __name__ == "__main__":
    import asyncio
    asyncio.run(run_tests())
    print("\n✅ All tests completed successfully!")