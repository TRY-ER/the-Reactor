import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

ALLOWED_MODELS = ["openai/gpt-oss-120b", 'openai/gpt-oss-20b']

async def run_agent(
    model: str,
    query: str,
    api_key: str
    ):
