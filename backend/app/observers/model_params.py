import dotenv
import os
import sys

sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../../..')))

dotenv.load_dotenv()


MODEL_PARAMS = {
    "groq": {
        "api_key": os.getenv("GROQ_API_KEY"),
        "models": [
            # "groq/moonshotai/kimi-k2-instruct",
            "meta-llama/llama-4-maverick-17b-128e-instruct",
            # "groq/openai/gpt-oss-120b",
            # "groq/openai/gpt-oss-20b",
        ]
    },
    # "openai": {
    #     "api_key": os.getenv("OPENAI_API_KEY"),
    #     "models": [
    #         "gpt-5-2025-08-07",
    #     ]
    # },
    # "xai": {
    #     "api_key": os.getenv("XAI_API_KEY"),
    #     "models": [
    #         "xai/grok-4-0709",
    #     ]
    # },
    # "google": {
    #     "api_key": os.getenv("GEMINI_API_KEY"),
    #     "models": [
    #         "gemini/gemini-2.5-pro",
    #     ]
    # },
    # "anthropic": {
    #     "api_key": os.getenv("ANTHROPIC_API_KEY"),
    #     "models": [
    #         "claude-opus-4-20250514",
    #     ]
    # }
}


MODEL_PROGRESS_PARAMS = {
    "max_model_type": len(MODEL_PARAMS),
}

MODEL_LIMITER_MAP = {
    "groq/openai/gpt-oss-120b": {
        "limit": 30,
        "time": 60
    },
    "groq/openai/gpt-oss-20b": {
        "limit": 30,
        "time": 60 
    },
}


