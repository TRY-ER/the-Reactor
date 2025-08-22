import os
import sys

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../../..')))

import dotenv
dotenv.load_dotenv()



MODEL_PARAMS = {
    "groq": {
        "api_key": os.getenv("GROQ_API_KEY"),
        "models": [
            "groq/openai/gpt-oss-120b",
        ]
    },
    "openai": {
        "api_key": os.getenv("OPENAI_API_KEY"),
        "models": [
            "gpt-5-2025-08-07",
        ]
    }
}
