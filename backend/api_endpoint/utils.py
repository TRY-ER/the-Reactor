from dotenv import load_dotenv
load_dotenv()

ALLOWED_MODEL_TYPES = [ 
    "groq"
]

def get_model_key(model_type: str):
    """Map model_type to model_key"""
    if model_type not in ALLOWED_MODEL_TYPES:
        raise ValueError(f"Unsupported model_type: {model_type}. Supported types: {ALLOWED_MODEL_TYPES}")
    if model_type == "groq":
        import os
        model_key = os.getenv("GROQ_API_KEY", "").strip()
        if not model_key:
            raise ValueError("GROQ_API_KEY environment variable is not set or empty")
        return model_key
