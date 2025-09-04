import json

def format_for_stream(type: str, content: str, aux: dict = {}) -> str:
    content = content.replace("\n", " ").replace("\r", " ")
    return f"data:{json.dumps({'type': type, 'content': content, 'aux' : aux})}\n\n"
