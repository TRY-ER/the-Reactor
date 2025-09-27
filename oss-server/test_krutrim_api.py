import os
from krutrim_cloud import KrutrimCloud
from dotenv import load_dotenv

# Load environment variables from .env file
load_dotenv()

# Initialize the Krutrim Cloud client
client = KrutrimCloud()

# Check if API key is available
api_key = os.getenv("KRUTRIM_CLOUD_API_KEY")
if not api_key:
    print("Error: KRUTRIM_CLOUD_API_KEY environment variable not found!")
    exit(1)

print(f"Using API key: {api_key[:10]}..." if len(api_key) > 10 else api_key)

# Define the model and messages
model_name = "gpt-oss-20b"
messages = [{"role": "user", "content": "how to get to Mars?"}]

print("\n=== Testing Streaming Response ===")
try:
    print("Starting stream...")
    
    stream = client.chat.completions.create(
        model=model_name,
        messages=messages,
        stream=True
    )
    
    print("Response: ", end='', flush=True)
    
    # Process each chunk as it arrives
    for chunk in stream:
        # Extract content from chunk (no debug output, just content)
        if type(chunk) is str:
            # print(chunk, end='', flush=True)
            if (chunk.startswith("data: ")):
                print("this if is getting triggerd !")
                dict_str = chunk[6:]
                if dict_str == "[DONE]":
                    break
                else:
                    print('dict_str >>', dict_str)
            else:
                print(chunk, end='', flush=True)
        # print("type of chunk >>",type(chunk))
        content = ""
        
        if hasattr(chunk, 'choices') and chunk.choices:
            choice = chunk.choices[0]
            
            # Check for delta content (streaming pattern)
            if hasattr(choice, 'delta') and choice.delta and hasattr(choice.delta, 'content'):
                content = choice.delta.content or ""
            
            # Check for message content (non-streaming pattern)
            elif hasattr(choice, 'message') and choice.message and hasattr(choice.message, 'content'):
                content = choice.message.content or ""
        
        # Print content exactly as received (preserve backend formatting)
        # if content:
        #     print(content, end='', flush=True)
    
    print("\n✅ Streaming completed!")
    
except Exception as exc:
    print(f"\n❌ Streaming failed: {exc}")

# print("\n=== Testing Non-Streaming Response ===")
# try:
#     print("Making regular API request...")
    
#     response = client.chat.completions.create(
#         model=model_name,
#         messages=messages,
#         stream=False
#     )
    
#     if hasattr(response, 'choices') and response.choices:
#         content = response.choices[0].message.content
#         print(f"Response: {content}")
#         print("✅ Regular API call successful!")
    
# except Exception as exc:
#     print(f"❌ Regular request failed: {exc}")
