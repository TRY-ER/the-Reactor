import os
from krutrim_cloud import KrutrimCloud
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Initialize client
client = KrutrimCloud()

# Check API key
api_key = os.getenv("KRUTRIM_CLOUD_API_KEY")
if not api_key:
    print("Error: KRUTRIM_CLOUD_API_KEY environment variable not found!")
    exit(1)

print(f"Using API key: {api_key[:10]}...")

# Define the model and messages
model_name = "gpt-oss-20b"
messages = [{"role": "user", "content": "Explain how to get to Mars in 2-3 sentences."}]

print("\n🌊 STREAMING RESPONSE:")
print("=" * 50)

try:
    # Create streaming response
    stream = client.chat.completions.create(
        model=model_name,
        messages=messages,
        stream=True
    )
    
    # Process stream chunks
    for chunk in stream:
        # Extract content from chunk
        content = ""
        
        # Try different content access patterns
        if hasattr(chunk, 'choices') and chunk.choices:
            choice = chunk.choices[0]
            
            # Check for delta content (typical for streaming)
            if hasattr(choice, 'delta') and choice.delta and hasattr(choice.delta, 'content'):
                content = choice.delta.content or ""
            
            # Check for message content (fallback)
            elif hasattr(choice, 'message') and choice.message and hasattr(choice.message, 'content'):
                content = choice.message.content or ""
        
        # Print content as-is (no extra newlines unless from backend)
        if content:
            print(content, end='', flush=True)
    
    print("\n")  # Single newline at the end
    print("=" * 50)
    print("✅ Streaming completed successfully!")
    
except Exception as e:
    print(f"\n❌ Streaming failed: {e}")
    import traceback
    traceback.print_exc()
