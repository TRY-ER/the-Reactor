import os
from dotenv import load_dotenv

print("Starting test...")

# Load environment variables
load_dotenv()
print("Environment loaded")

# Check API key
api_key = os.getenv("KRUTRIM_CLOUD_API_KEY")
print(f"API key exists: {bool(api_key)}")
if api_key:
    print(f"API key preview: {api_key[:10]}...")

# Try importing the SDK
try:
    from krutrim_cloud import KrutrimCloud
    print("✅ SDK imported successfully")
    
    # Try creating client
    try:
        client = KrutrimCloud()
        print("✅ Client created successfully")
        print(f"Client type: {type(client)}")
        
        # Try a simple request
        try:
            response = client.chat.completions.create(
                model="gpt-oss-20b",
                messages=[{"role": "user", "content": "Hello!"}]
            )
            print("✅ API call successful!")
            print(f"Response type: {type(response)}")
            
        except Exception as e:
            print(f"❌ API call failed: {e}")
            
    except Exception as e:
        print(f"❌ Client creation failed: {e}")
        
except Exception as e:
    print(f"❌ SDK import failed: {e}")
    import traceback
    traceback.print_exc()

print("Test completed")
