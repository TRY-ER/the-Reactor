#!/usr/bin/env python3

import os
import sys
from dotenv import load_dotenv

# Load environment variables
load_dotenv()

# Add debugging
print("Python version:", sys.version)
print("Starting streaming test...")

try:
    from krutrim_cloud import KrutrimCloud
    print("✅ Krutrim Cloud SDK imported")
    
    # Initialize client
    client = KrutrimCloud()
    print("✅ Client initialized")
    
    # Check API key
    api_key = os.getenv("KRUTRIM_CLOUD_API_KEY")
    if not api_key:
        print("❌ No API key found")
        sys.exit(1)
    
    print(f"✅ API key found: {api_key[:10]}...")
    
    # Test messages
    messages = [{"role": "user", "content": "Explain Mars in one sentence."}]
    model_name = "gpt-oss-20b"
    
    print("\n🔄 Testing NON-STREAMING response...")
    try:
        response = client.chat.completions.create(
            model=model_name,
            messages=messages,
            stream=False
        )
        print(f"✅ Non-streaming response: {type(response)}")
        if hasattr(response, 'choices') and response.choices:
            content = response.choices[0].message.content
            print(f"📄 Content: {content}")
    except Exception as e:
        print(f"❌ Non-streaming failed: {e}")
    
    print("\n🌊 Testing STREAMING response...")
    try:
        stream_response = client.chat.completions.create(
            model=model_name,
            messages=messages,
            stream=True
        )
        print(f"✅ Stream response: {type(stream_response)}")
        
        print("📡 Starting to stream...")
        content_buffer = ""
        chunk_count = 0
        
        for chunk in stream_response:
            chunk_count += 1
            print(f"\n🔸 Chunk {chunk_count}:")
            print(f"   Type: {type(chunk)}")
            print(f"   Raw: {chunk}")
            
            # Examine the chunk structure
            if hasattr(chunk, '__dict__'):
                print(f"   Attributes: {list(chunk.__dict__.keys())}")
            
            # Try different access patterns
            chunk_content = ""
            
            # Pattern 1: Standard OpenAI-style streaming
            if hasattr(chunk, 'choices'):
                try:
                    if chunk.choices and len(chunk.choices) > 0:
                        choice = chunk.choices[0]
                        if hasattr(choice, 'delta') and choice.delta:
                            if hasattr(choice.delta, 'content') and choice.delta.content:
                                chunk_content = choice.delta.content
                                print(f"   🎯 Delta content: '{chunk_content}'")
                        elif hasattr(choice, 'message') and choice.message:
                            if hasattr(choice.message, 'content') and choice.message.content:
                                chunk_content = choice.message.content
                                print(f"   🎯 Message content: '{chunk_content}'")
                except Exception as e:
                    print(f"   ⚠️ Error accessing choices: {e}")
            
            # Pattern 2: Direct content access
            if hasattr(chunk, 'content'):
                try:
                    chunk_content = chunk.content
                    print(f"   🎯 Direct content: '{chunk_content}'")
                except Exception as e:
                    print(f"   ⚠️ Error accessing content: {e}")
            
            # Pattern 3: Text field
            if hasattr(chunk, 'text'):
                try:
                    chunk_content = chunk.text
                    print(f"   🎯 Text field: '{chunk_content}'")
                except Exception as e:
                    print(f"   ⚠️ Error accessing text: {e}")
            
            # Add to buffer and display
            if chunk_content:
                content_buffer += chunk_content
                print(f"💬 ", end="", flush=True)
                print(chunk_content, end="", flush=True)
        
        print(f"\n\n✅ Streaming completed!")
        print(f"📊 Total chunks: {chunk_count}")
        print(f"📝 Full content: '{content_buffer}'")
        
    except Exception as e:
        print(f"❌ Streaming failed: {e}")
        import traceback
        traceback.print_exc()
    
except ImportError as e:
    print(f"❌ Import error: {e}")
except Exception as e:
    print(f"❌ General error: {e}")
    import traceback
    traceback.print_exc()

print("\n🏁 Test completed!")
