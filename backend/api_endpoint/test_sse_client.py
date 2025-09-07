#!/usr/bin/env python3
"""
Test SSE Client for Session Content Streaming

This script demonstrates how to connect to the SSE endpoint and receive
real-time updates from a running session.
"""

import asyncio
import aiohttp
import json
import uuid
from datetime import datetime


async def test_sse_stream(session_id: str, base_url: str = "http://localhost:8000"):
    """Test SSE streaming for a session"""
    client_id = str(uuid.uuid4())
    url = f"{base_url}/sessions/{session_id}/stream?client_id={client_id}"
    
    print(f"Connecting to SSE stream: {url}")
    print(f"Client ID: {client_id}")
    print(f"Timestamp: {datetime.now().isoformat()}")
    print("-" * 50)
    
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as response:
                if response.status != 200:
                    print(f"Error: HTTP {response.status}")
                    print(f"Response: {await response.text()}")
                    return
                
                print("Connected to SSE stream. Listening for events...")
                print("-" * 50)
                
                async for line in response.content:
                    line = line.decode('utf-8').strip()
                    
                    if line.startswith('data: '):
                        data_str = line[6:]  # Remove 'data: ' prefix
                        try:
                            data = json.loads(data_str)
                            timestamp = data.get('timestamp', datetime.now().isoformat())
                            
                            print(f"[{timestamp}] ", end="")
                            
                            if data.get('_type') == 'completion':
                                print(f"🏁 COMPLETION: {data.get('status')}")
                                if data.get('error'):
                                    print(f"   Error: {data.get('error')}")
                                break
                            elif data.get('type') == 'info':
                                print(f"ℹ️  INFO: {data.get('message', data.get('content', ''))}")
                                if 'aux' in data:
                                    aux = data['aux']
                                    if aux.get('info_type'):
                                        print(f"   Type: {aux['info_type']}")
                                    if aux.get('reaction_index') is not None:
                                        print(f"   Reaction: {aux['reaction_index']}")
                            elif data.get('type') == 'data':
                                print(f"📊 DATA: {data.get('message', '')}")
                                if 'aux' in data:
                                    aux = data['aux']
                                    if aux.get('returnable'):
                                        returnable = aux['returnable']
                                        print(f"   Values: {len(returnable.get('values', []))} items")
                                        print(f"   Invalid count: {returnable.get('invalid_count', 0)}")
                            elif data.get('type') == 'error':
                                print(f"❌ ERROR: {data.get('message', data.get('content', ''))}")
                            elif 'final_csv' in data:
                                print(f"✅ FINAL CSV: {len(data['final_csv'])} reactions processed")
                            else:
                                print(f"📝 RAW: {data}")
                                
                        except json.JSONDecodeError as e:
                            print(f"❌ JSON Error: {e}")
                            print(f"   Raw data: {data_str}")
                    
                    elif line:
                        print(f"📝 Non-data line: {line}")
                
                print("-" * 50)
                print("Stream ended.")
                
    except Exception as e:
        print(f"❌ Connection error: {e}")


async def get_streaming_status(session_id: str, base_url: str = "http://localhost:8000"):
    """Get streaming status for a session"""
    url = f"{base_url}/sessions/{session_id}/streaming-status"
    
    try:
        async with aiohttp.ClientSession() as session:
            async with session.get(url) as response:
                if response.status == 200:
                    data = await response.json()
                    print("Streaming Status:")
                    print(f"  Has active worker: {data.get('has_active_worker')}")
                    print(f"  Status: {data.get('status')}")
                    print(f"  Queue size: {data.get('queue_size')}")
                    print(f"  Is complete: {data.get('is_complete')}")
                    print(f"  Active clients: {data.get('streaming_clients', [])}")
                    if data.get('worker_id'):
                        print(f"  Worker ID: {data.get('worker_id')}")
                    return data
                else:
                    print(f"Error getting status: HTTP {response.status}")
                    return None
    except Exception as e:
        print(f"Error getting status: {e}")
        return None


async def main():
    """Main function to test SSE streaming"""
    if len(sys.argv) < 2:
        print("Usage: python test_sse_client.py <session_id> [base_url]")
        print("Example: python test_sse_client.py abc123 http://localhost:8000")
        return
    
    session_id = sys.argv[1]
    base_url = sys.argv[2] if len(sys.argv) > 2 else "http://localhost:8000"
    
    print(f"Testing SSE for session: {session_id}")
    
    # First check the streaming status
    status = await get_streaming_status(session_id, base_url)
    if status:
        print("-" * 50)
    
    # Then start streaming
    await test_sse_stream(session_id, base_url)


if __name__ == "__main__":
    import sys
    asyncio.run(main())
