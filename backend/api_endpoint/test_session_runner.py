#!/usr/bin/env python3
"""
Test script for session runner functionality
"""

import asyncio
import json
import requests
import time
from datetime import datetime

BASE_URL = "http://localhost:8000"

def create_test_project():
    """Create a test project"""
    project_data = {
        "name": "Test Project for Session Runner",
        "description": "A test project to demonstrate session runner functionality",
        "sessions": []
    }
    
    response = requests.post(f"{BASE_URL}/projects/", json=project_data)
    if response.status_code == 200:
        project = response.json()
        print(f"Created project: {project['id']}")
        return project['id']
    else:
        print(f"Failed to create project: {response.text}")
        return None

def create_test_session(project_id):
    """Create a test session with required project_id"""
    session_data = {
        "name": "Test Session",
        "project_id": project_id,
        "content": [],
        "file_paths": [],
        "query": ""  # Initially empty query
    }
    
    response = requests.post(f"{BASE_URL}/sessions/", json=session_data)
    if response.status_code == 201:  # Updated to 201 for creation
        session = response.json()
        print(f"Created session: {session['id']}")
        print(f"Using project: {project_id}")
        return session['id']
    else:
        print(f"Failed to create session: {response.text}")
        return None

def update_session_query(session_id, query):
    """Update the query for a session"""
    query_data = {"query": query}
    
    response = requests.put(f"{BASE_URL}/sessions/{session_id}/query", json=query_data)
    if response.status_code == 200:
        result = response.json()
        print(f"Updated query for session {session_id}")
        return True
    else:
        print(f"Failed to update query: {response.text}")
        return False

def get_session_query(session_id):
    """Get the current query for a session"""
    response = requests.get(f"{BASE_URL}/sessions/{session_id}/query")
    if response.status_code == 200:
        result = response.json()
        print(f"Current query: {result['query'][:100]}...")
        return result['query']
    else:
        print(f"Failed to get query: {response.text}")
        return None

def start_worker(session_id, with_query=None):
    """Start a worker for the session"""
    worker_data = {
        "model_type": "groq",
        "model_name": "groq/openai/gpt-oss-120b",
        "model_key": "your-api-key-here",  # Replace with actual API key
        "agent_name": "master_openai_oss_agent"
    }
    
    # Add query only if provided
    if with_query:
        worker_data["query"] = with_query
    
    response = requests.post(f"{BASE_URL}/sessions/{session_id}/start-worker", json=worker_data)
    if response.status_code == 200:
        result = response.json()
        print(f"Started worker: {result['worker_id']}")
        return result['worker_id']
    else:
        print(f"Failed to start worker: {response.text}")
        return None

def rerun_with_new_query(session_id, new_query):
    """Rerun session with a new query"""
    worker_data = {
        "model_type": "groq",
        "model_name": "groq/openai/gpt-oss-120b",
        "model_key": "your-api-key-here",  # Replace with actual API key
        "query": new_query,
        "agent_name": "master_openai_oss_agent"
    }
    
    response = requests.post(f"{BASE_URL}/sessions/{session_id}/rerun-with-query", json=worker_data)
    if response.status_code == 200:
        result = response.json()
        print(f"Restarted session with new query, worker: {result['worker_id']}")
        return result['worker_id']
    else:
        print(f"Failed to rerun session: {response.text}")
        return None

def check_worker_status(session_id):
    """Check worker status"""
    response = requests.get(f"{BASE_URL}/sessions/{session_id}/worker-status")
    if response.status_code == 200:
        status = response.json()
        print(f"Worker status: {status}")
        return status
    else:
        print(f"Failed to get worker status: {response.text}")
        return None

def get_session_content(session_id):
    """Get session content"""
    response = requests.get(f"{BASE_URL}/sessions/{session_id}/content")
    if response.status_code == 200:
        content = response.json()
        print(f"Session content has {len(content['content'])} items")
        return content['content']
    else:
        print(f"Failed to get session content: {response.text}")
        return None

def monitor_session_progress(session_id, max_wait_time=300):
    """Monitor session progress"""
    start_time = time.time()
    
    while time.time() - start_time < max_wait_time:
        status = check_worker_status(session_id)
        if not status:
            print("Worker not found, might have completed or errored")
            break
        
        if status['status'] in ['completed', 'error']:
            print(f"Worker finished with status: {status['status']}")
            break
        
        # Get current content
        content = get_session_content(session_id)
        if content:
            print(f"Progress: {len(content)} results so far")
        
        time.sleep(5)  # Wait 5 seconds before checking again
    
    # Get final content
    print("\nFinal session content:")
    final_content = get_session_content(session_id)
    if final_content:
        for i, item in enumerate(final_content[-5:]):  # Show last 5 items
            print(f"Item {i}: {json.dumps(item, indent=2)[:200]}...")

def main():
    """Main test function"""
    print("Testing Session Runner with Query Management...")
    
    # Create test project and session
    project_id = create_test_project()
    if not project_id:
        return
    
    session_id = create_test_session(project_id)
    if not session_id:
        return
    
    # Test 1: Set initial query
    print("\n=== Test 1: Setting initial query ===")
    initial_query = """
    A synthesis (or combination) reaction involves two or more simple substances combining to form a more complex compound, typically represented as A + B → AB. This type of reaction is characterized by multiple reactants yielding a single product. An example is the formation of iron(II) sulfide from iron and sulfur:

    8 Fe + S₈ → 8 FeS

    Another example includes the reaction of hydrogen and oxygen gases to produce water. Synthesis reactions are vital for forming complex molecules from simpler ones.
    """
    
    if not update_session_query(session_id, initial_query):
        return
    
    # Verify query was set
    get_session_query(session_id)
    
    # Test 2: Start worker using stored query
    print("\n=== Test 2: Starting worker with stored query ===")
    worker_id = start_worker(session_id)  # No query provided - should use stored query
    if not worker_id:
        return
    
    # Monitor progress
    print("\n=== Monitoring first run ===")
    monitor_session_progress(session_id, max_wait_time=60)
    
    # Test 3: Update query and rerun
    print("\n=== Test 3: Updating query and rerunning ===")
    new_query = """
    Decomposition reactions involve a single compound breaking down into two or more simpler substances, typically represented as AB → A + B. This is the opposite of synthesis reactions. An example is the decomposition of hydrogen peroxide:

    2 H₂O₂ → 2 H₂O + O₂

    Another example is the thermal decomposition of calcium carbonate to produce calcium oxide and carbon dioxide.
    """
    
    worker_id = rerun_with_new_query(session_id, new_query)
    if not worker_id:
        return
    
    # Monitor second run
    print("\n=== Monitoring second run ===")
    monitor_session_progress(session_id, max_wait_time=60)
    
    # Test 4: Start worker with inline query (should override stored query)
    print("\n=== Test 4: Starting worker with inline query ===")
    inline_query = """
    Single replacement reactions involve one element replacing another element in a compound, typically represented as A + BC → AC + B. An example is zinc replacing copper in copper sulfate:

    Zn + CuSO₄ → ZnSO₄ + Cu

    The more reactive metal (zinc) displaces the less reactive metal (copper) from its compound.
    """
    
    worker_id = start_worker(session_id, with_query=inline_query)
    if not worker_id:
        return
    
    # Monitor third run
    print("\n=== Monitoring third run ===")
    monitor_session_progress(session_id, max_wait_time=60)
    
    # Verify final query state
    print("\n=== Final query state ===")
    get_session_query(session_id)
    
    print("\nAll tests completed!")

if __name__ == "__main__":
    main()
