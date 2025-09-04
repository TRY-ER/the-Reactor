import sys
import os
import pytest
from fastapi.testclient import TestClient

# Ensure the app can be imported
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from db import Base, engine
from main import app


client = TestClient(app)

@pytest.fixture(autouse=True)
def setup_and_teardown_db():
    # Recreate tables before each test
    Base.metadata.drop_all(bind=engine)
    Base.metadata.create_all(bind=engine)
    yield
    Base.metadata.drop_all(bind=engine)


def test_create_session():
    # First create a project to get a project_id
    project_response = client.post("/projects/", json={"name": "Test Project", "description": "Test"})
    project_id = project_response.json()["id"]
    
    session_data = {
        "name": "sess1",
        "project_id": project_id
    }
    response = client.post("/sessions/", json=session_data)
    assert response.status_code == 200
    data = response.json()
    assert data["name"] == session_data["name"]
    assert data["project_id"] == project_id
    assert data["id"] is not None
    assert data["worker_id"] is not None
    assert data["content"] is None
    assert data["state"] == "active"  # Now defaults to "active" instead of None
    assert "last_updated" in data
    assert "file_paths" in data
    assert data["file_paths"] == []


def test_read_sessions():
    # First create a project to get a project_id
    project_response = client.post("/projects/", json={"name": "Test Project 2", "description": "Test"})
    project_id = project_response.json()["id"]
    
    # Create two sessions
    client.post("/sessions/", json={"name": "sess2", "project_id": project_id})
    client.post("/sessions/", json={"name": "sess3", "project_id": project_id})
    response = client.get("/sessions/")
    assert response.status_code == 200
    data = response.json()
    assert len(data) >= 2
    names = [s["name"] for s in data]
    assert "sess2" in names and "sess3" in names


def test_update_session():
    # First create a project to get a project_id
    project_response = client.post("/projects/", json={"name": "Test Project 3", "description": "Test"})
    project_id = project_response.json()["id"]
    
    # Create a session
    client.post("/sessions/", json={"name": "sess4", "project_id": project_id})
    # Get the session to find its id
    response = client.get("/sessions/")
    data = response.json()
    sess4 = next(s for s in data if s["name"] == "sess4")
    sess_id = sess4["id"]
    # Update it
    update_data = {
        "name": "sess4", 
        "content": "{\"updated\": true}", 
        "state": "active", 
        "worker_id": "w2",
        "project_id": project_id,
        "file_paths": ["test.txt"]
    }
    response = client.put(f"/sessions/{sess_id}", json=update_data)
    assert response.status_code == 200
    data = response.json()
    assert data["content"] == update_data["content"]
    assert data["state"] == update_data["state"]
    assert data["worker_id"] == update_data["worker_id"]
    assert data["name"] == update_data["name"]
    assert data["project_id"] == project_id
    assert data["file_paths"] == ["test.txt"]


def test_delete_session():
    # First create a project to get a project_id
    project_response = client.post("/projects/", json={"name": "Test Project 4", "description": "Test"})
    project_id = project_response.json()["id"]
    
    # Create a session
    client.post("/sessions/", json={"name": "sess5", "project_id": project_id})
    # Get the session to find its id
    response = client.get("/sessions/")
    data = response.json()
    sess5 = next(s for s in data if s["name"] == "sess5")
    sess_id = sess5["id"]
    # Delete it
    response = client.delete(f"/sessions/{sess_id}")
    assert response.status_code == 200
    assert response.json() == {"ok": True}
    # Ensure it's gone
    get_resp = client.get(f"/sessions/{sess_id}")
    assert get_resp.status_code == 404
