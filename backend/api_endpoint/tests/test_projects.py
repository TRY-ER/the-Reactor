import sys
import os
import pytest
from fastapi.testclient import TestClient

# Ensure the app can be imported
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from main import app
from db import SessionLocal, Base, engine

client = TestClient(app)

@pytest.fixture(autouse=True)
def setup_and_teardown_db():
    # Recreate tables before each test
    Base.metadata.drop_all(bind=engine)
    Base.metadata.create_all(bind=engine)
    yield
    Base.metadata.drop_all(bind=engine)


def test_create_project():
    response = client.post("/projects/", json={"name": "Test Project", "description": "A test project.", "sessions": ["sess1"]})
    assert response.status_code == 201  # Changed from 200 to 201 for creation
    data = response.json()
    assert data["name"] == "Test Project"
    assert data["description"] == "A test project."
    assert data["sessions"] == ["sess1"]
    assert "id" in data


def test_read_projects():
    # Create a project first
    client.post("/projects/", json={"name": "Proj1", "description": "Desc1", "sessions": ["sessA"]})
    client.post("/projects/", json={"name": "Proj2", "description": "Desc2", "sessions": ["sessB"]})
    response = client.get("/projects/")
    assert response.status_code == 200
    data = response.json()
    assert len(data) == 2
    assert data[0]["name"] == "Proj1"
    assert data[0]["sessions"] == ["sessA"]
    assert data[1]["name"] == "Proj2"
    assert data[1]["sessions"] == ["sessB"]


def test_read_project():
    # Create a project
    create_resp = client.post("/projects/", json={"name": "ProjX", "description": "DescX", "sessions": ["sessX"]})
    proj_id = create_resp.json()["id"]
    response = client.get(f"/projects/{proj_id}")
    assert response.status_code == 200
    data = response.json()
    assert data["name"] == "ProjX"
    assert data["description"] == "DescX"
    assert data["sessions"] == ["sessX"]


def test_update_project():
    # Create a project
    create_resp = client.post("/projects/", json={"name": "ProjY", "description": "DescY", "sessions": ["sessY"]})
    proj_id = create_resp.json()["id"]
    # Update it
    response = client.put(f"/projects/{proj_id}", json={"name": "ProjY-updated", "description": "DescY-updated", "sessions": ["sessY2"]})
    assert response.status_code == 200
    data = response.json()
    assert data["name"] == "ProjY-updated"
    assert data["description"] == "DescY-updated"
    assert data["sessions"] == ["sessY2"]


def test_delete_project():
    # Create a project
    create_resp = client.post("/projects/", json={"name": "ProjZ", "description": "DescZ", "sessions": ["sessZ"]})
    proj_id = create_resp.json()["id"]
    # Delete it
    response = client.delete(f"/projects/{proj_id}")
    assert response.status_code == 200
    assert response.json() == {"ok": True}
    # Ensure it's gone
    get_resp = client.get(f"/projects/{proj_id}")
    assert get_resp.status_code == 404
