import sys
import os
import pytest
from fastapi.testclient import TestClient

# Ensure the app can be imported
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from main import app

client = TestClient(app)

def test_hello_route():
    response = client.get("/hello")
    assert response.status_code == 200
    assert response.json() == {"message": "Hello from router!"}
