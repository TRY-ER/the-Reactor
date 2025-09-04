# API Endpoint Documentation

This API provides project management functionality using FastAPI and SQLite.

## Getting Started

1. Install dependencies:
   ```bash
   pip install -r requirements.txt
   ```
2. Run the server:
   ```bash
   uvicorn main:app --reload
   ```
3. The API will be available at `http://localhost:8000/`

---

## Environment Variables
- `API_KEY`: Your API key (not used in endpoints yet)
- `DEBUG`: Set to `True` for debug mode
- `ALLOWED_ORIGINS`: Comma-separated list of allowed CORS origins
- `DATABASE_URL`: (optional) SQLite DB URL, default is `sqlite:///./projects.db`

---

## Endpoints

### Dummy Route
- **GET /**
  - Returns: `{ "Hello": "World" }`

- **GET /hello**
  - Returns: `{ "message": "Hello from router!" }`


### Project CRUD

- **POST /projects/**
  - Create a new project
  - Request body (JSON):
    ```json
    {
      "name": "Project Name",
      "description": "Project description",
      "sessions": ["session_id_1", "session_id_2"]
    }
    ```
  - Response: Project object with `id`, `name`, `description`, `sessions`

- **GET /projects/**
  - List all projects (supports `skip` and `limit` query params)
  - Response: List of project objects, each with `id`, `name`, `description`, `sessions`

- **GET /projects/{project_id}**
  - Get a project by ID
  - Response: Project object with `id`, `name`, `description`, `sessions`
  - 404 if not found

- **PUT /projects/{project_id}**
  - Update a project by ID
  - Request body (JSON):
    ```json
    {
      "name": "New Name",
      "description": "New description",
      "sessions": ["session_id_1", "session_id_2"]
    }
    ```
  - Response: Updated project object with `id`, `name`, `description`, `sessions`
  - 404 if not found

- **DELETE /projects/{project_id}**
  - Delete a project by ID
  - Response: `{ "ok": true }`
  - 404 if not found

---

### Session CRUD

- **POST /sessions/**
  - Create a new session
  - Request body (JSON):
    ```json
    {
      "name": "My Session",
      "project_id": "123e4567-e89b-12d3-a456-426614174000"
    }
    ```
  - Response: Session object with `id`, `name`, `content`, `state`, `last_updated`, `worker_id`, `project_id`
  - Note: `id` is automatically generated, `worker_id` is assigned randomly, `content` and `state` are initially null. If no state is provided, it defaults to "active". **project_id is required** - you must create a project first.

- **GET /sessions/**
  - List all sessions (supports `skip` and `limit` query params)
  - Response: List of session objects, each with `id`, `name`, `content`, `state`, `last_updated`, `worker_id`

- **GET /sessions/{session_id}**
  - Get a session by ID
  - Response: Session object with `id`, `name`, `content`, `state`, `last_updated`, `worker_id`
  - 404 if not found

- **PUT /sessions/{session_id}**
  - Update a session by ID
  - Request body (JSON):
    ```json
    {
      "content": "{\"key\": \"new_value\"}",
      "state": "updated",
      "worker_id": "worker456"
    }
    ```
  - Response: Updated session object with `id`, `name`, `content`, `state`, `last_updated`, `worker_id`
  - Note: Only `content`, `state`, and `worker_id` can be updated; `name` remains unchanged
  - 404 if not found

- **DELETE /sessions/{session_id}**
  - Delete a session by ID
  - Response: `{ "ok": true }`
  - 404 if not found

---

## Testing

Run all tests:
```bash
pytest
```

---


## Notes
- The `state` field has been removed from the project schema. Projects now have a `sessions` field, which is a list of session IDs (strings) associated with the project.
- The API now supports a full Session object with fields: `id`, `name`, `content` (JSON as string), `state`, `last_updated`, and `worker_id`. Session CRUD endpoints are available.
- To associate sessions with a project, add their IDs to the project's `sessions` list.
- CORS is enabled for origins set in `.env`.
- Database is SQLite by default, file `projects.db` will be created in the project directory.
- All endpoints return JSON responses.
