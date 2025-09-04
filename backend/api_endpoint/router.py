from fastapi import APIRouter, HTTPException, Depends, status
from sqlalchemy.orm import Session
from datetime import datetime

from db import Project, SessionLocal, SessionObj
from models import SessionBase, SessionCreate, SessionRead, ProjectBase, ProjectCreate, ProjectRead
from file_system_handler import file_system_handler
from session_runner import session_runner
import json
import uuid
from utils import get_model_key

HEADERS = [
    "index",
    "reaction_type",
    "reactions_text",
    "reactions_composition_SMILES",
    "invalid_SMILE",
    "reactions_SMILES",
    "invalid_reaction_SMILES",
    "reactions_SMARTS",
    "invalid_reaction_SMARTS",
    "reactions_SMIRKS",
    "invalid_reaction_SMIRKS"
]

router = APIRouter(
    tags=["Chemical Reactor API"],
    responses={404: {"description": "Not found"}},
)

# Dependency to get DB session
def get_db():
    db = SessionLocal()
    try:
        yield db
    finally:
        db.close()

# CRUD endpoints

# Session CRUD endpoints
@router.post("/sessions/", 
             response_model=SessionRead,
             status_code=status.HTTP_201_CREATED,
             summary="Create a new session",
             description="""
             Create a new session within a project. A session represents a chemical reaction analysis workflow.
             
             **Example Request:**
             ```json
             {
               "name": "Iron Synthesis Analysis",
               "project_id": "123e4567-e89b-12d3-a456-426614174000",
               "query": "8 Fe + S₈ → 8 FeS",
               "content": [],
               "file_paths": []
             }
             ```
             
             **Note:** project_id is required. You must create a project first using the /projects/ endpoint, then use its ID when creating sessions.
             """)
def create_session(session: SessionCreate, db: Session = Depends(get_db)):
    # Verify the project exists
    existing_project = db.query(Project).filter(Project.id == session.project_id).first()
    if not existing_project:
        raise HTTPException(status_code=404, detail=f"Project with id {session.project_id} not found")
    
    # Create the session
    db_session = SessionObj(
        id=str(uuid.uuid4()),  # generate id automatically
        name=session.name,
        content=json.dumps(session.content or []),  # Store as JSON string in DB
        state=session.state or "active",  # Set initial state to "active" if not provided
        last_updated=datetime.utcnow(),
        worker_id=str(uuid.uuid4()),  # assign random worker_id
        project_id=session.project_id,
        file_paths=json.dumps(session.file_paths or []),  # Store as JSON string in DB
        query=session.query or ""  # Store query text
    )
    db.add(db_session)
    try:
        db.commit()
        db.refresh(db_session)
        # Create the session directory using string values
        file_system_handler.create_session_directory(str(session.project_id), str(db_session.id))
    except Exception as e:
        print('error >>', e)
        db.rollback()
        raise HTTPException(status_code=400, detail="Session with this id already exists or error occurred.")
    # Parse JSON strings back to lists for response  
    db_session.content = json.loads(str(db_session.content))
    db_session.file_paths = json.loads(str(db_session.file_paths))
    return db_session

@router.get("/sessions/", 
            response_model=list[SessionRead],
            summary="List all sessions",
            description="""
            Retrieve a paginated list of all sessions.
            
            **Parameters:**
            - **skip**: Number of sessions to skip (for pagination)
            - **limit**: Maximum number of sessions to return
            
            **Example Response:**
            Returns an array of session objects with their current state, content count, and metadata.
            """)
def read_sessions(skip: int = 0, limit: int = 10, db: Session = Depends(get_db)):
    sessions = db.query(SessionObj).offset(skip).limit(limit).all()
    for session in sessions:
        session.content = json.loads(session.content)
        session.file_paths = json.loads(session.file_paths)
    return sessions

@router.get("/sessions/{session_id}", 
            response_model=SessionRead,
            summary="Get a specific session",
            description="""
            Retrieve detailed information about a specific session.
            
            **Path Parameters:**
            - **session_id**: The unique identifier of the session
            
            **Returns:** Complete session data including content, file paths, query, and current state.
            """)
def read_session(session_id: str, db: Session = Depends(get_db)):
    session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
    if not session:
        raise HTTPException(status_code=404, detail="Session not found")
    session.content = json.loads(session.content)
    session.file_paths = json.loads(session.file_paths)
    return session

@router.put("/sessions/{session_id}", 
            response_model=SessionRead,
            summary="Update a session",
            description="""
            Update session information including content, state, and query.
            
            **Path Parameters:**
            - **session_id**: The unique identifier of the session
            
            **Example Request:**
            ```json
            {
              "name": "Updated Session Name",
              "project_id": "project-id",
              "query": "Updated chemical reaction",
              "state": "completed",
              "content": [],
              "file_paths": []
            }
            ```
            """)
def update_session(session_id: str, session: SessionBase, db: Session = Depends(get_db)):
    db_session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
    if not db_session:
        raise HTTPException(status_code=404, detail="Session not found")
    
    # Update all session fields
    db_session.name = session.name  # Update session name
    db_session.content = json.dumps(session.content or [])
    db_session.state = session.state
    # Only update worker_id if it's provided and not None
    if session.worker_id is not None:
        db_session.worker_id = session.worker_id
    db_session.project_id = session.project_id
    db_session.file_paths = json.dumps(session.file_paths or [])
    db_session.query = session.query or ""  # Update query field
    db_session.last_updated = datetime.utcnow()
    
    db.commit()
    db.refresh(db_session)
    # Parse JSON strings back to lists for response
    db_session.content = json.loads(db_session.content)
    db_session.file_paths = json.loads(db_session.file_paths)
    return db_session

@router.delete("/sessions/{session_id}",
               status_code=status.HTTP_200_OK,
               summary="Delete a session",
               description="""
               Permanently delete a session and all its associated data.
               
               **Path Parameters:**
               - **session_id**: The unique identifier of the session to delete
               
               **Warning:** This action cannot be undone. All session content and file paths will be lost.
               """)
def delete_session(session_id: str, db: Session = Depends(get_db)):
    db_session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
    if not db_session:
        raise HTTPException(status_code=404, detail="Session not found")
    db.delete(db_session)
    db.commit()
    return {"ok": True}

from typing import List

# Project schemas are now imported from models.py
# CRUD endpoints

@router.post("/projects/", 
             response_model=ProjectRead,
             status_code=status.HTTP_201_CREATED,
             summary="Create a new project",
             description="""
             Create a new project to organize chemical reaction analysis sessions.
             
             **Example Request:**
             ```json
             {
               "name": "Iron Chemistry Research",
               "description": "Analyzing various iron-based chemical reactions",
               "sessions": []
             }
             ```
             
             **Note:** A master.csv file will be automatically created in the project directory.
             """)
def create_project(project: ProjectCreate, db: Session = Depends(get_db)):
    project_id = str(uuid.uuid4())
    db_project = Project(
        id=project_id,
        name=project.name,
        description=project.description,
        sessions=json.dumps(project.sessions)
    )
    db.add(db_project)
    try:
        db.commit()
        db.refresh(db_project)
        # Create the project directory
        file_system_handler.create_project_directory(project_id)
        file_system_handler.create_master_csv(project_id, headers=HEADERS)
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=400, detail="Project with this name already exists.")
    db_project.sessions = json.loads(db_project.sessions)
    return db_project

@router.get("/projects/", 
            response_model=list[ProjectRead],
            summary="List all projects",
            description="""
            Retrieve a paginated list of all projects.
            
            **Parameters:**
            - **skip**: Number of projects to skip (for pagination)
            - **limit**: Maximum number of projects to return
            """)
def read_projects(skip: int = 0, limit: int = 10, db: Session = Depends(get_db)):
    projects = db.query(Project).offset(skip).limit(limit).all()
    for p in projects:
        p.sessions = json.loads(p.sessions)
    return projects

@router.get("/projects/{project_id}", 
            response_model=ProjectRead,
            summary="Get a specific project",
            description="""
            Retrieve detailed information about a specific project.
            
            **Path Parameters:**
            - **project_id**: The unique identifier of the project
            """)
def read_project(project_id: str, db: Session = Depends(get_db)):
    project = db.query(Project).filter(Project.id == project_id).first()
    if not project:
        raise HTTPException(status_code=404, detail="Project not found")
    project.sessions = json.loads(project.sessions)
    return project


@router.put("/projects/{project_id}", 
            response_model=ProjectRead,
            summary="Update a project",
            description="""
            Update project information including name, description, and associated sessions.
            
            **Path Parameters:**
            - **project_id**: The unique identifier of the project
            """)
def update_project(project_id: str, project: ProjectCreate, db: Session = Depends(get_db)):
    db_project = db.query(Project).filter(Project.id == project_id).first()
    if not db_project:
        raise HTTPException(status_code=404, detail="Project not found")
    db_project.name = project.name
    db_project.description = project.description
    db_project.sessions = json.dumps(project.sessions)
    db.commit()
    db.refresh(db_project)
    db_project.sessions = json.loads(db_project.sessions)
    return db_project

@router.delete("/projects/{project_id}",
               status_code=status.HTTP_200_OK,
               summary="Delete a project",
               description="""
               Permanently delete a project and all its associated files and directories.
               
               **Path Parameters:**
               - **project_id**: The unique identifier of the project to delete
               
               **Warning:** This action cannot be undone. All project data, sessions, and files will be permanently deleted.
               """)
def delete_project(project_id: str, db: Session = Depends(get_db)):
    db_project = db.query(Project).filter(Project.id == project_id).first()
    if not db_project:
        raise HTTPException(status_code=404, detail="Project not found")
    db.delete(db_project)
    db.commit()
    # Delete the project directory
    file_system_handler.delete_project_directory(project_id)
    return {"ok": True}

# New endpoints for session content and file paths management
@router.post("/sessions/{session_id}/content",
             status_code=status.HTTP_201_CREATED,
             summary="Add content to session",
             description="""
             Add a content dictionary to the session's content list. This is typically used by workers to save progress.
             
             **Path Parameters:**
             - **session_id**: The unique identifier of the session
             
             **Request Body:** Any valid JSON object that represents analysis results or progress updates.
             
             **Example Request:**
             ```json
             {
               "type": "reaction_analysis",
               "reaction_index": 0,
               "smiles_data": ["C", "O"],
               "timestamp": "2025-09-04T10:30:00Z"
             }
             ```
             """)
def add_content_to_session(session_id: str, content_dict: dict):
    """Add a content dictionary to the session's content list."""
    try:
        file_system_handler.add_content_to_session(session_id, content_dict)
        return {"message": "Content added successfully", "session_id": session_id}
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error adding content: {str(e)}")

@router.post("/sessions/{session_id}/file-paths",
             status_code=status.HTTP_201_CREATED,
             summary="Add file path to session",
             description="""
             Add a file path dictionary to the session's file_paths list.
             
             **Path Parameters:**
             - **session_id**: The unique identifier of the session
             
             **Example Request:**
             ```json
             {
               "file_type": "csv",
               "file_name": "reaction_results.csv",
               "file_path": "/path/to/file",
               "created_at": "2025-09-04T10:30:00Z"
             }
             ```
             """)
def add_file_path_to_session(session_id: str, file_path_dict: dict):
    """Add a file path dictionary to the session's file_paths list."""
    try:
        file_system_handler.add_file_path_to_session(session_id, file_path_dict)
        return {"message": "File path added successfully", "session_id": session_id}
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error adding file path: {str(e)}")

@router.get("/sessions/{session_id}/content",
            summary="Get session content",
            description="""
            Retrieve all content items for a session. This includes all analysis results, progress updates, and worker outputs.
            
            **Path Parameters:**
            - **session_id**: The unique identifier of the session
            
            **Response:** List of all content objects with timestamps and metadata.
            """)
def get_session_content(session_id: str):
    """Get the content list for a session."""
    try:
        content = file_system_handler.get_session_content(session_id)
        return {"session_id": session_id, "content": content}
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error retrieving content: {str(e)}")

@router.get("/sessions/{session_id}/file-paths",
            summary="Get session file paths",
            description="""
            Retrieve all file paths associated with a session.
            
            **Path Parameters:**
            - **session_id**: The unique identifier of the session
            
            **Response:** List of file path objects with metadata.
            """)
def get_session_file_paths(session_id: str):
    """Get the file paths list for a session."""
    try:
        file_paths = file_system_handler.get_session_file_paths(session_id)
        return {"session_id": session_id, "file_paths": file_paths}
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error retrieving file paths: {str(e)}")

# Health check endpoint
@router.get("/hello",
            summary="Health check",
            description="""
            Simple health check endpoint to verify the API is running.
            
            **Response:**
            ```json
            {
              "message": "Hello from router!"
            }
            ```
            """,
            tags=["Health Check"])
def hello():
    return {"message": "Hello from router!"}

# Session Worker Management Endpoints
from pydantic import BaseModel, Field
from typing import Optional

class StartWorkerRequest(BaseModel):
    """Request model for starting a worker
    
    Example:
    {
        "model_name": "groq/openai/gpt-oss-120b", 
        "query": "8 Fe + S₈ → 8 FeS",
        "agent_name": "master_openai_oss_agent"
    }
    """
    model_name: str  # e.g., "groq/openai/gpt-oss-120b"
    query: Optional[str] = None  # Optional - uses session's stored query if not provided
    agent_name: str = "master_openai_oss_agent"

class UpdateQueryRequest(BaseModel):
    """Request model for updating session query
    
    Example:
    {
        "query": "2 H₂O₂ → 2 H₂O + O₂"
    }
    """
    query: str  # Chemical reaction text

@router.post("/sessions/{session_id}/start-worker",
             status_code=status.HTTP_201_CREATED,
             summary="Start AI worker for session",
             description="""
             Start an asynchronous AI worker to analyze chemical reactions in a session.
             
             **How it works:**
             1. Creates a background worker with unique ID
             2. Uses provided query OR session's stored query
             3. Processes reactions through AI agents (SMILES, SMARTS, SMIRKS analysis)
             4. Saves all results incrementally to the database
             5. Continues even if client disconnects
             
             **Path Parameters:**
             - **session_id**: The unique identifier of the session
             
             **Request Body Examples:**
             
             *Using session's stored query:*
             ```json
             {
               "model_type": "groq",
               "model_name": "groq/openai/gpt-oss-120b",
               "model_key": "your-api-key-here"
             }
             ```
             
             *Override with new query:*
             ```json
             {
               "model_type": "groq", 
               "model_name": "groq/openai/gpt-oss-120b",
               "model_key": "your-api-key-here",
               "query": "8 Fe + S₈ → 8 FeS"
             }
             ```
             
             **Response:** Worker ID and confirmation message
             """)
async def start_session_worker(session_id: str, request: StartWorkerRequest):
    """Start a worker for a session"""
    try:
        model_name = request.model_name.strip()
        model_type = model_name.split("/")[0]
        assert model_type != "", "The model_name must contain the model_type as the prefix before the first '/'"

        model_key = get_model_key(model_type)

        assert model_key is not None, f"Unsupported model_type '{model_type}'. Supported types: 'groq', 'openai', 'azure'."

        worker_id = await session_runner.start_session_worker(
            session_id=session_id,
            model_type=model_type,
            model_name=model_name,
            model_key=model_key,
            query=request.query,  # Optional - will use session's query if None
            agent_name=request.agent_name
        )
        return {
            "message": "Worker started successfully",
            "session_id": session_id,
            "worker_id": worker_id
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error starting worker: {str(e)}")

@router.get("/sessions/{session_id}/worker-status",
            summary="Get worker status",
            description="""
            Check the current status of a worker running for a session.
            
            **Path Parameters:**
            - **session_id**: The unique identifier of the session
            
            **Response Statuses:**
            - `idle`: Worker created but not started
            - `running`: Worker actively processing reactions
            - `completed`: Worker finished successfully  
            - `error`: Worker encountered an error
            - `cancelled`: Worker was manually stopped
            
            **Example Response:**
            ```json
            {
              "session_id": "session-uuid",
              "worker_id": "worker-uuid",
              "status": "running", 
              "created_at": "2025-09-04T10:30:00"
            }
            ```
            """)
async def get_session_worker_status(session_id: str):
    """Get the status of a session's worker"""
    status = session_runner.get_worker_status(session_id)
    if status:
        return status
    else:
        raise HTTPException(status_code=404, detail="No active worker found for this session")

@router.post("/sessions/{session_id}/stop-worker",
             summary="Stop worker",
             description="""
             Stop a running worker for a session. This will gracefully cancel the background processing.
             
             **Path Parameters:**
             - **session_id**: The unique identifier of the session
             
             **Note:** Stopping a worker will preserve all progress made up to that point.
             """)
async def stop_session_worker(session_id: str):
    """Stop a worker for a session"""
    try:
        stopped = await session_runner.stop_worker(session_id)
        if stopped:
            return {
                "message": "Worker stopped successfully",
                "session_id": session_id
            }
        else:
            raise HTTPException(status_code=404, detail="No active worker found for this session")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error stopping worker: {str(e)}")

@router.get("/workers/status",
            summary="Get all workers status",
            description="""
            Get the status of all currently active workers across all sessions.
            
            **Response:** Overview of all active workers with their session IDs, worker IDs, and current status.
            
            **Example Response:**
            ```json
            {
              "active_workers": {
                "session-uuid-1": {
                  "worker_id": "worker-uuid-1",
                  "status": "running",
                  "created_at": "2025-09-04T10:30:00"
                }
              },
              "total_active": 1
            }
            ```
            """)
async def get_all_workers_status():
    """Get status of all active workers"""
    return {
        "active_workers": session_runner.get_all_workers_status(),
        "total_active": len(session_runner.active_workers)
    }

# Query Management Endpoints
@router.put("/sessions/{session_id}/query",
            summary="Update session query",
            description="""
            Update the chemical reaction query text for a session. This query will be used by future workers.
            
            **Path Parameters:**
            - **session_id**: The unique identifier of the session
            
            **Request Body:**
            ```json
            {
              "query": "8 Fe + S₈ → 8 FeS"
            }
            ```
            
            **Use Cases:**
            - Set initial query before running analysis
            - Update query for re-running with different reactions
            - Modify query based on previous results
            """)
async def update_session_query(session_id: str, request: UpdateQueryRequest, db: Session = Depends(get_db)):
    """Update the query for a session"""
    try:
        db_session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
        if not db_session:
            raise HTTPException(status_code=404, detail="Session not found")
        
        db_session.query = request.query
        db_session.last_updated = datetime.utcnow()
        db.commit()
        
        return {
            "message": "Query updated successfully",
            "session_id": session_id,
            "query": request.query
        }
    except Exception as e:
        db.rollback()
        raise HTTPException(status_code=500, detail=f"Error updating query: {str(e)}")

@router.get("/sessions/{session_id}/query",
            summary="Get session query",
            description="""
            Retrieve the current chemical reaction query text stored for a session.
            
            **Path Parameters:**
            - **session_id**: The unique identifier of the session
            
            **Response:**
            ```json
            {
              "session_id": "session-uuid",
              "query": "8 Fe + S₈ → 8 FeS"
            }
            ```
            """)
async def get_session_query(session_id: str, db: Session = Depends(get_db)):
    """Get the current query for a session"""
    try:
        db_session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
        if not db_session:
            raise HTTPException(status_code=404, detail="Session not found")
        
        return {
            "session_id": session_id,
            "query": db_session.query or ""
        }
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error retrieving query: {str(e)}")

@router.post("/sessions/{session_id}/rerun-with-query",
             status_code=status.HTTP_201_CREATED,
             summary="Rerun session with new query",
             description="""
             Update the session's query and restart the worker in a single operation. This is a convenience endpoint.
             
             **What it does:**
             1. Stops any currently running worker
             2. Updates the session's stored query
             3. Starts a new worker with the updated query
             
             **Path Parameters:**
             - **session_id**: The unique identifier of the session
             
             **Request Body:**
             ```json
             {
               "model_type": "groq",
               "model_name": "groq/openai/gpt-oss-120b",
               "model_key": "your-api-key-here",
               "query": "2 H₂O₂ → 2 H₂O + O₂",
               "agent_name": "master_openai_oss_agent"
             }
             ```
             
             **Use Cases:**
             - Quick iteration on different chemical reactions
             - Restarting analysis with corrected reaction text
             - Testing multiple reaction scenarios
             """)
async def rerun_session_with_new_query(session_id: str, request: StartWorkerRequest):
    """Update query and start/restart worker for a session"""
    try:
        # Stop existing worker if running
        await session_runner.stop_worker(session_id)
        
        # Start new worker with updated query
        worker_id = await session_runner.start_session_worker(
            session_id=session_id,
            model_type=request.model_type,
            model_name=request.model_name,
            model_key=request.model_key,
            query=request.query,  # This will update the session's query
            agent_name=request.agent_name
        )
        
        return {
            "message": "Session restarted with new query",
            "session_id": session_id,
            "worker_id": worker_id,
            "query": request.query
        }
    except ValueError as e:
        raise HTTPException(status_code=400, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error rerunning session: {str(e)}")