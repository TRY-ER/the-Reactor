from fastapi import APIRouter, HTTPException, Depends, status
from fastapi.responses import StreamingResponse
from sqlalchemy.orm import Session
from datetime import datetime
import asyncio
import json
import uuid
import os
import csv
from typing import Optional, List, Dict, Any

from db import Project, SessionLocal, SessionObj
from models import SessionBase, SessionCreate, SessionRead, ProjectBase, ProjectCreate, ProjectRead
from file_system_handler import file_system_handler
from session_runner import SessionRunner 
from viewer import Viewer
import json
import uuid
from utils import get_model_key

session_runner = SessionRunner("./file_system")

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
             
             **Note:** 
             - project_id is required. You must create a project first using the /projects/ endpoint, then use its ID when creating sessions.
             - New sessions start in "init" state. They change to "active" when a worker is started.
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
        state=session.state or "init",  # Set initial state to "init" if not provided
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
        
        # Add session ID to the project's sessions list
        sessions_list = json.loads(existing_project.sessions)
        if str(db_session.id) not in sessions_list:
            sessions_list.append(str(db_session.id))
            existing_project.sessions = json.dumps(sessions_list)
            db.commit()
        
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
    
    # Get the project_id before deleting the session
    project_id = db_session.project_id
    
    # Delete the session from sessions table
    db.delete(db_session)
    
    # Remove session_id from the project's sessions list
    db_project = db.query(Project).filter(Project.id == project_id).first()
    if db_project:
        # Parse the sessions JSON, remove the session_id, and update
        sessions_list = json.loads(db_project.sessions)
        if session_id in sessions_list:
            sessions_list.remove(session_id)
            db_project.sessions = json.dumps(sessions_list)
    
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
    
    # Delete all sessions associated with this project
    db.query(SessionObj).filter(SessionObj.project_id == project_id).delete()
    
    # Delete the project
    db.delete(db_project)
    db.commit()
    
    # Delete the project directory
    file_system_handler.delete_project_directory(project_id)
    return {"ok": True}

@router.get("/projects/{project_id}/master-csv",
            summary="Get project master CSV file",
            description="""
            Retrieve the content of a project's master CSV file as structured data.
            
            **Path Parameters:**
            - **project_id**: The unique identifier of the project
            
            **Query Parameters:**
            - **file_name**: Optional name of the CSV file (default: "master.csv")
            
            **Response:**
            ```json
            {
                "project_id": "project-uuid",
                "file_path": "/path/to/master.csv",
                "data": [
                    {"index": "0", "reactions_text": "Fe + S -> FeS", ...},
                    {"index": "1", "reactions_text": "H2 + O2 -> H2O", ...}
                ],
                "headers": ["index", "reactions_text", "reactions_composition_SMILES", ...],
                "exists": true,
                "row_count": 2
            }
            ```
            
            **When file doesn't exist:**
            ```json
            {
                "project_id": "project-uuid",
                "file_path": "/path/to/master.csv",
                "data": null,
                "headers": null,
                "exists": false,
                "row_count": 0,
                "message": "Master CSV file not found. You need to create the file first or merge session data.",
                "detail": "The master.csv file will be created when sessions are merged or data is saved to the project."
            }
            ```
            
            **Error Cases:**
            - 404: Project not found
            - 500: File system or CSV parsing error
            """)
def get_project_master_csv(project_id: str, file_name: str = "master.csv", db: Session = Depends(get_db)):
    """Get the master CSV file content for a project as structured data."""
    try:
        # Verify project exists
        db_project = db.query(Project).filter(Project.id == project_id).first()
        if not db_project:
            raise HTTPException(status_code=404, detail="Project not found")
        
        # Get project directory path
        project_dir = file_system_handler.get_project_directory(project_id)
        
        # Check if project directory exists
        if not os.path.exists(project_dir):
            raise HTTPException(status_code=404, detail=f"Project directory for project {project_id} does not exist")
        
        # Create the full file path
        file_path = os.path.join(project_dir, file_name)
        
        # Check if file exists and read content
        if os.path.exists(file_path):
            try:
                data = []
                headers = []
                
                with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
                    reader = csv.DictReader(csvfile)
                    headers = reader.fieldnames if reader.fieldnames else []
                    
                    for row in reader:
                        data.append(dict(row))
                
                return {
                    "project_id": project_id,
                    "file_path": file_path,
                    "data": data,
                    "headers": headers,
                    "exists": True,
                    "row_count": len(data)
                }
            except Exception as e:
                raise HTTPException(status_code=500, detail=f"Failed to read CSV file {file_name}: {str(e)}")
        else:
            return {
                "project_id": project_id,
                "file_path": file_path,
                "data": None,
                "headers": None,
                "exists": False,
                "row_count": 0,
                "message": "Master CSV file not found. You need to create the file first or merge session data.",
                "detail": f"The {file_name} file will be created when sessions are merged or data is saved to the project."
            }
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error retrieving project CSV file: {str(e)}")

@router.post("/projects/{project_id}/save-master-csv",
             status_code=status.HTTP_201_CREATED,
             summary="Save data to project master CSV",
             description="""
             Save structured CSV data directly to a project's master CSV file. This endpoint allows
             direct manipulation of the project's master dataset.
             
             **Path Parameters:**
             - **project_id**: The unique identifier of the project
             
             **Request Body:**
             ```json
             {
               "data": [
                 {
                   "index": "0",
                   "reactions_text": "8 Fe + S₈ → 8 FeS",
                   "reactions_composition_SMILES": "C=C",
                   "invalid_SMILE": false,
                   "reactions_SMILES": "[Fe+2].[S-2]",
                   "invalid_reaction_SMILES": false,
                   "reactions_SMARTS": "[Fe:1]>>[Fe:1][S]",
                   "invalid_reaction_SMARTS": false,
                   "reactions_SMIRKS": "[Fe:1]>>[Fe:1][S:2]",
                   "invalid_reaction_SMIRKS": false
                 }
               ],
               "file_name": "master.csv",
               "mode": "overwrite"
             }
             ```
             
             **Parameters:**
             - **data**: List of dictionaries containing reaction data
             - **file_name**: Name of the CSV file (default: "master.csv")
             - **mode**: Write mode - "overwrite" (default) or "append"
             
             **Response:**
             ```json
             {
               "message": "Master CSV data saved successfully",
               "project_id": "project-uuid",
               "file_path": "/path/to/master.csv",
               "rows_saved": 1,
               "file_name": "master.csv",
               "mode": "overwrite"
             }
             ```
             
             **Use Cases:**
             - Initialize project master CSV with seed data
             - Bulk update project datasets
             - Import external data into project
             - Reset or overwrite corrupted data
             - Append additional curated data
             """)
def save_project_master_csv(project_id: str, request: dict, db: Session = Depends(get_db)):
    """Save CSV data directly to a project's master CSV file."""
    try:
        # Verify project exists
        db_project = db.query(Project).filter(Project.id == project_id).first()
        if not db_project:
            raise HTTPException(status_code=404, detail="Project not found")
        
        # Extract request parameters
        data = request.get("data", [])
        file_name = request.get("file_name", "master.csv")
        mode = request.get("mode", "overwrite")  # "overwrite" or "append"
        
        # Validate data
        if not data:
            raise HTTPException(status_code=400, detail="Data array cannot be empty")
        
        if not isinstance(data, list):
            raise HTTPException(status_code=400, detail="Data must be a list of dictionaries")
        
        if not isinstance(data[0], dict):
            raise HTTPException(status_code=400, detail="Data must contain dictionaries")
        
        if mode not in ["overwrite", "append"]:
            raise HTTPException(status_code=400, detail="Mode must be 'overwrite' or 'append'")
        
        # Get project directory
        project_dir = file_system_handler.get_project_directory(project_id)
        
        # Ensure project directory exists
        os.makedirs(project_dir, exist_ok=True)
        
        # Create the full file path
        file_path = os.path.join(project_dir, file_name)
        
        # Define standard headers (same as used in session CSV)
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
        
        # For append mode, read existing data first
        existing_data = []
        if mode == "append" and os.path.exists(file_path):
            try:
                with open(file_path, 'r', newline='', encoding='utf-8') as csvfile:
                    reader = csv.DictReader(csvfile)
                    for row in reader:
                        existing_data.append(dict(row))
            except Exception as e:
                raise HTTPException(status_code=500, detail=f"Failed to read existing CSV file: {str(e)}")
        
        # Combine data for append mode
        if mode == "append":
            all_data = existing_data + data
        else:
            all_data = data
        
        # Write CSV file
        try:
            with open(file_path, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                
                # Write headers
                writer.writerow(HEADERS)
                
                # Write data rows
                for row_dict in all_data:
                    row = []
                    for header in HEADERS:
                        value = row_dict.get(header, '')
                        row.append(value)
                    writer.writerow(row)
            
            return {
                "message": "Master CSV data saved successfully",
                "project_id": project_id,
                "file_path": file_path,
                "rows_saved": len(data),
                "total_rows": len(all_data),
                "file_name": file_name,
                "mode": mode
            }
            
        except Exception as e:
            raise HTTPException(status_code=500, detail=f"Failed to write CSV file: {str(e)}")
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error saving project CSV data: {str(e)}")

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

@router.post("/sessions/{session_id}/save-csv",
             status_code=status.HTTP_201_CREATED,
             summary="Save CSV data to session",
             description="""
             Save structured CSV data to a session's file system. This endpoint allows frontend clients
             to save processed reaction data directly to the session's CSV file.
             
             **Path Parameters:**
             - **session_id**: The unique identifier of the session
             
             **Request Body:**
             ```json
             {
               "data": [
                 {
                   "index": "0",
                   "reactions_text": "8 Fe + S₈ → 8 FeS",
                   "reactions_composition_SMILES": "C=C",
                   "invalid_SMILE": false,
                   "reactions_SMILES": "[Fe+2].[S-2]",
                   "invalid_reaction_SMILES": false,
                   "reactions_SMARTS": "[Fe:1]>>[Fe:1][S]",
                   "invalid_reaction_SMARTS": false,
                   "reactions_SMIRKS": "[Fe:1]>>[Fe:1][S:2]",
                   "invalid_reaction_SMIRKS": false
                 }
               ],
               "file_name": "custom_data.csv"
             }
             ```
             
             **Expected CSV Headers:**
             - index, reaction_type, reactions_text, reactions_composition_SMILES
             - invalid_SMILE, reactions_SMILES, invalid_reaction_SMILES
             - reactions_SMARTS, invalid_reaction_SMARTS, reactions_SMIRKS, invalid_reaction_SMIRKS
             
             **Response:**
             ```json
             {
               "message": "CSV data saved successfully",
               "session_id": "session-uuid",
               "project_id": "project-uuid",
               "file_path": "/path/to/data.csv",
               "rows_saved": 1,
               "file_name": "data.csv"
             }
             ```
             
             **Use Cases:**
             - Save manually processed reaction data from frontend
             - Store results from frontend-based analysis tools
             - Create custom CSV files with specific reaction datasets
             - Backup or export processed data
             """)
def save_csv_to_session(session_id: str, request: dict, db: Session = Depends(get_db)):
    """Save CSV data to a session's file system."""
    try:
        # Get session to retrieve project_id
        db_session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
        if not db_session:
            raise HTTPException(status_code=404, detail="Session not found")
        
        project_id = str(db_session.project_id)
        
        # Validate that data is not empty
        if not request.get("data"):
            raise HTTPException(status_code=400, detail="Data array cannot be empty")
        
        # Get file_name with default
        file_name = request.get("file_name", "data.csv")
        
        # Save the CSV data using file_system_handler
        file_path = file_system_handler.save_data_to_session_csv(
            project_id=project_id,
            session_id=session_id,
            data=request["data"],
            file_name=file_name
        )
        
        return {
            "message": "CSV data saved successfully",
            "session_id": session_id,
            "project_id": project_id,
            "file_path": file_path,
            "rows_saved": len(request["data"]),
            "file_name": file_name
        }
        
    except HTTPException:
        raise
    except TypeError as e:
        raise HTTPException(status_code=400, detail=f"Invalid data format: {str(e)}")
    except ValueError as e:
        raise HTTPException(status_code=400, detail=f"Validation error: {str(e)}")
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error saving CSV data: {str(e)}")

@router.post("/sessions/{session_id}/merge-csv",
             status_code=status.HTTP_200_OK,
             summary="Merge session CSV with project CSV",
             description="""
             Merge a session's CSV data with the project's master CSV file. This endpoint consolidates
             individual session results into the project's cumulative dataset.
             
             **Path Parameters:**
             - **session_id**: The unique identifier of the session
             
             **Query Parameters:**
             - **session_file_name**: Optional name of the session CSV file (default: "data.csv")
             - **project_file_name**: Optional name of the project CSV file (default: "master.csv")
             
             **Process:**
             1. Reads the session's CSV data
             2. Reads the existing project master CSV
             3. Merges data while preventing duplicates
             4. Combines headers from both files
             5. Writes the merged data back to the project CSV
             
             **Response:**
             ```json
             {
               "message": "CSV merge completed successfully",
               "session_id": "session-uuid",
               "project_id": "project-uuid",
               "merge_results": {
                 "success": true,
                 "session_file_path": "/path/to/session/data.csv",
                 "project_file_path": "/path/to/project/master.csv",
                 "rows_read": 10,
                 "rows_added": 8,
                 "rows_skipped": 2,
                 "total_rows_after_merge": 25,
                 "message": "Successfully merged 8 new rows from session CSV to project CSV"
               }
             }
             ```
             
             **Use Cases:**
             - Consolidate session analysis results into project dataset
             - Build cumulative databases of chemical reaction analysis
             - Combine multiple session results for comprehensive analysis
             - Create master datasets for research and reporting
             
             **Error Cases:**
             - 404: Session not found
             - 404: Session CSV file doesn't exist
             - 400: Invalid file names or data format
             - 500: File system or merge operation errors
             """)
def merge_session_csv_to_project(
    session_id: str, 
    session_file_name: str = "data.csv",
    project_file_name: str = "master.csv", 
    db: Session = Depends(get_db)
):
    """Merge session CSV data with project master CSV."""
    try:
        # Get session to retrieve project_id
        db_session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
        if not db_session:
            raise HTTPException(status_code=404, detail="Session not found")
        
        project_id = str(db_session.project_id)
        
        # Perform the merge operation using file_system_handler
        merge_results = file_system_handler.merge_session_csv_to_project(
            project_id=project_id,
            session_id=session_id,
            session_file_name=session_file_name,
            project_file_name=project_file_name
        )
        
        return {
            "message": "CSV merge completed successfully",
            "session_id": session_id,
            "project_id": project_id,
            "merge_results": merge_results
        }
        
    except HTTPException:
        raise
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error during CSV merge: {str(e)}")

@router.post("/projects/{project_id}/merge-all-sessions",
             status_code=status.HTTP_200_OK,
             summary="Merge all session CSVs with project CSV",
             description="""
             Merge CSV data from all sessions in a project with the project's master CSV file.
             This is a batch operation that consolidates all session results at once.
             
             **Path Parameters:**
             - **project_id**: The unique identifier of the project
             
             **Query Parameters:**
             - **session_file_name**: Optional name of the session CSV files (default: "data.csv")
             - **project_file_name**: Optional name of the project CSV file (default: "master.csv")
             
             **Response:**
             ```json
             {
               "message": "Batch CSV merge completed",
               "project_id": "project-uuid",
               "sessions_processed": 3,
               "sessions_successful": 2,
               "sessions_failed": 1,
               "total_rows_added": 25,
               "merge_details": [
                 {
                   "session_id": "session-1",
                   "success": true,
                   "rows_added": 10,
                   "message": "Successfully merged"
                 },
                 {
                   "session_id": "session-2", 
                   "success": false,
                   "error": "CSV file not found"
                 }
               ]
             }
             ```
             
             **Use Cases:**
             - Consolidate all session results into master dataset
             - Bulk data processing for project completion
             - Generate comprehensive project reports
             - Prepare data for analysis or export
             """)
def merge_all_sessions_to_project(
    project_id: str,
    session_file_name: str = "data.csv",
    project_file_name: str = "master.csv",
    db: Session = Depends(get_db)
):
    """Merge CSV data from all sessions in a project with the project master CSV."""
    try:
        # Verify project exists and get all its sessions
        db_project = db.query(Project).filter(Project.id == project_id).first()
        if not db_project:
            raise HTTPException(status_code=404, detail="Project not found")
        
        # Parse sessions list from project
        sessions_list = json.loads(str(db_project.sessions)) if str(db_project.sessions) else []
        
        if not sessions_list:
            return {
                "message": "No sessions found in project",
                "project_id": project_id,
                "sessions_processed": 0,
                "sessions_successful": 0,
                "sessions_failed": 0,
                "total_rows_added": 0,
                "merge_details": []
            }
        
        merge_details = []
        sessions_successful = 0
        sessions_failed = 0
        total_rows_added = 0
        
        # Process each session
        for session_id in sessions_list:
            try:
                merge_result = file_system_handler.merge_session_csv_to_project(
                    project_id=project_id,
                    session_id=session_id,
                    session_file_name=session_file_name,
                    project_file_name=project_file_name
                )
                
                if merge_result['success']:
                    sessions_successful += 1
                    total_rows_added += merge_result['rows_added']
                    merge_details.append({
                        "session_id": session_id,
                        "success": True,
                        "rows_added": merge_result['rows_added'],
                        "rows_skipped": merge_result['rows_skipped'],
                        "message": merge_result['message']
                    })
                else:
                    sessions_failed += 1
                    merge_details.append({
                        "session_id": session_id,
                        "success": False,
                        "error": "Merge operation failed"
                    })
                    
            except Exception as e:
                sessions_failed += 1
                merge_details.append({
                    "session_id": session_id,
                    "success": False,
                    "error": str(e)
                })
        
        return {
            "message": "Batch CSV merge completed",
            "project_id": project_id,
            "sessions_processed": len(sessions_list),
            "sessions_successful": sessions_successful,
            "sessions_failed": sessions_failed,
            "total_rows_added": total_rows_added,
            "merge_details": merge_details
        }
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error during batch merge: {str(e)}")

@router.get("/sessions/{session_id}/prompt-file",
            summary="Get session prompt text file",
            description="""
            Retrieve the prompt text file content for a session.
            
            **Path Parameters:**
            - **session_id**: The unique identifier of the session
            
            **Query Parameters:**
            - **file_name**: Optional name of the text file (default: "master_prompt.txt")
            
            **Response:**
            ```json
            {
                "session_id": "session-uuid",
                "project_id": "project-uuid", 
                "file_path": "/path/to/file.txt",
                "content": "text content of the file",
                "exists": true
            }
            ```
            
            **When file doesn't exist:**
            ```json
            {
                "session_id": "session-uuid",
                "project_id": "project-uuid",
                "file_path": "/path/to/file.txt", 
                "content": null,
                "exists": false,
                "message": "Prompt text file not found. You need to complete a full session run with a worker to generate this file.",
                "detail": "The file 'master_prompt.txt' will be created when a session worker processes the query and generates results."
            }
            ```
            
            **Error Cases:**
            - 404: Session not found
            - 404: Session directory or file doesn't exist
            - 500: File system error
            """)
def get_session_prompt_file(session_id: str, file_name: str = "master_prompt.txt", db: Session = Depends(get_db)):
    """Get the prompt text file content for a session."""
    try:
        # Get session to retrieve project_id
        db_session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
        if not db_session:
            raise HTTPException(status_code=404, detail="Session not found")
        
        project_id = str(db_session.project_id)
        
        # Get the prompt text file content using file_system_handler
        result = file_system_handler.get_session_prompt_text_file(project_id, session_id, file_name)
        
        # If file doesn't exist, provide helpful message
        if not result["exists"]:
            return {
                "session_id": session_id,
                "project_id": project_id,
                "file_path": result["file_path"],
                "content": None,
                "exists": False,
                "message": "Prompt text file not found. You need to complete a full session run with a worker to generate this file.",
                "detail": f"The file '{file_name}' will be created when a session worker processes the query and generates results."
            }
        
        return {
            "session_id": session_id,
            "project_id": project_id,
            "file_path": result["file_path"],
            "content": result["content"],
            "exists": result["exists"]
        }
        
    except HTTPException:
        raise
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error retrieving prompt file: {str(e)}")

@router.get("/sessions/{session_id}/csv-file",
            summary="Get session CSV file",
            description="""
            Retrieve the CSV file content for a session as structured data.
            
            **Path Parameters:**
            - **session_id**: The unique identifier of the session
            
            **Query Parameters:**
            - **file_name**: Optional name of the CSV file (default: "data.csv")
            
            **Response:**
            ```json
            {
                "session_id": "session-uuid",
                "project_id": "project-uuid",
                "file_path": "/path/to/file.csv",
                "data": [
                    {"index": "0", "reactions_text": "Fe + S -> FeS", ...},
                    {"index": "1", "reactions_text": "H2 + O2 -> H2O", ...}
                ],
                "headers": ["index", "reactions_text", "reactions_composition_SMILES", ...],
                "exists": true,
                "row_count": 2
            }
            ```
            
            **When file doesn't exist:**
            ```json
            {
                "session_id": "session-uuid",
                "project_id": "project-uuid",
                "file_path": "/path/to/file.csv",
                "data": null,
                "headers": null,
                "exists": false,
                "row_count": 0,
                "message": "CSV data file not found. You need to complete a full session run with a worker to generate this file.",
                "detail": "The file 'data.csv' will be created when a session worker processes chemical reactions and generates SMILES, SMARTS, and SMIRKS data."
            }
            ```
            
            **Error Cases:**
            - 404: Session not found
            - 404: Session directory or file doesn't exist
            - 500: File system or CSV parsing error
            """)
def get_session_csv_file(session_id: str, file_name: str = "data.csv", db: Session = Depends(get_db)):
    """Get the CSV file content for a session as structured data."""
    try:
        # Get session to retrieve project_id
        db_session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
        if not db_session:
            raise HTTPException(status_code=404, detail="Session not found")
        
        project_id = str(db_session.project_id)
        
        # Get the CSV file content using file_system_handler
        result = file_system_handler.get_session_csv_file(project_id, session_id, file_name)
        
        # If file doesn't exist, provide helpful message
        if not result["exists"]:
            return {
                "session_id": session_id,
                "project_id": project_id,
                "file_path": result["file_path"],
                "data": None,
                "headers": None,
                "exists": False,
                "row_count": 0,
                "message": "CSV data file not found. You need to complete a full session run with a worker to generate this file.",
                "detail": f"The file '{file_name}' will be created when a session worker processes chemical reactions and generates SMILES, SMARTS, and SMIRKS data."
            }
        
        return {
            "session_id": session_id,
            "project_id": project_id,
            "file_path": result["file_path"],
            "data": result["data"],
            "headers": result["headers"],
            "exists": result["exists"],
            "row_count": len(result["data"]) if result["data"] else 0
        }
        
    except HTTPException:
        raise
    except ValueError as e:
        raise HTTPException(status_code=404, detail=str(e))
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Error retrieving CSV file: {str(e)}")

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


# SSE Content Streaming Endpoint
@router.get("/sessions/{session_id}/stream",
            summary="Stream session content via SSE",
            description="""
            Stream session content in real-time using Server-Sent Events (SSE).
            
            This endpoint allows clients to receive real-time updates as a session worker
            processes chemical reactions and generates results.
            
            **Path Parameters:**
            - **session_id**: The unique identifier of the session to stream
            
            **Query Parameters:**
            - **client_id**: Optional unique identifier for the client connection
            
            **Response:** Server-Sent Events stream with content updates.
            
            **Event Format:**
            ```
            data: {"type": "info", "message": "Processing...", "timestamp": "..."}
            
            data: {"type": "data", "content": {...}, "timestamp": "..."}
            
            data: {"_type": "completion", "status": "completed"}
            ```
            
            **Usage:**
            ```javascript
            const eventSource = new EventSource('/sessions/{session_id}/stream?client_id=client123');
            eventSource.onmessage = function(event) {
                const data = JSON.parse(event.data);
                console.log('Received:', data);
            };
            ```
            """,
            tags=["Session Streaming"])
async def stream_session_content(
    session_id: str, 
    client_id: Optional[str] = None
):
    """Stream session content via Server-Sent Events"""
    
    async def generate_sse():
        """Generate SSE formatted data"""
        try:
            async for content in session_runner.stream_session_content(session_id, client_id):
                # Format as SSE
                data = json.dumps(content)
                yield f"data: {data}\n\n"
                
                # If this is a completion event, stop streaming
                if content.get("_type") == "completion":
                    break
                    
        except Exception as e:
            # Send error as final event
            error_data = {
                "_type": "completion",
                "status": "error", 
                "error": f"Stream error: {str(e)}"
            }
            yield f"data: {json.dumps(error_data)}\n\n"
    
    return StreamingResponse(
        generate_sse(),
        media_type="text/event-stream",
        headers={
            "Cache-Control": "no-cache",
            "Connection": "keep-alive",
            "X-Accel-Buffering": "no"  # Disable nginx buffering
        }
    )


@router.get("/sessions/{session_id}/streaming-status",
            summary="Get session streaming status",
            description="""
            Get the streaming status for a session including active streaming clients
            and current worker state.
            
            **Path Parameters:**
            - **session_id**: The unique identifier of the session
            
            **Response:** Streaming status information including worker state and active clients.
            """,
            tags=["Session Streaming"])
async def get_session_streaming_status(session_id: str):
    """Get streaming status for a session"""
    status = session_runner.get_streaming_status(session_id)
    return {"session_id": session_id, **status}


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

class RenderRequest(BaseModel):
    """Request model for rendering chemical structures
    
    Example:
    {
        "input_type": "SMILES",
        "data": "CCO"
    }
    or
    {
        "input_type": "SMARTS", 
        "data": "[C:1]=[O:2]>>[C:1][OH:2]"
    }
    """
    input_type: str = Field(..., description="Type of input data (SMILES or SMARTS)")
    data: str = Field(..., description="Chemical structure data to render")

# Chemical Structure Rendering Endpoint
@router.post("/render",
            summary="Render chemical structures",
            description="""
            Render chemical structures (SMILES or SMARTS) as base64-encoded PNG images.
            
            **Supported Input Types:**
            - **SMILES**: Simplified Molecular Input Line Entry System
            - **SMARTS**: SMiles ARbitrary Target Specification (for reactions)
            
            **Request Body:**
            ```json
            {
                "input_type": "SMILES",
                "data": "CCO"
            }
            ```
            
            **Response:**
            ```json
            {
                "input_type": "SMILES",
                "data": "CCO", 
                "image_base64": "iVBORw0KGgoAAAANSUhEUgAA...",
                "success": true
            }
            ```
            
            **Example Usage:**
            - SMILES: `{"input_type": "SMILES", "data": "CCO"}` (ethanol)
            - SMARTS: `{"input_type": "SMARTS", "data": "[C:1]=[O:2]>>[C:1][OH:2]"}` (carbonyl to alcohol)
            """,
            tags=["Chemical Rendering"])
async def render_chemical_structure(request: RenderRequest):
    """Render chemical structures as base64 images"""
    try:
        # Validate input type
        if request.input_type not in ["SMILES", "SMARTS"]:
            raise HTTPException(
                status_code=400, 
                detail=f"Invalid input_type '{request.input_type}'. Must be 'SMILES' or 'SMARTS'."
            )
        
        # Create viewer and render
        viewer = Viewer(request.input_type)
        image_base64 = viewer.render(request.data)
        
        # Check if rendering failed
        if isinstance(image_base64, str) and ("Invalid" in image_base64 or "Error" in image_base64):
            raise HTTPException(
                status_code=400,
                detail=f"Rendering failed: {image_base64}"
            )
        
        return {
            "input_type": request.input_type,
            "data": request.data,
            "image_base64": image_base64,
            "success": True
        }
        
    except HTTPException:
        raise
    except Exception as e:
        raise HTTPException(
            status_code=500,
            detail=f"Error rendering chemical structure: {str(e)}"
        )

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
               "model_name": "groq/openai/gpt-oss-120b",
             }
             ```
             
             *Override with new query:*
             ```json
             {
               "model_name": "groq/openai/gpt-oss-120b",
               "query": "8 Fe + S₈ → 8 FeS"
             }
             ```
             
             **Response:** Worker ID and confirmation message
             """)
async def start_session_worker(session_id: str, request: StartWorkerRequest, db: Session = Depends(get_db)):
    """Start a worker for a session"""
    try:
        model_name = request.model_name.strip()
        model_type = model_name.split("/")[0]
        assert model_type != "", "The model_name must contain the model_type as the prefix before the first '/'"

        model_key = get_model_key(model_type)

        assert model_key is not None, f"Unsupported model_type '{model_type}'. Supported types: 'groq', 'openai', 'azure'."

        # Update session state to "active" when starting worker
        db_session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
        if db_session:
            db_session.state = "active"
            db_session.last_updated = datetime.utcnow()
            db.commit()

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

async def rerun_session_with_new_query(session_id: str, request: StartWorkerRequest, db: Session = Depends(get_db)):
    """Update query and start/restart worker for a session"""
    try:
        # Stop existing worker if running
        await session_runner.stop_worker(session_id)
        
        # Update session state to "active" when restarting worker
        db_session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
        if db_session:
            db_session.state = "active"
            db_session.last_updated = datetime.utcnow()
            db.commit()
        
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