import os
import shutil
import csv
import json
from sqlalchemy.orm import Session
from db import SessionLocal, SessionObj
from datetime import datetime

class FileSystemHandler:
    def __init__(self, base_directory: str = "./file_system"):
        self.base_directory = base_directory
        # Ensure base directory exists
        os.makedirs(self.base_directory, exist_ok=True)

    def create_project_directory(self, project_id: str):
        """Create a directory for the project."""
        project_dir = os.path.join(self.base_directory, f"project_{project_id}")
        os.makedirs(project_dir, exist_ok=True)
        # Create master.csv with default headers
        return project_dir

    def create_master_csv(self, project_id: str, headers: list):
        """Create master.csv with specified headers."""
        project_dir = self.get_project_directory(project_id)
        master_csv_path = os.path.join(project_dir, "master.csv")
        with open(master_csv_path, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            writer.writerow(headers)

    def create_session_directory(self, project_id: str, session_id: str):
        """Create a directory for the session under the project."""
        project_dir = self.get_project_directory(project_id)
        sessions_dir = os.path.join(project_dir, "sessions")
        os.makedirs(sessions_dir, exist_ok=True)
        session_dir = os.path.join(sessions_dir, f"session_{session_id}")
        os.makedirs(session_dir, exist_ok=True)
        return session_dir

    def delete_project_directory(self, project_id: str):
        """Delete the project directory."""
        project_dir = os.path.join(self.base_directory, f"project_{project_id}")
        if os.path.exists(project_dir):
            shutil.rmtree(project_dir)

    def get_project_directory(self, project_id: str):
        """Get the path to the project directory."""
        return os.path.join(self.base_directory, f"project_{project_id}")

    def add_content_to_session(self, session_id: str, content_dict: dict):
        """Add content dictionary to the session's content list."""
        db: Session = SessionLocal()
        try:
            # Get the session
            session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
            if not session:
                raise ValueError(f"Session with id {session_id} not found")
            
            # Parse current content list
            current_content = json.loads(session.content) if session.content else []
            
            # Add new content dictionary
            current_content.append(content_dict)
            
            # Update the session
            session.content = json.dumps(current_content)
            session.last_updated = datetime.utcnow()
            
            db.commit()
            return True
        except Exception as e:
            db.rollback()
            raise e
        finally:
            db.close()

    def add_file_path_to_session(self, session_id: str, file_path_dict: dict):
        """Add file path dictionary to the session's file_paths list."""
        db: Session = SessionLocal()
        try:
            # Get the session
            session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
            if not session:
                raise ValueError(f"Session with id {session_id} not found")
            
            # Parse current file paths list
            current_file_paths = json.loads(session.file_paths) if session.file_paths else []
            
            # Add new file path dictionary
            current_file_paths.append(file_path_dict)
            
            # Update the session
            session.file_paths = json.dumps(current_file_paths)
            session.last_updated = datetime.utcnow()
            
            db.commit()
            return True
        except Exception as e:
            db.rollback()
            raise e
        finally:
            db.close()

    def get_session_content(self, session_id: str):
        """Get the content list for a session."""
        db: Session = SessionLocal()
        try:
            session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
            if not session:
                raise ValueError(f"Session with id {session_id} not found")
            
            return json.loads(session.content) if session.content else []
        finally:
            db.close()

    def get_session_file_paths(self, session_id: str):
        """Get the file paths list for a session."""
        db: Session = SessionLocal()
        try:
            session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
            if not session:
                raise ValueError(f"Session with id {session_id} not found")
            
            return json.loads(session.file_paths) if session.file_paths else []
        finally:
            db.close()

# Singleton instance
file_system_handler = FileSystemHandler()
