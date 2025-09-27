from pydantic import BaseModel
from datetime import datetime
from typing import Optional, List


# Session Pydantic Schemas
class SessionBase(BaseModel):
    name: str
    content: Optional[List[dict]] = []  # List of content dictionaries instead of JSON string
    state: Optional[str] = None
    worker_id: Optional[str] = None
    project_id: str
    file_paths: Optional[List[dict]] = []  # List of file path dictionaries instead of strings
    query: Optional[str] = ""  # Query text for the session


class SessionCreate(BaseModel):
    name: str
    content: Optional[List[dict]] = []  # List of content dictionaries
    state: Optional[str] = None
    project_id: str  # Required - must provide project_id when creating sessions
    file_paths: Optional[List[dict]] = []  # List of file path dictionaries
    query: Optional[str] = ""  # Query text for the session


class SessionRead(SessionBase):
    id: str
    last_updated: Optional[datetime]

    class Config:
        from_attributes = True


# Project Pydantic Schemas
class ProjectBase(BaseModel):
    id: str
    name: str
    description: Optional[str] = None
    sessions: List[str] = []


class ProjectCreate(BaseModel):
    name: str
    description: Optional[str] = None
    sessions: List[str] = []


class ProjectRead(ProjectBase):
    id: str

    class Config:
        from_attributes = True
