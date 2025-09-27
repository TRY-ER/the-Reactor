from sqlalchemy import create_engine, Column, Integer, String, DateTime
from sqlalchemy.types import Text
from datetime import datetime
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
import os

DATABASE_URL = os.getenv("DATABASE_URL", "sqlite:///./projects.db")

engine = create_engine(DATABASE_URL, connect_args={"check_same_thread": False})
SessionLocal = sessionmaker(autocommit=False, autoflush=False, bind=engine)
Base = declarative_base()


from sqlalchemy.types import TEXT
import json

class Project(Base):
    __tablename__ = "projects"
    id = Column(String, primary_key=True, index=True)
    name = Column(String, unique=True, index=True, nullable=False)
    description = Column(String, nullable=True)
    sessions = Column(TEXT, nullable=False, default="[]")  # Store session IDs as JSON string


# Session model

class SessionObj(Base):
    __tablename__ = "sessions"
    id = Column(String, primary_key=True, index=True)  # session id (string)
    name = Column(String, nullable=False)  # session name
    content = Column(Text, nullable=False, default="[]")  # JSON array as string to store list of content dictionaries
    state = Column(String, nullable=True)  # can be null initially
    last_updated = Column(DateTime, default=datetime.utcnow, onupdate=datetime.utcnow)
    worker_id = Column(String, nullable=True)  # can be null - workers are assigned later
    project_id = Column(String, nullable=False)  # project id (string)
    file_paths = Column(TEXT, nullable=False, default="[]")  # Store file paths as JSON string
    query = Column(TEXT, nullable=True, default="")  # Query text for the session - initially empty

# Function to safely initialize database with migrations
def init_database():
    """Initialize database with proper migration handling"""
    import sqlite3
    import os
    
    # Extract database path from DATABASE_URL
    if DATABASE_URL.startswith("sqlite:///"):
        db_path = DATABASE_URL.replace("sqlite:///./", "")
    else:
        db_path = "projects.db"
    
    # Create tables if database doesn't exist
    if not os.path.exists(db_path):
        print("Creating new database...")
        Base.metadata.create_all(bind=engine)
        print("Database created successfully!")
        return
    
    # Check if query column exists in existing database and fix worker_id nullable constraint
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Check if sessions table exists and has query column
        cursor.execute("PRAGMA table_info(sessions)")
        columns = cursor.fetchall()
        column_names = [col[1] for col in columns]
        
        # Add query column if missing
        if 'query' not in column_names:
            print("Adding missing query column to sessions table...")
            cursor.execute("ALTER TABLE sessions ADD COLUMN query TEXT DEFAULT ''")
            conn.commit()
            print("Query column added successfully!")
        
        # Check if worker_id is NOT NULL and fix it
        worker_id_col = next((col for col in columns if col[1] == 'worker_id'), None)
        if worker_id_col and worker_id_col[3] == 1:  # notnull=1 means NOT NULL
            print("Fixing worker_id column to allow NULL values...")
            
            # Create new table with correct schema
            cursor.execute("""
                CREATE TABLE sessions_new (
                    id VARCHAR PRIMARY KEY,
                    name VARCHAR NOT NULL,
                    content TEXT NOT NULL DEFAULT '[]',
                    state VARCHAR,
                    last_updated DATETIME,
                    worker_id VARCHAR,
                    project_id VARCHAR NOT NULL,
                    file_paths TEXT NOT NULL DEFAULT '[]',
                    query TEXT DEFAULT ''
                )
            """)
            
            # Copy data from old table to new table
            cursor.execute("""
                INSERT INTO sessions_new (id, name, content, state, last_updated, worker_id, project_id, file_paths, query)
                SELECT id, name, content, state, last_updated, worker_id, project_id, file_paths, 
                       CASE WHEN query IS NULL THEN '' ELSE query END
                FROM sessions
            """)
            
            # Drop old table and rename new one
            cursor.execute("DROP TABLE sessions")
            cursor.execute("ALTER TABLE sessions_new RENAME TO sessions")
            
            conn.commit()
            print("worker_id column fixed to allow NULL values!")
        
        conn.close()
        
    except Exception as e:
        print(f"Warning: Could not verify/migrate columns: {e}")
        print("You may need to run the migration script manually.")
    
    # Ensure all tables are created/updated
    Base.metadata.create_all(bind=engine)

# Initialize the database safely
init_database()
