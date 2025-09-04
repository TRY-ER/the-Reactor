#!/usr/bin/env python3
"""
Database migration utility to add the query column to existing sessions table.
Run this script to update your existing database schema.
"""

import sqlite3
import os
from sqlalchemy import create_engine, inspect
from db import DATABASE_URL, SessionObj, Base, engine

def check_column_exists(table_name, column_name, db_path="projects.db"):
    """Check if a column exists in a table"""
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Get table info
        cursor.execute(f"PRAGMA table_info({table_name})")
        columns = cursor.fetchall()
        
        # Check if column exists
        column_names = [col[1] for col in columns]
        exists = column_name in column_names
        
        conn.close()
        return exists, column_names
    except Exception as e:
        print(f"Error checking column: {e}")
        return False, []

def add_query_column():
    """Add the query column to the sessions table if it doesn't exist"""
    
    # Extract database path from DATABASE_URL
    if DATABASE_URL.startswith("sqlite:///"):
        db_path = DATABASE_URL.replace("sqlite:///./", "")
    else:
        db_path = "projects.db"
    
    if not os.path.exists(db_path):
        print(f"Database file {db_path} not found. Creating new database...")
        Base.metadata.create_all(bind=engine)
        print("New database created with all required columns.")
        return
    
    # Check if query column exists
    exists, columns = check_column_exists("sessions", "query", db_path)
    
    if exists:
        print("Query column already exists in sessions table. No migration needed.")
        return
    
    print(f"Current columns in sessions table: {columns}")
    print("Query column not found. Adding query column...")
    
    try:
        # Connect to database and add column
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Add the query column with default empty string
        cursor.execute("ALTER TABLE sessions ADD COLUMN query TEXT DEFAULT ''")
        
        # Commit changes
        conn.commit()
        conn.close()
        
        print("Successfully added query column to sessions table.")
        
        # Verify the change
        exists_after, columns_after = check_column_exists("sessions", "query", db_path)
        if exists_after:
            print("Migration completed successfully!")
            print(f"Updated columns in sessions table: {columns_after}")
        else:
            print("Warning: Column addition may have failed.")
            
    except Exception as e:
        print(f"Error during migration: {e}")
        print("You may need to manually add the column or recreate the database.")

def migrate_worker_id_column():
    """Fix worker_id column to allow NULL values"""
    
    # Extract database path from DATABASE_URL
    if DATABASE_URL.startswith("sqlite:///"):
        db_path = DATABASE_URL.replace("sqlite:///./", "")
    else:
        db_path = "projects.db"
    
    if not os.path.exists(db_path):
        print(f"Database file {db_path} not found. Creating new database...")
        Base.metadata.create_all(bind=engine)
        print("New database created with all required columns.")
        return
    
    try:
        conn = sqlite3.connect(db_path)
        cursor = conn.cursor()
        
        # Check current worker_id column constraints
        cursor.execute("PRAGMA table_info(sessions)")
        columns = cursor.fetchall()
        
        worker_id_col = next((col for col in columns if col[1] == 'worker_id'), None)
        
        if worker_id_col and worker_id_col[3] == 1:  # notnull=1 means NOT NULL
            print("worker_id column is NOT NULL. Fixing to allow NULL values...")
            
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
                       CASE 
                           WHEN query IS NULL THEN '' 
                           ELSE query 
                       END as query
                FROM sessions
            """)
            
            # Drop old table and rename new one
            cursor.execute("DROP TABLE sessions")
            cursor.execute("ALTER TABLE sessions_new RENAME TO sessions")
            
            conn.commit()
            print("Successfully fixed worker_id column to allow NULL values!")
        else:
            print("worker_id column already allows NULL values. No migration needed.")
        
        conn.close()
        
    except Exception as e:
        print(f"Error during worker_id migration: {e}")
        print("You may need to manually fix the schema.")

def migrate_database():
    """Main migration function"""
    print("Starting database migration...")
    print("=" * 50)
    
    add_query_column()
    migrate_worker_id_column()
    
    print("=" * 50)
    print("Migration complete!")

if __name__ == "__main__":
    migrate_database()
