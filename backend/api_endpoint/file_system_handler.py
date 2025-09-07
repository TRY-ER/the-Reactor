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

    def clear_session_content(self, session_id: str):
        """Clear all content for a session by setting it to an empty list."""
        db: Session = SessionLocal()
        try:
            session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
            if not session:
                raise ValueError(f"Session with id {session_id} not found")
            
            # Set content to empty list
            session.content = json.dumps([])
            session.last_updated = datetime.utcnow()
            
            db.commit()
            return True
        except Exception as e:
            db.rollback()
            raise e
        finally:
            db.close()

    def add_session_prompt_text_file(self, project_id: str, session_id: str, file_data: dict):
        """
        Add a text file to the session directory.
        
        Args:
            project_id (str): The project ID
            session_id (str): The session ID
            file_data (dict): Dictionary containing 'file_name' and 'text_obj' keys
        
        Returns:
            str: The full path to the created file
        
        Raises:
            ValueError: If file_data doesn't contain required keys
            Exception: If file writing fails
        """
        # Validate input
        if not isinstance(file_data, dict):
            raise ValueError("file_data must be a dictionary")
        
        if 'file_name' not in file_data or 'text_obj' not in file_data:
            raise ValueError("file_data must contain 'file_name' and 'text_obj' keys")
        
        file_name = file_data['file_name']
        text_obj = file_data['text_obj']
        
        # Get the session directory (create if it doesn't exist)
        session_dir = self.create_session_directory(project_id, session_id)
        
        # Create the full file path
        file_path = os.path.join(session_dir, file_name)
        
        try:
            # Write the text content to the file
            with open(file_path, 'w', encoding='utf-8') as f:
                f.write(str(text_obj))
            
            return file_path
        except Exception as e:
            raise Exception(f"Failed to write file {file_name}: {str(e)}")

    def save_data_to_session_csv(self, project_id: str, session_id: str, data: list, file_name: str = "data.csv"):
        """
        Save a list of dictionaries to a CSV file in the session directory.
        The dictionaries should have keys that match the predefined headers.
        
        Args:
            project_id (str): The project ID
            session_id (str): The session ID
            data (list): List of dictionaries with keys matching the headers
            file_name (str): The name of the CSV file (default: "data.csv")
        
        Returns:
            str: The full path to the created CSV file
        
        Raises:
            TypeError: If data is not a list or doesn't contain dictionaries
            Exception: If file writing fails
        """
        # Define the standard headers for reaction data
        HEADERS = [
            "index",
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
        
        # Validate input
        if not isinstance(data, list):
            raise TypeError("data must be a list")
        
        if len(data) == 0:
            raise ValueError("data list cannot be empty")
        
        if not isinstance(data[0], dict):
            raise TypeError("data must contain dictionaries")
        
        # Ensure file name has .csv extension
        if not file_name.endswith('.csv'):
            file_name += '.csv'
        
        # Get the session directory (create if it doesn't exist)
        session_dir = self.create_session_directory(project_id, session_id)
        
        # Create the full file path
        file_path = os.path.join(session_dir, file_name)
        
        try:
            with open(file_path, 'w', newline='', encoding='utf-8') as csvfile:
                writer = csv.writer(csvfile)
                
                # Write the headers first
                writer.writerow(HEADERS)
                
                # Write each dictionary as a row, aligned with headers
                for row_dict in data:
                    row = []
                    for header in HEADERS:
                        # Get the value for this header, default to empty string if not found
                        value = row_dict.get(header, '')
                        row.append(value)
                    writer.writerow(row)
            
            return file_path
        except Exception as e:
            raise Exception(f"Failed to save CSV file {file_name}: {str(e)}")

    def get_session_prompt_text_file(self, project_id: str, session_id: str, file_name: str = "master_prompt.txt"):
        """
        Get the content of a session prompt text file.
        
        Args:
            project_id (str): The project ID
            session_id (str): The session ID
            file_name (str): The name of the text file (default: "master_prompt.txt")
        
        Returns:
            dict: Dictionary containing 'file_path', 'content', and 'exists' keys
        
        Raises:
            ValueError: If project or session directory doesn't exist
        """
        # Get the session directory path
        project_dir = self.get_project_directory(project_id)
        sessions_dir = os.path.join(project_dir, "sessions")
        session_dir = os.path.join(sessions_dir, f"session_{session_id}")
        
        # Check if session directory exists
        if not os.path.exists(session_dir):
            raise ValueError(f"Session directory for session {session_id} in project {project_id} does not exist")
        
        # Create the full file path
        file_path = os.path.join(session_dir, file_name)
        
        # Check if file exists and read content
        if os.path.exists(file_path):
            try:
                with open(file_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                return {
                    'file_path': file_path,
                    'content': content,
                    'exists': True
                }
            except Exception as e:
                raise Exception(f"Failed to read file {file_name}: {str(e)}")
        else:
            return {
                'file_path': file_path,
                'content': None,
                'exists': False
            }

    def get_session_csv_file(self, project_id: str, session_id: str, file_name: str = "data.csv"):
        """
        Get the content of a session CSV file as a list of dictionaries.
        
        Args:
            project_id (str): The project ID
            session_id (str): The session ID
            file_name (str): The name of the CSV file (default: "data.csv")
        
        Returns:
            dict: Dictionary containing 'file_path', 'data', 'headers', and 'exists' keys
        
        Raises:
            ValueError: If project or session directory doesn't exist
        """
        # Get the session directory path
        project_dir = self.get_project_directory(project_id)
        sessions_dir = os.path.join(project_dir, "sessions")
        session_dir = os.path.join(sessions_dir, f"session_{session_id}")
        
        # Check if session directory exists
        if not os.path.exists(session_dir):
            raise ValueError(f"Session directory for session {session_id} in project {project_id} does not exist")
        
        # Create the full file path
        file_path = os.path.join(session_dir, file_name)
        
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
                    'file_path': file_path,
                    'data': data,
                    'headers': headers,
                    'exists': True
                }
            except Exception as e:
                raise Exception(f"Failed to read CSV file {file_name}: {str(e)}")
        else:
            return {
                'file_path': file_path,
                'data': None,
                'headers': None,
                'exists': False
            }

    def merge_session_csv_to_project(self, project_id: str, session_id: str, session_file_name: str = "data.csv", project_file_name: str = "master.csv"):
        """
        Merge a session CSV file with the project's master CSV file.
        
        This method reads the session CSV file and appends its data to the project's master CSV file,
        ensuring no duplicate rows are added and maintaining data integrity.
        
        Args:
            project_id (str): The project ID
            session_id (str): The session ID
            session_file_name (str): The name of the session CSV file (default: "data.csv")
            project_file_name (str): The name of the project CSV file (default: "master.csv")
        
        Returns:
            dict: Dictionary containing merge results with keys:
                - 'success': Boolean indicating if merge was successful
                - 'session_file_path': Path to the session CSV file
                - 'project_file_path': Path to the project CSV file
                - 'rows_read': Number of rows read from session file
                - 'rows_added': Number of new rows added to project file
                - 'rows_skipped': Number of duplicate rows skipped
                - 'total_rows_after_merge': Total rows in project file after merge
        
        Raises:
            ValueError: If project or session directory doesn't exist
            Exception: If file operations fail
        """
        try:
            # Get session CSV data
            session_csv_result = self.get_session_csv_file(project_id, session_id, session_file_name)
            
            if not session_csv_result['exists']:
                raise ValueError(f"Session CSV file '{session_file_name}' does not exist for session {session_id}")
            
            session_data = session_csv_result['data']
            session_headers = session_csv_result['headers']
            
            if not session_data:
                return {
                    'success': True,
                    'session_file_path': session_csv_result['file_path'],
                    'project_file_path': None,
                    'rows_read': 0,
                    'rows_added': 0,
                    'rows_skipped': 0,
                    'total_rows_after_merge': 0,
                    'message': 'Session CSV file is empty, nothing to merge'
                }
            
            # Get project directory and master CSV path
            project_dir = self.get_project_directory(project_id)
            project_csv_path = os.path.join(project_dir, project_file_name)
            
            # Read existing project CSV data
            existing_project_data = []
            project_headers = session_headers  # Use session headers as default
            
            if os.path.exists(project_csv_path):
                try:
                    with open(project_csv_path, 'r', newline='', encoding='utf-8') as csvfile:
                        reader = csv.DictReader(csvfile)
                        project_headers = reader.fieldnames if reader.fieldnames else session_headers
                        
                        for row in reader:
                            existing_project_data.append(dict(row))
                except Exception as e:
                    raise Exception(f"Failed to read existing project CSV file: {str(e)}")
            
            # Create a set of existing rows for duplicate checking
            # Use a tuple of key-value pairs for hashable comparison
            existing_rows_set = set()
            for row in existing_project_data:
                # Create a hashable representation of the row
                row_tuple = tuple(sorted(row.items()))
                existing_rows_set.add(row_tuple)
            
            # Merge headers (union of both sets, preserving order)
            merged_headers = list(project_headers) if project_headers else []
            for header in session_headers:
                if header not in merged_headers:
                    merged_headers.append(header)
            
            # Filter out duplicate rows from session data
            new_rows = []
            rows_skipped = 0
            
            for session_row in session_data:
                # Ensure session row has all headers with empty string defaults
                normalized_row = {}
                for header in merged_headers:
                    normalized_row[header] = session_row.get(header, '')
                
                # Create hashable representation for comparison
                row_tuple = tuple(sorted(normalized_row.items()))
                
                if row_tuple not in existing_rows_set:
                    new_rows.append(normalized_row)
                    existing_rows_set.add(row_tuple)  # Add to set to prevent duplicates within session data
                else:
                    rows_skipped += 1
            
            # Normalize existing project data to have all headers
            normalized_existing_data = []
            for row in existing_project_data:
                normalized_row = {}
                for header in merged_headers:
                    normalized_row[header] = row.get(header, '')
                normalized_existing_data.append(normalized_row)
            
            # Calculate the next index value for new rows
            next_index = 0
            if normalized_existing_data:
                # Find the highest existing index
                max_index = -1
                for row in normalized_existing_data:
                    try:
                        current_index = int(row.get('index', '0'))
                        max_index = max(max_index, current_index)
                    except (ValueError, TypeError):
                        # If index is not a valid integer, treat as 0
                        pass
                next_index = max_index + 1
            
            # Update index values for new rows
            for i, new_row in enumerate(new_rows):
                new_row['index'] = str(next_index + i)
            
            # Combine existing and new data
            all_data = normalized_existing_data + new_rows
            
            # Write merged data back to project CSV
            try:
                with open(project_csv_path, 'w', newline='', encoding='utf-8') as csvfile:
                    writer = csv.writer(csvfile)
                    
                    # Write headers
                    writer.writerow(merged_headers)
                    
                    # Write all data rows
                    for row_dict in all_data:
                        row = []
                        for header in merged_headers:
                            value = row_dict.get(header, '')
                            row.append(value)
                        writer.writerow(row)
                
                return {
                    'success': True,
                    'session_file_path': session_csv_result['file_path'],
                    'project_file_path': project_csv_path,
                    'rows_read': len(session_data),
                    'rows_added': len(new_rows),
                    'rows_skipped': rows_skipped,
                    'total_rows_after_merge': len(all_data),
                    'message': f'Successfully merged {len(new_rows)} new rows from session CSV to project CSV'
                }
                
            except Exception as e:
                raise Exception(f"Failed to write merged data to project CSV: {str(e)}")
                
        except ValueError:
            raise
        except Exception as e:
            raise Exception(f"Error during CSV merge operation: {str(e)}")

# Singleton instance
file_system_handler = FileSystemHandler()
