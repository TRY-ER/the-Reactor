import asyncio
import json
import uuid
import logging
from datetime import datetime
from typing import Dict, Optional, AsyncGenerator
from sqlalchemy.orm import Session
from collections import deque
from file_system_handler import FileSystemHandler

from db import SessionLocal, SessionObj
from file_system_handler import file_system_handler
import sys
import os
import time

# Set up logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Add the path to access the openai_runner_service
# Current file is in: backend/api_endpoint/session_runner.py
# Target import is in: backend/app/observers/master/openai_runner_service.py
# We need to go up one level to backend/, then into app/observers/master/
backend_root = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, backend_root)

try:
    from app.observers.master.openai_runner_service import run_service_agent
    logger.info("Successfully imported openai_runner_service")
except ImportError as e:
    logger.warning(f"Could not import openai_runner_service: {e}")
    logger.info("Using mock function for testing")
    # Mock function for testing
    async def run_service_agent(**kwargs):
        yield {"type": "info", "message": "Mock service started"}
        await asyncio.sleep(1)
        yield {"type": "data", "message": "Mock result 1"}
        await asyncio.sleep(1) 
        yield {"type": "data", "message": "Mock result 2"}
        yield {"final_csv": [{"mock": "data"}]}

class SessionWorker:
    """Individual worker that runs a session"""
    
    def __init__(self, session_id: str, worker_id: str, file_system_handler: FileSystemHandler):
        self.session_id = session_id
        self.worker_id = worker_id
        self.task: Optional[asyncio.Task] = None
        self.status = "idle"  # idle, running, completed, error
        self.created_at = datetime.utcnow()
        # get the project id from the session id from database
        db: Session = SessionLocal()
        try:
            session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
            if not session:
                raise ValueError(f"Session with id {session_id} not found")
            self.project_id = session.project_id
        finally:
            db.close()
        
        # Queue for SSE streaming
        self.content_queue: asyncio.Queue = asyncio.Queue()
        self.is_complete = False
        self.error_state = None
        self.streaming_clients = set()  # Track active streaming clients
        self.file_system_handler = file_system_handler
        
    async def run(self, **service_params):
        """Run the service agent and save results to database and queue"""
        self.status = "running"
         
        try:
            # Update session status in database
            await self._update_session_status("running")
            
            # Run the service agent and save each yield to database and queue
            async for result in run_service_agent(**service_params):
                if type(result) == str:
                    # Save to database (existing functionality)
                    await self._save_result_to_session(result)
                    
                    # Add to queue for streaming
                    await self._add_to_queue(result)
                if type(result) == dict:
                    if "final_csv" in result:
                        self.file_system_handler.save_data_to_session_csv(
                            self.project_id, self.session_id, result["final_csv"]
                        )
                   
                    if "master_prompt" in result:
                        self.file_system_handler.add_session_prompt_text_file(
                            self.project_id, self.session_id, {
                                "file_name": "master_prompt.txt",
                                "text_obj": result["master_prompt"]
                            } 
                        )
                    self.file_system_handler.add_file_path_to_session(
                        self.session_id, result 
                    )
                
            self.status = "completed"
            self.is_complete = True
            await self._update_session_status("completed")
            
            # Send completion signal to queue
            await self.content_queue.put({"_type": "completion", "status": "completed"})
            
        except Exception as e:
            # Debug block starts

            import traceback
            traceback.print_exc()

            # Debug block ends

            self.status = "error"
            self.error_state = str(e)
            self.is_complete = True
            
            error_result = {
                "type": "error",
                "message": str(e),
                "timestamp": datetime.utcnow().isoformat()
            }
            await self._save_result_to_session(error_result)
            await self._update_session_status("error")
            
            # Send error signal to queue
            await self.content_queue.put({"_type": "completion", "status": "error", "error": str(e)})
            raise
    
    async def _add_to_queue(self, result):
        """Add result to the content queue for streaming"""
        try:
            # Ensure result is a dictionary
            if not isinstance(result, dict):
                result = {
                    "type": "info", 
                    "content": str(result),
                    "timestamp": datetime.utcnow().isoformat()
                }
            
            # Add timestamp and worker_id if not present
            if "timestamp" not in result:
                result["timestamp"] = datetime.utcnow().isoformat()
            if "worker_id" not in result:
                result["worker_id"] = self.worker_id
            
            # Put result in queue for streaming
            await self.content_queue.put(result)
            
        except Exception as e:
            print(f"Error adding result to queue for session {self.session_id}: {e}")

    async def stream_content(self, client_id: Optional[str] = None) -> AsyncGenerator[dict, None]:
        """Stream content from the queue for SSE"""
        if client_id:
            self.streaming_clients.add(client_id)
        
        try:
            while True:
                try:
                    # Wait for content with timeout to allow for graceful shutdown
                    content = await asyncio.wait_for(self.content_queue.get(), timeout=1.0)
                    yield content
                    
                    # If completion event, stop streaming
                    if content.get("_type") == "completion":
                        break
                        
                except asyncio.TimeoutError:
                    # Check if worker is complete and queue is empty
                    if self.is_complete and self.content_queue.empty():
                        # Send final completion event if not already sent
                        yield {
                            "_type": "completion", 
                            "status": self.status,
                            "timestamp": datetime.utcnow().isoformat()
                        }
                        break
                    continue
                    
        except Exception as e:
            # Send error as final event
            yield {
                "_type": "completion",
                "status": "error",
                "error": f"Streaming error: {str(e)}",
                "timestamp": datetime.utcnow().isoformat()
            }
        finally:
            if client_id:
                self.streaming_clients.discard(client_id)
    
    async def _save_result_to_session(self, result):
        """Save a result to the session's content in database"""
        try:
            # Parse the result if it's a string from the stream
            if isinstance(result, str):
                # Check if it's a stream format like "data:{json}\n\n"
                if result.startswith("data:") and result.endswith("\n\n"):
                    # Extract JSON from the stream format
                    json_str = result[5:-2]  # Remove "data:" prefix and "\n\n" suffix
                    try:
                        result = json.loads(json_str)
                    except json.JSONDecodeError as e:
                        print(f"Failed to parse stream data: {e}")
                        # Create a fallback dict
                        result = {
                            "type": "error",
                            "content": f"Failed to parse stream data: {result}",
                            "timestamp": datetime.utcnow().isoformat()
                        }
                else:
                    # If it's a plain string, wrap it in a dict
                    result = {
                        "type": "info",
                        "content": str(result),
                        "timestamp": datetime.utcnow().isoformat()
                    }
            
            # Ensure result is a dictionary
            if not isinstance(result, dict):
                result = {
                    "type": "info", 
                    "content": str(result),
                    "timestamp": datetime.utcnow().isoformat()
                }
            
            # Add timestamp and worker_id to result
            result["timestamp"] = datetime.utcnow().isoformat()
            result["worker_id"] = self.worker_id
            
            # Use the file system handler to add content
            file_system_handler.add_content_to_session(self.session_id, result)
            
        except Exception as e:
            print(f"Error saving result to session {self.session_id}: {e}")
    
    async def _update_session_status(self, status: str):
        """Update session status in database"""
        db: Session = SessionLocal()
        try:
            session = db.query(SessionObj).filter(SessionObj.id == self.session_id).first()
            if session:
                session.state = status
                session.last_updated = datetime.utcnow()
                db.commit()
        except Exception as e:
            db.rollback()
            print(f"Error updating session status: {e}")
        finally:
            db.close()


class SessionRunner:
    """Manages session workers and their lifecycle"""
    
    def __init__(self, base_file_path: str):
        self.active_workers: Dict[str, SessionWorker] = {}
        self.worker_tasks: Dict[str, asyncio.Task] = {}
        self.file_system_handler = FileSystemHandler(base_file_path)

    async def start_session_worker(
        self, 
        session_id: str, 
        model_type: str,
        model_name: str,
        model_key: str,
        agent_name: str = "master_openai_oss_agent",
        query: Optional[str] = None  # Optional - will use session's query if not provided
    ) -> str:
        """Start a new worker for a session"""
        
        # Check if session exists and get query
        db: Session = SessionLocal()
        try:
            session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
            if not session:
                raise ValueError(f"Session with id {session_id} not found")
            
            # Use provided query or session's stored query
            if query is not None:
                # Update session with new query
                session.query = query
                db.commit()
                final_query = query
            else:
                # Use existing session query
                final_query = session.query or ""
                
            if not final_query.strip():
                raise ValueError("No query provided and session has no stored query")
                
        finally:
            db.close()
        
        # Clear previous session content before starting new worker
        logger.info(f"Clearing previous content for session {session_id}")
        try:
            file_system_handler.clear_session_content(session_id)
            logger.info(f"Successfully cleared content for session {session_id}")
        except Exception as e:
            logger.warning(f"Error clearing session content for {session_id}: {e}")
        
        # Check if worker is already running for this session
        if session_id in self.active_workers:
            existing_worker = self.active_workers[session_id]
            if existing_worker.status == "running":
                raise ValueError(f"Worker already running for session {session_id}")
        
        # Create new worker
        worker_id = str(uuid.uuid4())
        worker = SessionWorker(session_id, worker_id, self.file_system_handler)
        
        # Prepare service parameters
        service_params = {
            "model_type": model_type,
            "model_name": model_name,
            "model_key": model_key,
            "agent_name": agent_name,
            "query": final_query
        }
        
        # Update session with worker_id in database
        await self._assign_worker_to_session(session_id, worker_id)
        
        # Start the worker task
        task = asyncio.create_task(worker.run(**service_params))
        worker.task = task
        
        # Store worker and task
        self.active_workers[session_id] = worker
        self.worker_tasks[worker_id] = task
        
        # Add cleanup callback
        task.add_done_callback(lambda t: self._cleanup_worker(session_id, worker_id))
        
        return worker_id
    
    async def _assign_worker_to_session(self, session_id: str, worker_id: str):
        """Assign worker_id to session in database"""
        db: Session = SessionLocal()
        try:
            session = db.query(SessionObj).filter(SessionObj.id == session_id).first()
            if session:
                session.worker_id = worker_id
                session.state = "assigned"
                session.last_updated = datetime.utcnow()
                db.commit()
        except Exception as e:
            db.rollback()
            print(f"Error assigning worker to session: {e}")
        finally:
            db.close()

    async def stream_session_content(self, session_id: str, client_id: Optional[str] = None) -> AsyncGenerator[dict, None]:
        """Stream content for a session via Server-Sent Events"""
        logger.info(f"Starting SSE stream for session {session_id}, client {client_id}")
        
        # Check if there's an active worker for this session
        if session_id in self.active_workers:
            worker = self.active_workers[session_id]
            logger.info(f"Found active worker {worker.worker_id} for session {session_id}")
            
            # Stream from the active worker's queue
            async for content in worker.stream_content(client_id):
                yield content
        else:
            # No active worker - send any existing content from database and complete
            logger.info(f"No active worker for session {session_id}, sending existing content")
            
            try:
                # Get existing content from database
                content_list = file_system_handler.get_session_content(session_id)
                
                # Send existing content items
                for item in content_list:
                    if isinstance(item, dict):
                        # Add streaming metadata
                        item["_streaming"] = "historical"
                        yield item
                
                # Send completion event
                yield {
                    "_type": "completion",
                    "status": "no_active_worker",
                    "message": f"No active worker found for session {session_id}",
                    "timestamp": datetime.utcnow().isoformat()
                }
                
            except Exception as e:
                logger.error(f"Error streaming content for session {session_id}: {e}")
                yield {
                    "_type": "completion",
                    "status": "error",
                    "error": f"Error retrieving content: {str(e)}",
                    "timestamp": datetime.utcnow().isoformat()
                }
    
    def _cleanup_worker(self, session_id: str, worker_id: str):
        """Clean up completed worker"""
        if session_id in self.active_workers:
            del self.active_workers[session_id]
        
        if worker_id in self.worker_tasks:
            del self.worker_tasks[worker_id]
    
    def get_worker_status(self, session_id: str) -> Optional[dict]:
        """Get status of worker for a session"""
        if session_id in self.active_workers:
            worker = self.active_workers[session_id]
            return {
                "session_id": session_id,
                "worker_id": worker.worker_id,
                "status": worker.status,
                "created_at": worker.created_at.isoformat()
            }
        return None
    
    def get_all_workers_status(self) -> Dict[str, dict]:
        """Get status of all active workers"""
        return {
            session_id: {
                "worker_id": worker.worker_id,
                "status": worker.status,
                "created_at": worker.created_at.isoformat()
            }
            for session_id, worker in self.active_workers.items()
        }
    
    async def stop_worker(self, session_id: str) -> bool:
        """Stop a worker for a session"""
        if session_id in self.active_workers:
            worker = self.active_workers[session_id]
            if worker.task and not worker.task.done():
                worker.task.cancel()
                try:
                    await worker.task
                except asyncio.CancelledError:
                    pass
                
                # Update session status
                await worker._update_session_status("cancelled")
                return True
        return False

    def get_streaming_status(self, session_id: str) -> dict:
        """Get streaming status for a session"""
        if session_id in self.active_workers:
            worker = self.active_workers[session_id]
            return {
                "has_active_worker": True,
                "worker_id": worker.worker_id,
                "status": worker.status,
                "is_complete": worker.is_complete,
                "queue_size": worker.content_queue.qsize(),
                "streaming_clients": list(worker.streaming_clients),
                "created_at": worker.created_at.isoformat()
            }
        else:
            return {
                "has_active_worker": False,
                "status": "no_worker",
                "is_complete": True,
                "queue_size": 0,
                "streaming_clients": [],
                "message": f"No active worker found for session {session_id}"
            }
