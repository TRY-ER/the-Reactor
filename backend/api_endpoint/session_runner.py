import asyncio
import json
import uuid
import logging
from datetime import datetime
from typing import Dict, Optional
from sqlalchemy.orm import Session

from db import SessionLocal, SessionObj
from file_system_handler import file_system_handler
import sys
import os

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
    
    def __init__(self, session_id: str, worker_id: str):
        self.session_id = session_id
        self.worker_id = worker_id
        self.task: Optional[asyncio.Task] = None
        self.status = "idle"  # idle, running, completed, error
        self.created_at = datetime.utcnow()
        
    async def run(self, **service_params):
        """Run the service agent and save results to database"""
        self.status = "running"
        
        try:
            # Update session status in database
            await self._update_session_status("running")
            
            # Run the service agent and save each yield to database
            async for result in run_service_agent(**service_params):
                await self._save_result_to_session(result)
                
            self.status = "completed"
            await self._update_session_status("completed")
            
        except Exception as e:
            self.status = "error"
            error_result = {
                "type": "error",
                "message": str(e),
                "timestamp": datetime.utcnow().isoformat()
            }
            await self._save_result_to_session(error_result)
            await self._update_session_status("error")
            raise
    
    async def _save_result_to_session(self, result: dict):
        """Save a result to the session's content in database"""
        try:
            # Add timestamp to result
            if isinstance(result, dict):
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
    
    def __init__(self):
        self.active_workers: Dict[str, SessionWorker] = {}
        self.worker_tasks: Dict[str, asyncio.Task] = {}
    
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
        
        # Check if worker is already running for this session
        if session_id in self.active_workers:
            existing_worker = self.active_workers[session_id]
            if existing_worker.status == "running":
                raise ValueError(f"Worker already running for session {session_id}")
        
        # Create new worker
        worker_id = str(uuid.uuid4())
        worker = SessionWorker(session_id, worker_id)
        
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


# Global session runner instance
session_runner = SessionRunner()