import React, { useEffect, useState } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import Sidebar from '../Sidebar';
import { sessionService, type Session } from '../../services';
import './SessionDetails.css';

const SessionDetails: React.FC = () => {
  const { sessionId } = useParams<{ sessionId: string }>();
  const navigate = useNavigate();
  const [session, setSession] = useState<Session | null>(null);
  const [loading, setLoading] = useState(true);
  const [error, setError] = useState('');
  const [sidebarCollapsed, setSidebarCollapsed] = useState(false);
  const [isEditing, setIsEditing] = useState(false);
  const [editForm, setEditForm] = useState({
    content: '',
    state: '',
    worker_id: ''
  });

  useEffect(() => {
    if (sessionId) loadSession(sessionId);
  }, [sessionId]);

  const loadSession = async (id: string) => {
    try {
      setLoading(true);
      setError('');
      const data = await sessionService.getSession(id);
      setSession(data);
      setEditForm({
        content: data.content || '',
        state: data.state || '',
        worker_id: data.worker_id
      });
    } catch (err) {
      setError('Session not found');
    } finally {
      setLoading(false);
    }
  };

  const handleMenuSelect = (menu: string) => {
    if (menu === 'projects') navigate('/');
    else if (menu === 'dashboard') navigate('/?section=dashboard');
    else if (menu === 'settings') navigate('/?section=settings');
  };
  const handleSidebarToggle = (collapsed: boolean) => setSidebarCollapsed(collapsed);

  const handleEdit = () => setIsEditing(true);
  const handleCancelEdit = () => {
    setIsEditing(false);
    if (session) setEditForm({ content: session.content || '', state: session.state || '', worker_id: session.worker_id });
  };
  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement | HTMLSelectElement>) => {
    const { name, value } = e.target;
    setEditForm(prev => ({ ...prev, [name]: value }));
  };
  const handleSaveEdit = async () => {
    if (!session || !sessionId) return;
    try {
      setLoading(true);
      const updated = await sessionService.updateSession(sessionId, editForm);
      setSession(updated);
      setIsEditing(false);
      setError('');
    } catch {
      setError('Failed to update session');
    } finally {
      setLoading(false);
    }
  };
  const handleDelete = async () => {
    if (!session || !sessionId) return;
    if (!window.confirm('Delete this session?')) return;
    try {
      setLoading(true);
      await sessionService.deleteSession(sessionId);
      navigate(-1);
    } catch {
      setError('Failed to delete session');
      setLoading(false);
    }
  };

  return (
    <div className="session-details-page">
      <Sidebar onMenuSelect={handleMenuSelect} onToggle={handleSidebarToggle} activeMenu="projects" />
      <div className={`session-details-content ${sidebarCollapsed ? 'sidebar-collapsed' : ''}`}>
        <div className="breadcrumb">
          <span onClick={() => navigate('/')} className="breadcrumb-link">Projects</span>
          <span className="breadcrumb-separator">/</span>
          <span className="breadcrumb-link" onClick={() => navigate(-1)}>Project</span>
          <span className="breadcrumb-separator">/</span>
          <span className="breadcrumb-current">{session?.name || sessionId}</span>
        </div>
        <div className="session-details">
          {loading ? (
            <div className="session-details-loading"><div className="loading-spinner"></div><p>Loading session...</p></div>
          ) : error ? (
            <div className="session-details-error"><h2>Error</h2><p>{error}</p></div>
          ) : !session ? (
            <div className="session-details-error"><h2>Not Found</h2><p>Session not found.</p></div>
          ) : (
            <div className="session-details-content-inner">
              <div className="session-details-header">
                <button onClick={() => navigate(-1)} className="back-button">← Back</button>
                <div className="session-actions">
                  {!isEditing && <><button onClick={handleEdit} className="btn-primary">Edit</button><button onClick={handleDelete} className="btn-danger">Delete</button></>}
                </div>
              </div>
              {isEditing ? (
                <div className="edit-form">
                  <h1>Edit Session</h1>
                  <div className="form-group">
                    <label htmlFor="content">Content (JSON)</label>
                    <textarea id="content" name="content" value={editForm.content} onChange={handleInputChange} rows={6} placeholder="Enter JSON content" />
                  </div>
                  <div className="form-group">
                    <label htmlFor="state">State</label>
                    <input id="state" name="state" value={editForm.state} onChange={handleInputChange} placeholder="Enter state" />
                  </div>
                  <div className="form-group">
                    <label htmlFor="worker_id">Worker ID</label>
                    <input id="worker_id" name="worker_id" value={editForm.worker_id} onChange={handleInputChange} placeholder="Enter worker ID" />
                  </div>
                  <div className="form-actions">
                    <button onClick={handleCancelEdit} className="btn-secondary">Cancel</button>
                    <button onClick={handleSaveEdit} className="btn-primary">Save Changes</button>
                  </div>
                </div>
              ) : (
                <div className="session-info">
                  <h1>Session: {session.name}</h1>
                  <div className="session-meta">
                    <div><b>State:</b> {session.state || 'No state'}</div>
                    <div><b>Worker:</b> {session.worker_id}</div>
                    <div><b>Last Updated:</b> {session.last_updated}</div>
                  </div>
                  <div className="session-content">
                    <h2>Content</h2>
                    <pre>{session.content || 'No content'}</pre>
                  </div>
                </div>
              )}
            </div>
          )}
        </div>
      </div>
    </div>
  );
};

export default SessionDetails;
