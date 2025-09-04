import React, { useState, useEffect } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import Sidebar from '../Sidebar';
import { projectService, sessionService, type ProjectResponse, type Session } from '../../services';
import './ProjectDetails.css';

const ProjectDetails: React.FC = () => {
  const { id } = useParams<{ id: string }>();
  const navigate = useNavigate();
  const [project, setProject] = useState<ProjectResponse | null>(null);
  const [sessions, setSessions] = useState<Session[]>([]);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string>('');
  const [isEditing, setIsEditing] = useState<boolean>(false);
  const [sidebarCollapsed, setSidebarCollapsed] = useState<boolean>(false);
  const [editForm, setEditForm] = useState({
    name: '',
    description: ''
  });
  const [newSession, setNewSession] = useState({
    name: ''
  });
  const [showSessionModal, setShowSessionModal] = useState(false);
  const [sessionError, setSessionError] = useState('');

  useEffect(() => {
    if (id) {
      loadProject(id);
    }
  }, [id]);

  // Fetch sessions when project data is loaded
  useEffect(() => {
    if (project) {
      fetchSessions();
    }
  }, [project]);

  const handleMenuSelect = (menu: string) => {
    if (menu === 'projects') {
      navigate('/');
    } else if (menu === 'dashboard') {
      navigate('/?section=dashboard');
    } else if (menu === 'settings') {
      navigate('/?section=settings');
    }
  };

  const handleSidebarToggle = (collapsed: boolean) => {
    setSidebarCollapsed(collapsed);
  };

  const loadProject = async (projectId: string) => {
    try {
      setLoading(true);
      setError('');
      const projectData = await projectService.getProject(projectId);
      setProject(projectData);
      setEditForm({
        name: projectData.name,
        description: projectData.description
      });
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to load project');
    } finally {
      setLoading(false);
    }
  };

  const fetchSessions = async () => {
    if (!project || !project.sessions) return setSessions([]);
    try {
      const all = await Promise.all(
        project.sessions.map((sid: string) => sessionService.getSession(sid).catch(() => null))
      );
      setSessions(all.filter(Boolean) as Session[]);
    } catch {
      setSessions([]);
    }
  };

  const handleEdit = () => {
    setIsEditing(true);
  };

  const handleCancelEdit = () => {
    setIsEditing(false);
    if (project) {
      setEditForm({
        name: project.name,
        description: project.description
      });
    }
  };

  const handleSaveEdit = async () => {
    if (!project || !id) return;

    try {
      setLoading(true);
      const updatedProject = await projectService.updateProject(id, {
        name: editForm.name,
        description: editForm.description,
        sessions: project.sessions // keep sessions unchanged
      });
      setProject(updatedProject);
      setIsEditing(false);
      setError('');
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to update project');
    } finally {
      setLoading(false);
    }
  };

  const handleDelete = async () => {
    if (!project || !id) return;

    const confirmed = window.confirm(`Are you sure you want to delete "${project.name}"? This action cannot be undone.`);
    if (!confirmed) return;

    try {
      setLoading(true);
      await projectService.deleteProject(id);
      navigate('/');
    } catch (err) {
      setError(err instanceof Error ? err.message : 'Failed to delete project');
      setLoading(false);
    }
  };

  const handleNewSessionChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement>) => {
    const { name, value } = e.target;
    setNewSession(prev => ({ ...prev, [name]: value }));
  };

  const handleCreateSession = async (e: React.FormEvent) => {
    e.preventDefault();
    setSessionError('');
    try {
      const created = await sessionService.createSession(newSession);
      setSessions(prev => [...prev, created]);
      setNewSession({ name: '' });
      setShowSessionModal(false);

      // Update the project's sessions array to include the new session
      if (project && id) {
        const updatedSessions = [...project.sessions, created.id];
        await projectService.updateProject(id, {
          name: project.name,
          description: project.description,
          sessions: updatedSessions
        });
        // Update local project state
        setProject(prev => prev ? { ...prev, sessions: updatedSessions } : null);
      }
    } catch (error) {
      console.error('Error creating session:', error);
      setSessionError('Failed to create session');
    }
  };

  const handleInputChange = (e: React.ChangeEvent<HTMLInputElement | HTMLTextAreaElement | HTMLSelectElement>) => {
    const { name, value } = e.target;
    setEditForm(prev => ({
      ...prev,
      [name]: value
    }));
  };

  if (loading) {
    return (
      <div className="project-details-page">
        <Sidebar onMenuSelect={handleMenuSelect} onToggle={handleSidebarToggle} activeMenu="projects" />
        <div className={`project-details-content ${sidebarCollapsed ? 'sidebar-collapsed' : ''}`}>
          <div className="breadcrumb">
            <span onClick={() => navigate('/')} className="breadcrumb-link">Projects</span>
            <span className="breadcrumb-separator">/</span>
            <span className="breadcrumb-current">Loading...</span>
          </div>
          <div className="project-details">
            <div className="project-details-loading">
              <div className="loading-spinner"></div>
              <p>Loading project...</p>
            </div>
          </div>
        </div>
      </div>
    );
  }

  if (error) {
    return (
      <div className="project-details-page">
        <Sidebar onMenuSelect={handleMenuSelect} onToggle={handleSidebarToggle} activeMenu="projects" />
        <div className={`project-details-content ${sidebarCollapsed ? 'sidebar-collapsed' : ''}`}>
          <div className="breadcrumb">
            <span onClick={() => navigate('/')} className="breadcrumb-link">Projects</span>
            <span className="breadcrumb-separator">/</span>
            <span className="breadcrumb-current">Error</span>
          </div>
          <div className="project-details">
            <div className="project-details-error">
              <h2>Error</h2>
              <p>{error}</p>
              <button onClick={() => navigate('/')} className="btn-secondary">
                Back to Projects
              </button>
            </div>
          </div>
        </div>
      </div>
    );
  }

  if (!project) {
    return (
      <div className="project-details-page">
        <Sidebar onMenuSelect={handleMenuSelect} onToggle={handleSidebarToggle} activeMenu="projects" />
        <div className={`project-details-content ${sidebarCollapsed ? 'sidebar-collapsed' : ''}`}>
          <div className="breadcrumb">
            <span onClick={() => navigate('/')} className="breadcrumb-link">Projects</span>
            <span className="breadcrumb-separator">/</span>
            <span className="breadcrumb-current">Not Found</span>
          </div>
          <div className="project-details">
            <div className="project-details-error">
              <h2>Project Not Found</h2>
              <p>The requested project could not be found.</p>
              <button onClick={() => navigate('/')} className="btn-secondary">
                Back to Projects
              </button>
            </div>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="project-details-page">
      <Sidebar onMenuSelect={handleMenuSelect} onToggle={handleSidebarToggle} activeMenu="projects" />
      <div className={`project-details-content ${sidebarCollapsed ? 'sidebar-collapsed' : ''}`}>
        <div className="breadcrumb">
          <span onClick={() => navigate('/')} className="breadcrumb-link">Projects</span>
          <span className="breadcrumb-separator">/</span>
          <span className="breadcrumb-current">{project.id}</span>
        </div>
        <div className="project-details">
          <div className="project-details-header">
            <button onClick={() => navigate('/')} className="back-button">
              ← Back to Projects
            </button>
            <div className="project-actions">
              {!isEditing && (
                <>
                  <button onClick={handleEdit} className="btn-primary">
                    Edit Project
                  </button>
                  <button onClick={handleDelete} className="btn-danger">
                    Delete Project
                  </button>
                </>
              )}
            </div>
          </div>

          <div className="project-details-content-inner">
            {isEditing ? (
              <div className="edit-form">
                <h1>Edit Project</h1>
                <div className="form-group">
                  <label htmlFor="name">Project Name</label>
                  <input
                    type="text"
                    id="name"
                    name="name"
                    value={editForm.name}
                    onChange={handleInputChange}
                    placeholder="Enter project name"
                  />
                </div>
                <div className="form-group">
                  <label htmlFor="description">Description</label>
                  <textarea
                    id="description"
                    name="description"
                    value={editForm.description}
                    onChange={handleInputChange}
                    placeholder="Enter project description"
                    rows={6}
                  />
                </div>
                <div className="form-actions">
                  <button onClick={handleCancelEdit} className="btn-secondary">
                    Cancel
                  </button>
                  <button onClick={handleSaveEdit} className="btn-primary">
                    Save Changes
                  </button>
                </div>
              </div>
            ) : (
              <div className="project-info">
                <div className="project-header">
                  <h1>{project.name}</h1>
                  <div className="project-sessions-badge">
                    {project.sessions ? `${project.sessions.length} session${project.sessions.length !== 1 ? 's' : ''}` : '0 sessions'}
                  </div>
                </div>
                <div className="project-description">
                  <h2>About this project</h2>
                  <p>{project.description}</p>
                </div>
                <div className="project-metadata">
                  <div className="metadata-item">
                    <label>Project ID:</label>
                    <span>{project.id}</span>
                  </div>
                </div>
                {!isEditing && (
                  <div className="project-sessions-list">
                    <div className="sessions-header">
                      <h2>Sessions</h2>
                      <button className="add-session-btn" onClick={() => setShowSessionModal(true)} title="Add Session">
                        +
                      </button>
                    </div>
                    {sessions.length === 0 ? <p>No sessions for this project.</p> : (
                      <ul>
                        {sessions.map(sess => (
                          <li key={sess.id}>
                            <span className="session-link" onClick={() => navigate(`/project/${project.id}/session/${sess.id}`)}>{sess.name}</span>
                            <span className="session-id">{sess.id}</span>
                          </li>
                        ))}
                      </ul>
                    )}
                  </div>
                )}
              </div>
            )}
          </div>
        </div>
      </div>
      
      {/* Session Creation Modal */}
      {showSessionModal && (
        <div className="modal-overlay" onClick={() => setShowSessionModal(false)}>
          <div className="modal-content" onClick={(e) => e.stopPropagation()}>
            <div className="modal-header">
              <h2>Create New Session</h2>
              <button className="modal-close" onClick={() => setShowSessionModal(false)}>×</button>
            </div>
            <form onSubmit={handleCreateSession}>
              <div className="modal-body">
                <div className="form-group">
                  <label htmlFor="sessionName">Session Name</label>
                  <input
                    type="text"
                    id="sessionName"
                    name="name"
                    value={newSession.name}
                    onChange={handleNewSessionChange}
                    placeholder="Enter session name"
                    required
                  />
                </div>
              </div>
              <div className="modal-footer">
                <button type="button" className="btn-secondary" onClick={() => setShowSessionModal(false)}>
                  Cancel
                </button>
                <button type="submit" className="btn-primary">
                  Create Session
                </button>
              </div>
              {sessionError && <div className="error-banner">{sessionError}</div>}
            </form>
          </div>
        </div>
      )}
    </div>
  );
};

export default ProjectDetails;
