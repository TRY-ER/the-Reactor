import React, { useState, useEffect } from 'react';
import { useParams, useNavigate } from 'react-router-dom';
import { MdTableChart, MdClose, MdDelete } from "react-icons/md";
import type { IconBaseProps } from 'react-icons';
import Sidebar from '../Sidebar';
import { projectService, sessionService, type ProjectResponse, type Session, type SessionCreate } from '../../services';
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
    name: '',
    query: '',
    content: [] as Array<Record<string, any>>,
    file_paths: [] as Array<Record<string, any>>
  });
  const [showSessionModal, setShowSessionModal] = useState(false);
  const [sessionError, setSessionError] = useState('');
  
  // Master CSV state
  const [masterCsvFile, setMasterCsvFile] = useState<{
    data: Array<Record<string, any>> | null;
    headers: string[] | null;
    exists: boolean;
    file_path: string;
    row_count: number;
    message?: string;
  } | null>(null);
  const [loadingMasterCsv, setLoadingMasterCsv] = useState(false);
  
  // Master CSV modal states
  const [masterCsvModal, setMasterCsvModal] = useState<{
    isOpen: boolean;
    data: Array<Record<string, any>>;
    headers: string[];
    isEditing: boolean;
  }>({ isOpen: false, data: [], headers: [], isEditing: false });
  
  // Master CSV editing state
  const [editedMasterCsvData, setEditedMasterCsvData] = useState<Array<Record<string, any>>>([]);
  const [isSavingMasterCsv, setIsSavingMasterCsv] = useState(false);

  useEffect(() => {
    if (id) {
      loadProject(id);
    }
  }, [id]);

  // Fetch sessions when project data is loaded
  useEffect(() => {
    if (project) {
      fetchSessions();
      fetchMasterCsv();
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
      const validSessions = all.filter(Boolean) as Session[];
      
      // Sort sessions by last_updated date in descending order (most recent first)
      const sortedSessions = validSessions.sort((a, b) => {
        const dateA = a.last_updated ? new Date(a.last_updated).getTime() : 0;
        const dateB = b.last_updated ? new Date(b.last_updated).getTime() : 0;
        return dateB - dateA; // Descending order (newest first)
      });
      
      setSessions(sortedSessions);
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
    
    if (!project?.id) {
      setSessionError('Project ID is required');
      return;
    }
    
    try {
      const sessionData = {
        ...newSession,
        project_id: project.id
      };
      const created = await sessionService.createSession(sessionData);
      setSessions(prev => [created, ...prev]); // Add new session at the beginning
      setNewSession({ 
        name: '', 
        query: '', 
        content: [], 
        file_paths: [] 
      });
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

  // Master CSV functions
  const fetchMasterCsv = async () => {
    if (!project?.id) return;
    
    try {
      setLoadingMasterCsv(true);
      const result = await projectService.getProjectMasterCsv(project.id);
      setMasterCsvFile({
        data: result.data,
        headers: result.headers,
        exists: result.exists,
        file_path: result.file_path,
        row_count: result.row_count,
        message: result.message
      });
    } catch (error) {
      console.warn('Failed to fetch master CSV:', error);
      setMasterCsvFile(null);
    } finally {
      setLoadingMasterCsv(false);
    }
  };

  const openMasterCsvModal = (data: Array<Record<string, any>>, headers: string[]) => {
    setMasterCsvModal({ isOpen: true, data, headers, isEditing: false });
    setEditedMasterCsvData(JSON.parse(JSON.stringify(data))); // Deep copy for editing
  };

  const closeMasterCsvModal = () => {
    setMasterCsvModal({ isOpen: false, data: [], headers: [], isEditing: false });
    setEditedMasterCsvData([]);
    setIsSavingMasterCsv(false);
  };

  const toggleMasterCsvEditing = () => {
    setMasterCsvModal(prev => ({ ...prev, isEditing: !prev.isEditing }));
  };

  const handleMasterCsvCellChange = (rowIndex: number, column: string, value: string) => {
    setEditedMasterCsvData(prev => {
      const newData = [...prev];
      newData[rowIndex] = { ...newData[rowIndex], [column]: value };
      return newData;
    });
  };

  const deleteMasterCsvRow = (rowIndex: number) => {
    setEditedMasterCsvData(prev => {
      const newData = [...prev];
      newData.splice(rowIndex, 1);
      return newData;
    });
  };

  const saveMasterCsvChanges = async () => {
    if (!project?.id) return;
    
    try {
      setIsSavingMasterCsv(true);
      
      // Save the edited data to the backend
      const result = await projectService.saveProjectMasterCsv(project.id, editedMasterCsvData);
      
      // Update the local state with the saved data
      setMasterCsvModal(prev => ({ ...prev, data: editedMasterCsvData, isEditing: false }));
      setMasterCsvFile(prev => prev ? { 
        ...prev, 
        data: editedMasterCsvData,
        row_count: editedMasterCsvData.length 
      } : null);
      
      console.log('Master CSV data saved successfully:', result);
      
    } catch (error) {
      console.error('Failed to save master CSV data:', error);
      alert('Failed to save master CSV data. Please try again.');
    } finally {
      setIsSavingMasterCsv(false);
    }
  };

  const resetMasterCsvChanges = () => {
    setEditedMasterCsvData(JSON.parse(JSON.stringify(masterCsvModal.data)));
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
                  <button 
                    onClick={handleEdit} 
                    className="btn-icon btn-primary" 
                    data-tooltip="Edit Project"
                  >
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                      <path d="M11 4H4a2 2 0 0 0-2 2v14a2 2 0 0 0 2 2h14a2 2 0 0 0 2-2v-7"/>
                      <path d="m18.5 2.5a2.121 2.121 0 0 1 3 3L12 15l-4 1 1-4 9.5-9.5z"/>
                    </svg>
                  </button>
                  <button 
                    onClick={handleDelete} 
                    className="btn-icon btn-danger" 
                    data-tooltip="Delete Project"
                  >
                    <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2" strokeLinecap="round" strokeLinejoin="round">
                      <polyline points="3,6 5,6 21,6"/>
                      <path d="m19,6v14a2,2 0 0,1 -2,2H7a2,2 0 0,1 -2,-2V6m3,0V4a2,2 0 0,1 2,-2h4a2,2 0 0,1 2,2v2"/>
                      <line x1="10" y1="11" x2="10" y2="17"/>
                      <line x1="14" y1="11" x2="14" y2="17"/>
                    </svg>
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
                  <p>{project.description}</p>
                </div>
                <div className="project-metadata">
                  <div className="metadata-item">
                    <label>Project ID:</label>
                    <span>{project.id}</span>
                  </div>
                </div>
                {!isEditing && (
                  <>
                    {/* Master CSV Section */}
                    <div className="project-master-csv-section">
                      <div className="master-csv-header">
                        <h2>Project Master CSV</h2>
                        <span className="master-csv-subtitle">Consolidated project data</span>
                      </div>
                      {loadingMasterCsv ? (
                        <div className="master-csv-loading">
                          <p>Loading master CSV...</p>
                        </div>
                      ) : masterCsvFile?.exists ? (
                        <div className="master-csv-item">
                          <div className="master-csv-icon">
                            {React.createElement(MdTableChart as React.ComponentType<IconBaseProps>, { 
                              size: 24,
                              title: "Master CSV Data File" 
                            })}
                          </div>
                          <div className="master-csv-info">
                            <h4>Master Data</h4>
                            <div className="master-csv-details">
                              <span className="master-csv-size">
                                {masterCsvFile.row_count} rows
                              </span>
                              {masterCsvFile.headers && (
                                <span className="master-csv-columns">
                                  {masterCsvFile.headers.length} columns
                                </span>
                              )}
                            </div>
                          </div>
                          <div className="master-csv-actions">
                            {masterCsvFile.data && masterCsvFile.headers && (
                              <>
                                <button 
                                  className="master-csv-btn view-btn"
                                  onClick={() => openMasterCsvModal(masterCsvFile.data!, masterCsvFile.headers!)}
                                  title="View master CSV data table"
                                >
                                  View Table
                                </button>
                                <button 
                                  className="master-csv-btn download-btn"
                                  onClick={() => {
                                    // Convert data back to CSV format
                                    const csvContent = [
                                      masterCsvFile.headers!.join(','),
                                      ...masterCsvFile.data!.map(row => 
                                        masterCsvFile.headers!.map(header => 
                                          JSON.stringify(row[header] || '')
                                        ).join(',')
                                      )
                                    ].join('\n');
                                    
                                    const blob = new Blob([csvContent], { type: 'text/csv' });
                                    const url = URL.createObjectURL(blob);
                                    const a = document.createElement('a');
                                    a.href = url;
                                    a.download = 'master.csv';
                                    document.body.appendChild(a);
                                    a.click();
                                    document.body.removeChild(a);
                                    URL.revokeObjectURL(url);
                                  }}
                                  title="Download master CSV file"
                                >
                                  Download
                                </button>
                              </>
                            )}
                          </div>
                        </div>
                      ) : (
                        <div className="master-csv-empty">
                          <p>No master CSV file found.</p>
                          <p className="master-csv-hint">
                            {masterCsvFile?.message || "Merge session data or create a master CSV to get started."}
                          </p>
                        </div>
                      )}
                    </div>
                    
                    {/* Sessions Section */}
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
                            <li key={sess.id} className="session-item">
                              <div 
                                className="session-container" 
                                onClick={() => navigate(`/project/${project.id}/session/${sess.id}`)}
                              >
                                <div className="session-main-content">
                                  <span className="session-name">{sess.name}</span>
                                  <span className={`session-state-badge ${sess.state || 'unknown'}`}>
                                    {sess.state || 'unknown'}
                                  </span>
                                </div>
                                <div className="session-details">
                                  <span className="session-id">ID: {sess.id}</span>
                                  {sess.last_updated && (
                                    <span className="session-updated">
                                      Updated: {new Date(sess.last_updated).toLocaleDateString()}
                                    </span>
                                  )}
                                </div>
                              </div>
                            </li>
                          ))}
                        </ul>
                      )}
                    </div>
                  </>
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
                <div className="form-group">
                  <label htmlFor="sessionQuery">Query (Optional)</label>
                  <textarea
                    id="sessionQuery"
                    name="query"
                    value={newSession.query}
                    onChange={handleNewSessionChange}
                    placeholder="Enter chemical reaction query (e.g., 8 Fe + S₈ → 8 FeS)"
                    rows={3}
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

      {/* Master CSV Table Modal */}
      {masterCsvModal.isOpen && (
        <div className="modal-overlay" onClick={closeMasterCsvModal}>
          <div className="modal-content csv-modal" onClick={(e) => e.stopPropagation()}>
            <div className="modal-header">
              <div className="modal-title">
                <h2>Master CSV Data Table</h2>
                <span className="csv-info">
                  {masterCsvModal.data.length} rows × {masterCsvModal.headers.length} columns
                </span>
              </div>
              <div className="modal-actions">
                <button
                  className={`modal-btn edit-toggle-btn ${masterCsvModal.isEditing ? 'active' : ''}`}
                  onClick={toggleMasterCsvEditing}
                  disabled={isSavingMasterCsv}
                >
                  {masterCsvModal.isEditing ? 'View Mode' : 'Edit Mode'}
                </button>
                <button 
                  className="modal-close-btn"
                  onClick={closeMasterCsvModal}
                  title="Close"
                >
                  {React.createElement(MdClose as React.ComponentType<IconBaseProps>, { size: 24 })}
                </button>
              </div>
            </div>
            <div className="modal-body csv-modal-body">
              <div className="csv-table-container">
                <table className="csv-table">
                  <thead>
                    <tr>
                      {masterCsvModal.headers.map((header, index) => (
                        <th key={index} className="csv-header">
                          {header}
                        </th>
                      ))}
                      {masterCsvModal.isEditing && (
                        <th className="csv-header csv-actions-header">
                          Actions
                        </th>
                      )}
                    </tr>
                  </thead>
                  <tbody>
                    {(masterCsvModal.isEditing ? editedMasterCsvData : masterCsvModal.data).map((row, rowIndex) => (
                      <tr key={rowIndex} className="csv-row">
                        {masterCsvModal.headers.map((header, colIndex) => (
                          <td key={colIndex} className="csv-cell">
                            {masterCsvModal.isEditing ? (
                              <input
                                type="text"
                                className="csv-cell-input"
                                value={editedMasterCsvData[rowIndex]?.[header] || ''}
                                onChange={(e) => handleMasterCsvCellChange(rowIndex, header, e.target.value)}
                              />
                            ) : (
                              <span className="csv-cell-value">
                                {row[header] || ''}
                              </span>
                            )}
                          </td>
                        ))}
                        {masterCsvModal.isEditing && (
                          <td className="csv-cell csv-actions-cell">
                            <button
                              className="csv-delete-btn"
                              onClick={() => deleteMasterCsvRow(rowIndex)}
                              title="Delete this row"
                              disabled={isSavingMasterCsv}
                            >
                              {React.createElement(MdDelete as React.ComponentType<IconBaseProps>, { size: 16 })}
                            </button>
                          </td>
                        )}
                      </tr>
                    ))}
                  </tbody>
                </table>
              </div>
            </div>
            <div className="modal-footer">
              {masterCsvModal.isEditing && (
                <div className="edit-actions">
                  <button
                    className="modal-btn save-btn"
                    onClick={saveMasterCsvChanges}
                    disabled={isSavingMasterCsv}
                  >
                    {isSavingMasterCsv ? 'Saving...' : 'Save Changes'}
                  </button>
                  <button
                    className="modal-btn reset-btn"
                    onClick={resetMasterCsvChanges}
                    disabled={isSavingMasterCsv}
                  >
                    Reset
                  </button>
                </div>
              )}
              <button
                className="modal-btn download-btn"
                onClick={() => {
                  const dataToDownload = masterCsvModal.isEditing ? editedMasterCsvData : masterCsvModal.data;
                  const csvContent = [
                    masterCsvModal.headers.join(','),
                    ...dataToDownload.map(row => 
                      masterCsvModal.headers.map(header => 
                        JSON.stringify(row[header] || '')
                      ).join(',')
                    )
                  ].join('\n');
                  
                  const blob = new Blob([csvContent], { type: 'text/csv' });
                  const url = URL.createObjectURL(blob);
                  const a = document.createElement('a');
                  a.href = url;
                  a.download = 'master.csv';
                  document.body.appendChild(a);
                  a.click();
                  document.body.removeChild(a);
                  URL.revokeObjectURL(url);
                }}
              >
                Download CSV
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default ProjectDetails;
