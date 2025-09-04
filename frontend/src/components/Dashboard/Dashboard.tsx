import React, { useState, useEffect } from 'react';
import { useNavigate, useSearchParams } from 'react-router-dom';
import Sidebar from '../Sidebar';
import { projectService, type ProjectResponse } from '../../services';
import './Dashboard.css';

interface Project {
  id: string;
  name: string;
  description: string;
  sessions: string[];
}

const Dashboard: React.FC = () => {
  const navigate = useNavigate();
  const [searchParams] = useSearchParams();
  const [activeSection, setActiveSection] = useState<string>('projects');
  const [sidebarCollapsed, setSidebarCollapsed] = useState<boolean>(false);
  const [showProjectModal, setShowProjectModal] = useState<boolean>(false);
  const [newProjectName, setNewProjectName] = useState<string>('');
  const [newProjectDescription, setNewProjectDescription] = useState<string>('');
  const [searchTerm, setSearchTerm] = useState<string>('');
  const [projects, setProjects] = useState<Project[]>([]);
  const [loading, setLoading] = useState<boolean>(true);
  const [error, setError] = useState<string>('');

  // Load projects on component mount
  useEffect(() => {
    loadProjects();
  }, []);

  // Handle URL search params for section navigation
  useEffect(() => {
    const section = searchParams.get('section');
    if (section && ['dashboard', 'settings', 'projects'].includes(section)) {
      setActiveSection(section);
    }
  }, [searchParams]);

  const loadProjects = async () => {
    try {
      setLoading(true);
      setError('');
      const projectsData = await projectService.getProjects();
      setProjects(projectsData);
    } catch (err) {
      setError('Failed to load projects. Please check if the backend is running.');
     
      console.error('Error loading projects:', err);
      // Fallback to demo data if backend is not available
      setProjects([
        {
          id: '1',
          name: 'Project Alpha',
          description: 'A cutting-edge application built with React',
          sessions: ['session1', 'session2']
        },
        {
          id: '2',
          name: 'Project Beta',
          description: 'Machine learning pipeline for data analysis',
          sessions: []
        },
        {
          id: '3',
          name: 'Project Gamma',
          description: 'Mobile application development',
          sessions: ['session3']
        }
      ]);
    } finally {
      setLoading(false);
    }
  };

  const handleMenuSelect = (menu: string) => {
    setActiveSection(menu);
    // Update URL to reflect the current section
    if (menu === 'projects') {
      navigate('/');
    } else {
      navigate(`/?section=${menu}`);
    }
  };

  const handleSidebarToggle = (collapsed: boolean) => {
    setSidebarCollapsed(collapsed);
  };

  const openProjectModal = () => {
    setShowProjectModal(true);
  };

  const closeProjectModal = () => {
    setShowProjectModal(false);
    setNewProjectName('');
    setNewProjectDescription('');
  };

  const handleCreateProject = async () => {
    if (newProjectName.trim() && newProjectDescription.trim()) {
      try {
        const newProject = await projectService.createProject({
          name: newProjectName.trim(),
          description: newProjectDescription.trim(),
          sessions: []
        });
        setProjects([...projects, newProject]);
        closeProjectModal();
        setError('');
      } catch (err) {
        setError('Failed to create project. Please try again.');
        console.error('Error creating project:', err);
      }
    }
  };

  // Filter projects based on search term
  const filteredProjects = projects.filter(project =>
    project.name.toLowerCase().includes(searchTerm.toLowerCase()) ||
    project.description.toLowerCase().includes(searchTerm.toLowerCase()) ||
    (project.sessions && project.sessions.join(',').toLowerCase().includes(searchTerm.toLowerCase()))
  );

  const handleCardClick = (projectId: string) => {
    navigate(`/project/${projectId}`);
  };

  const handleDeleteProject = async (projectId: string, projectName: string, e: React.MouseEvent) => {
    e.stopPropagation(); // Prevent card click
    const confirmed = window.confirm(`Are you sure you want to delete "${projectName}"?`);
    if (!confirmed) return;

    try {
      await projectService.deleteProject(projectId);
      setProjects(projects.filter(p => p.id !== projectId));
      setError('');
    } catch (err) {
      setError('Failed to delete project. Please try again.');
      console.error('Error deleting project:', err);
    }
  };

  const renderContent = () => {
    switch (activeSection) {
      case 'projects':
        return (
          <div className="content-section">
            <div className="section-header">
              <h1>Projects</h1>
              <button className="add-button" onClick={openProjectModal}>
                <span className="add-icon">+</span>
              </button>
            </div>
            <p>Welcome to the Projects section. Here you can manage all your projects.</p>
            {error && (
              <div className="error-banner">
                <p>{error}</p>
              </div>
            )}
            <div className="search-container">
              <input
                type="text"
                placeholder="Search projects by name, description, or state..."
                value={searchTerm}
                onChange={(e) => setSearchTerm(e.target.value)}
                className="search-input"
              />
              <span className="search-icon">⌕</span>
            </div>
            {loading ? (
              <div className="loading-container">
                <div className="loading-spinner"></div>
                <p>Loading projects...</p>
              </div>
            ) : (
              <div className="projects-grid">
                {filteredProjects.map((project) => (
                  <div 
                    key={project.id} 
                    className="project-card"
                    onClick={() => handleCardClick(project.id)}
                  >
                    <div className="project-card-header">
                      <h3>{project.name}</h3>
                      <div className="view-hint">◎</div>
                    </div>
                    <p>{project.description}</p>
                    <div className="project-card-footer">
                      <div className="project-sessions">
                        {project.sessions ? `${project.sessions.length} session${project.sessions.length !== 1 ? 's' : ''}` : '0 sessions'}
                      </div>
                      <div className="project-actions">
                        <button
                          className="delete-btn"
                          onClick={(e) => handleDeleteProject(project.id, project.name, e)}
                          title="Delete Project"
                        >
                          ⌦
                        </button>
                      </div>
                    </div>
                  </div>
                ))}
                {filteredProjects.length === 0 && searchTerm && (
                  <div className="no-results">
                    <p>No projects found matching "{searchTerm}"</p>
                  </div>
                )}
                {filteredProjects.length === 0 && !searchTerm && !loading && (
                  <div className="no-results">
                    <p>No projects available. Create your first project!</p>
                  </div>
                )}
              </div>
            )}
          </div>
        );
      case 'dashboard':
        return (
          <div className="content-section">
            <h1>Dashboard</h1>
            <p>Overview of your system metrics and analytics.</p>
            <div className="dashboard-widgets">
              <div className="widget">
                <h3>Total Projects</h3>
                <div className="widget-value">12</div>
              </div>
              <div className="widget">
                <h3>Active Tasks</h3>
                <div className="widget-value">8</div>
              </div>
              <div className="widget">
                <h3>Completed</h3>
                <div className="widget-value">24</div>
              </div>
            </div>
          </div>
        );
      case 'settings':
        return (
          <div className="content-section">
            <h1>Settings</h1>
            <p>Configure your application preferences.</p>
            <div className="settings-form">
              <div className="setting-item">
                <label>Theme</label>
                <select>
                  <option>Dark Green</option>
                  <option>Classic Dark</option>
                </select>
              </div>
              <div className="setting-item">
                <label>Notifications</label>
                <input type="checkbox" defaultChecked />
              </div>
            </div>
          </div>
        );
      default:
        return (
          <div className="content-section">
            <h1>Welcome</h1>
            <p>Select an option from the sidebar to get started.</p>
          </div>
        );
    }
  };

  return (
    <div className="dashboard">
      <Sidebar onMenuSelect={handleMenuSelect} onToggle={handleSidebarToggle} activeMenu={activeSection} />
      <div className={`dashboard-content ${sidebarCollapsed ? 'sidebar-collapsed' : ''}`}>
        {renderContent()}
      </div>
      
      {/* Project Modal */}
      {showProjectModal && (
        <div className="modal-overlay" onClick={closeProjectModal}>
          <div className="modal-content" onClick={(e) => e.stopPropagation()}>
            <div className="modal-header">
              <h2>Create New Project</h2>
              <button className="modal-close" onClick={closeProjectModal}>×</button>
            </div>
            <div className="modal-body">
              <div className="form-group">
                <label htmlFor="projectName">Project Name</label>
                <input
                  type="text"
                  id="projectName"
                  value={newProjectName}
                  onChange={(e) => setNewProjectName(e.target.value)}
                  placeholder="Enter project name"
                />
              </div>
              <div className="form-group">
                <label htmlFor="projectDescription">Description</label>
                <textarea
                  id="projectDescription"
                  value={newProjectDescription}
                  onChange={(e) => setNewProjectDescription(e.target.value)}
                  placeholder="Enter project description"
                  rows={4}
                />
              </div>
            </div>
            <div className="modal-footer">
              <button className="btn-secondary" onClick={closeProjectModal}>
                Cancel
              </button>
              <button className="btn-primary" onClick={handleCreateProject}>
                Create Project
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default Dashboard;
