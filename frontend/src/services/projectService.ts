// Project interface to match the backend API
export interface Project {
  id?: string;
  name: string;
  description: string;
  sessions: string[];
}

// API response interfaces
export interface ProjectResponse extends Project {
  id: string;
}

export interface ProjectListResponse {
  projects: ProjectResponse[];
  total?: number;
}

class ProjectService {
  private baseUrl: string;

  constructor() {
    this.baseUrl = process.env.REACT_APP_BACKEND_URL || 'http://localhost:8000';
  }

  // Create a new project
  async createProject(project: Omit<Project, 'id'>): Promise<ProjectResponse> {
    try {
      const response = await fetch(`${this.baseUrl}/projects/`, {
        method: 'POST',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(project),
      });

      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      return await response.json();
    } catch (error) {
      console.error('Error creating project:', error);
      throw error;
    }
  }

  // Get all projects with optional pagination
  async getProjects(skip: number = 0, limit: number = 100): Promise<ProjectResponse[]> {
    try {
      const url = new URL(`${this.baseUrl}/projects/`);
      if (skip > 0) url.searchParams.append('skip', skip.toString());
      if (limit !== 100) url.searchParams.append('limit', limit.toString());

      const response = await fetch(url.toString());

      if (!response.ok) {
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      const data = await response.json();
      // Handle both array response and object with projects array
      return Array.isArray(data) ? data : data.projects || [];
    } catch (error) {
      console.error('Error fetching projects:', error);
      throw error;
    }
  }

  // Get a single project by ID
  async getProject(projectId: string): Promise<ProjectResponse> {
    try {
      const response = await fetch(`${this.baseUrl}/projects/${projectId}`);

      if (!response.ok) {
        if (response.status === 404) {
          throw new Error('Project not found');
        }
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      return await response.json();
    } catch (error) {
      console.error('Error fetching project:', error);
      throw error;
    }
  }

  // Update a project
  async updateProject(projectId: string, project: Omit<Project, 'id'>): Promise<ProjectResponse> {
    try {
      const response = await fetch(`${this.baseUrl}/projects/${projectId}`, {
        method: 'PUT',
        headers: {
          'Content-Type': 'application/json',
        },
        body: JSON.stringify(project),
      });

      if (!response.ok) {
        if (response.status === 404) {
          throw new Error('Project not found');
        }
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      return await response.json();
    } catch (error) {
      console.error('Error updating project:', error);
      throw error;
    }
  }

  // Delete a project
  async deleteProject(projectId: string): Promise<{ ok: boolean }> {
    try {
      const response = await fetch(`${this.baseUrl}/projects/${projectId}`, {
        method: 'DELETE',
      });

      if (!response.ok) {
        if (response.status === 404) {
          throw new Error('Project not found');
        }
        throw new Error(`HTTP error! status: ${response.status}`);
      }

      return await response.json();
    } catch (error) {
      console.error('Error deleting project:', error);
      throw error;
    }
  }
}

// Export a singleton instance
export const projectService = new ProjectService();
export default projectService;
