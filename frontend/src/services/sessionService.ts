export interface Session {
  id: string;
  name: string;
  content: Array<Record<string, any>> | null;
  state: string | null;
  last_updated: string;
  worker_id: string;
  project_id: string;
  file_paths: Array<Record<string, any>> | null;
  query: string | null;
}

export interface SessionCreate {
  name: string;
  content?: Array<Record<string, any>>;
  state?: string;
  project_id: string;
  file_paths?: Array<Record<string, any>>;
  query?: string;
}

class SessionService {
  private baseUrl: string;
  constructor() {
    this.baseUrl = process.env.REACT_APP_BACKEND_URL || 'http://localhost:8000';
  }

  async createSession(session: SessionCreate): Promise<Session> {
    const response = await fetch(`${this.baseUrl}/sessions/`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(session),
    });
    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to create session: ${response.status} ${errorText}`);
    }
    return response.json();
  }

  async getSessions(skip = 0, limit = 100): Promise<Session[]> {
    const url = new URL(`${this.baseUrl}/sessions/`);
    if (skip > 0) url.searchParams.append('skip', skip.toString());
    if (limit !== 100) url.searchParams.append('limit', limit.toString());
    const response = await fetch(url.toString());
    if (!response.ok) throw new Error('Failed to fetch sessions');
    return response.json();
  }

  async getSessionsByProject(projectId: string, skip = 0, limit = 100): Promise<Session[]> {
    const url = new URL(`${this.baseUrl}/sessions/`);
    url.searchParams.append('project_id', projectId);
    if (skip > 0) url.searchParams.append('skip', skip.toString());
    if (limit !== 100) url.searchParams.append('limit', limit.toString());
    const response = await fetch(url.toString());
    if (!response.ok) throw new Error('Failed to fetch sessions for project');
    return response.json();
  }

  async getSession(sessionId: string): Promise<Session> {
    const response = await fetch(`${this.baseUrl}/sessions/${sessionId}`);
    if (!response.ok) throw new Error('Session not found');
    return response.json();
  }

  async updateSession(sessionId: string, session: { 
    content?: Array<Record<string, any>>; 
    state?: string; 
    worker_id?: string;
    file_paths?: Array<Record<string, any>>;
    query?: string;
  }): Promise<Session> {
    const response = await fetch(`${this.baseUrl}/sessions/${sessionId}`, {
      method: 'PUT',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(session),
    });
    if (!response.ok) throw new Error('Failed to update session');
    return response.json();
  }

  async deleteSession(sessionId: string): Promise<{ ok: boolean }> {
    const response = await fetch(`${this.baseUrl}/sessions/${sessionId}`, { method: 'DELETE' });
    if (!response.ok) throw new Error('Failed to delete session');
    return response.json();
  }

  async updateSessionQuery(sessionId: string, query: string): Promise<{ message: string; session_id: string; query: string }> {
    const response = await fetch(`${this.baseUrl}/sessions/${sessionId}/query`, {
      method: 'PUT',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({ query }),
    });
    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to update query: ${response.status} ${errorText}`);
    }
    return response.json();
  }

  async getSessionQuery(sessionId: string): Promise<{ session_id: string; query: string }> {
    const response = await fetch(`${this.baseUrl}/sessions/${sessionId}/query`);
    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to get query: ${response.status} ${errorText}`);
    }
    return response.json();
  }

  // Stop/cancel the worker for a session
  async stopSessionWorker(sessionId: string): Promise<{ message: string; session_id: string }> {
    const response = await fetch(`${this.baseUrl}/sessions/${sessionId}/stop-worker`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
    });
    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to stop worker: ${response.status} ${errorText}`);
    }
    return response.json();
  }

  // Stream session content via SSE
  async streamSessionContent(
    sessionId: string, 
    clientId?: string,
    onMessage?: (data: any) => void,
    onError?: (error: Error) => void,
    onComplete?: () => void
  ): Promise<EventSource> {
    const url = new URL(`${this.baseUrl}/sessions/${sessionId}/stream`);
    if (clientId) {
      url.searchParams.append('client_id', clientId);
    }

    const eventSource = new EventSource(url.toString());

    eventSource.onmessage = (event) => {
      try {
        const data = JSON.parse(event.data);
        
        // Handle completion event
        if (data._type === 'completion') {
          eventSource.close();
          if (data.status === 'error' && onError) {
            onError(new Error(data.error || 'Stream completed with error'));
          } else if (onComplete) {
            onComplete();
          }
          return;
        }

        // Call the message handler with parsed data
        if (onMessage) {
          onMessage(data);
        }
      } catch (error) {
        console.error('Error parsing SSE data:', error);
        if (onError) {
          onError(error as Error);
        }
      }
    };

    eventSource.onerror = (error) => {
      console.error('SSE connection error:', error);
      if (onError) {
        onError(new Error('SSE connection error'));
      }
    };

    return eventSource;
  }

  // Get session prompt text file
  async getSessionPromptFile(sessionId: string, fileName: string = "master_prompt.txt"): Promise<{
    session_id: string;
    project_id: string;
    file_path: string;
    content: string | null;
    exists: boolean;
    message?: string;
    detail?: string;
  }> {
    const url = new URL(`${this.baseUrl}/sessions/${sessionId}/prompt-file`);
    if (fileName !== "master_prompt.txt") {
      url.searchParams.append('file_name', fileName);
    }
    
    const response = await fetch(url.toString());
    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to get prompt file: ${response.status} ${errorText}`);
    }
    return response.json();
  }

  // Get session CSV file
  async getSessionCsvFile(sessionId: string, fileName: string = "data.csv"): Promise<{
    session_id: string;
    project_id: string;
    file_path: string;
    data: Array<Record<string, any>> | null;
    headers: string[] | null;
    exists: boolean;
    row_count: number;
    message?: string;
    detail?: string;
  }> {
    const url = new URL(`${this.baseUrl}/sessions/${sessionId}/csv-file`);
    if (fileName !== "data.csv") {
      url.searchParams.append('file_name', fileName);
    }
    
    const response = await fetch(url.toString());
    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to get CSV file: ${response.status} ${errorText}`);
    }
    return response.json();
  }

  // Save CSV data to session
  async saveSessionCsvData(sessionId: string, data: Array<Record<string, any>>, fileName: string = "data.csv"): Promise<{
    message: string;
    session_id: string;
    project_id: string;
    file_path: string;
    rows_saved: number;
    file_name: string;
  }> {
    const response = await fetch(`${this.baseUrl}/sessions/${sessionId}/save-csv`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify({
        data: data,
        file_name: fileName
      }),
    });
    
    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to save CSV data: ${response.status} ${errorText}`);
    }
    return response.json();
  }

  // Merge session CSV with project master CSV
  async mergeSessionCsvToProject(
    sessionId: string, 
    sessionFileName: string = "data.csv", 
    projectFileName: string = "master.csv"
  ): Promise<{
    message: string;
    session_id: string;
    project_id: string;
    merge_results: {
      success: boolean;
      session_file_path: string;
      project_file_path: string;
      rows_read: number;
      rows_added: number;
      rows_skipped: number;
      total_rows_after_merge: number;
      message: string;
    };
  }> {
    const url = new URL(`${this.baseUrl}/sessions/${sessionId}/merge-csv`);
    if (sessionFileName !== "data.csv") {
      url.searchParams.append('session_file_name', sessionFileName);
    }
    if (projectFileName !== "master.csv") {
      url.searchParams.append('project_file_name', projectFileName);
    }
    
    const response = await fetch(url.toString(), {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
    });
    
    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to merge CSV files: ${response.status} ${errorText}`);
    }
    return response.json();
  }
}

export const sessionService = new SessionService();
export default sessionService;
