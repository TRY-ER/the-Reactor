export interface Session {
  id: string;
  name: string;
  content: string | null;
  state: string | null;
  last_updated: string;
  worker_id: string;
}

class SessionService {
  private baseUrl: string;
  constructor() {
    this.baseUrl = process.env.REACT_APP_BACKEND_URL || 'http://localhost:8000';
  }

  async createSession(session: { name: string }): Promise<Session> {
    // Generate a unique ID for the session
    const sessionId = `session-${Date.now()}-${Math.random().toString(36).substr(2, 9)}`;
    const sessionData = {
      id: sessionId,
      name: session.name
    };

    const response = await fetch(`${this.baseUrl}/sessions/`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(sessionData),
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

  async getSession(sessionId: string): Promise<Session> {
    const response = await fetch(`${this.baseUrl}/sessions/${sessionId}`);
    if (!response.ok) throw new Error('Session not found');
    return response.json();
  }

  async updateSession(sessionId: string, session: { content?: string; state?: string; worker_id?: string }): Promise<Session> {
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
}

export const sessionService = new SessionService();
export default sessionService;
