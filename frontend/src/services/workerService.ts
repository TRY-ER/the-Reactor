export interface StartWorkerRequest {
  model_name: string;
  query?: string;
  agent_name?: string;
}

export interface StartWorkerResponse {
  message: string;
  session_id: string;
  worker_id: string;
}

export interface WorkerStatus {
  session_id: string;
  worker_id: string;
  status: 'idle' | 'running' | 'completed' | 'error' | 'cancelled';
  created_at: string;
}

class WorkerService {
  private baseUrl: string;
  private modelMapping: Record<string, string> = {
    'KIMI-K2-Instruct': 'groq/moonshotai/kimi-k2-instruct',
    'GPT-OSS-120B': 'groq/openai/gpt-oss-120b',
    'GPT-OSS-20B': 'groq/openai/gpt-oss-20b'
  };

  constructor() {
    this.baseUrl = process.env.REACT_APP_BACKEND_URL || 'http://localhost:8000';
  }

  /**
   * Start an AI worker for a session to analyze chemical reactions
   * @param sessionId - The unique identifier of the session
   * @param request - Worker configuration including model name and optional query
   * @returns Promise containing worker ID and confirmation message
   */
  async startSessionWorker(sessionId: string, request: StartWorkerRequest): Promise<StartWorkerResponse> {
    // Convert display model name to backend model name if needed
    const backendModelName = this.getBackendModelName(request.model_name);
    
    const requestBody = {
      ...request,
      model_name: backendModelName
    };

    const response = await fetch(`${this.baseUrl}/sessions/${sessionId}/start-worker`, {
      method: 'POST',
      headers: { 'Content-Type': 'application/json' },
      body: JSON.stringify(requestBody),
    });

    if (!response.ok) {
      const errorText = await response.text();
      throw new Error(`Failed to start worker: ${response.status} ${errorText}`);
    }

    return response.json();
  }

  /**
   * Get the current status of a worker running for a session
   * @param sessionId - The unique identifier of the session
   * @returns Promise containing worker status information
   */
  async getSessionWorkerStatus(sessionId: string): Promise<WorkerStatus> {
    const response = await fetch(`${this.baseUrl}/sessions/${sessionId}/worker-status`);
    
    if (!response.ok) {
      if (response.status === 404) {
        throw new Error('No active worker found for this session');
      }
      const errorText = await response.text();
      throw new Error(`Failed to get worker status: ${response.status} ${errorText}`);
    }

    return response.json();
  }

  /**
   * Helper method to check if a worker is currently active (running or idle)
   * @param sessionId - The unique identifier of the session
   * @returns Promise<boolean> - true if worker is active, false otherwise
   */
  async isWorkerActive(sessionId: string): Promise<boolean> {
    try {
      const status = await this.getSessionWorkerStatus(sessionId);
      return status.status === 'running' || status.status === 'idle';
    } catch (error) {
      // If no worker found, it's not active
      return false;
    }
  }

  /**
   * Helper method to get supported model types
   * @returns Array of supported model types
   */
  getSupportedModelTypes(): string[] {
    return ['groq', 'openai', 'azure'];
  }

  /**
   * Get available model options for dropdown
   * @returns Array of model display names
   */
  getAvailableModels(): string[] {
    return Object.keys(this.modelMapping);
  }

  /**
   * Convert display model name to backend model name
   * @param displayName - Model display name (e.g., "GPT-OSS-120B")
   * @returns Backend model name (e.g., "groq/openai/gpt-oss-120b")
   */
  getBackendModelName(displayName: string): string {
    return this.modelMapping[displayName] || displayName;
  }

  /**
   * Convert backend model name to display name
   * @param backendName - Backend model name (e.g., "groq/openai/gpt-oss-120b")
   * @returns Display model name (e.g., "GPT-OSS-120B")
   */
  getDisplayModelName(backendName: string): string {
    const entry = Object.entries(this.modelMapping).find(([, value]) => value === backendName);
    return entry ? entry[0] : backendName;
  }

  /**
   * Validate model name format
   * @param modelName - Model name to validate (display name or backend format)
   * @returns boolean - true if valid format
   */
  validateModelName(modelName: string): boolean {
    // Check if it's a valid display name
    if (Object.keys(this.modelMapping).includes(modelName)) {
      return true;
    }
    
    // Check if it's a valid backend format
    const parts = modelName.split('/');
    if (parts.length < 2 || parts[0].trim() === '') {
      return false;
    }
    
    const modelType = parts[0].trim();
    return this.getSupportedModelTypes().includes(modelType);
  }
}

export const workerService = new WorkerService();
export default workerService;
