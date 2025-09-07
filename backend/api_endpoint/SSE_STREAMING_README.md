# Session Content Streaming with Server-Sent Events (SSE)

This implementation provides real-time streaming of session content using Server-Sent Events (SSE). When a session worker processes chemical reactions, the results are queued and streamed to connected clients in real-time.

## Architecture Overview

### Queue-Based Streaming
- Each `SessionWorker` maintains an `asyncio.Queue` for content streaming
- Results from `run_service_agent` are added to both the database and the streaming queue
- Multiple clients can connect and receive the same stream of events

### Components

1. **SessionWorker** (`session_runner.py`)
   - Manages individual session processing
   - Maintains content queue for streaming
   - Tracks connected streaming clients
   - Provides completion signals

2. **SessionRunner** (`session_runner.py`)
   - Manages multiple session workers
   - Provides streaming interface for sessions
   - Handles fallback to historical content

3. **SSE Endpoints** (`router.py`)
   - `/sessions/{session_id}/stream` - Main SSE streaming endpoint
   - `/sessions/{session_id}/streaming-status` - Get streaming status

## API Endpoints

### Stream Session Content
```
GET /sessions/{session_id}/stream?client_id={optional_client_id}
```

**Response**: Server-Sent Events stream with the following event types:

#### Info Events
```json
{
  "type": "info",
  "message": "Starting SMILES Agent for reaction 0",
  "timestamp": "2024-01-01T12:00:00.000Z",
  "worker_id": "worker-uuid",
  "aux": {
    "info_type": "smiles_init",
    "reaction_index": 0
  }
}
```

#### Data Events
```json
{
  "type": "data",
  "message": "Got SMILES values for reaction 0",
  "timestamp": "2024-01-01T12:00:00.000Z",
  "worker_id": "worker-uuid",
  "aux": {
    "info_type": "smiles_update",
    "reaction_index": 0,
    "returnable": {
      "values": ["CCO", "C=O"],
      "invalid_count": 1
    }
  }
}
```

#### Completion Events
```json
{
  "_type": "completion",
  "status": "completed"
}
```

#### Error Events
```json
{
  "_type": "completion",
  "status": "error",
  "error": "Error message"
}
```

### Get Streaming Status
```
GET /sessions/{session_id}/streaming-status
```

**Response**:
```json
{
  "session_id": "session-uuid",
  "has_active_worker": true,
  "streaming_clients": ["client-1", "client-2"],
  "queue_size": 5,
  "is_complete": false,
  "status": "running",
  "worker_id": "worker-uuid",
  "created_at": "2024-01-01T12:00:00.000Z"
}
```

## Usage Examples

### JavaScript (Browser)
```javascript
const eventSource = new EventSource('/sessions/session-id/stream?client_id=my-client');

eventSource.onmessage = function(event) {
    const data = JSON.parse(event.data);
    
    if (data._type === 'completion') {
        console.log('Stream completed:', data.status);
        eventSource.close();
        return;
    }
    
    switch (data.type) {
        case 'info':
            console.log('Info:', data.message);
            break;
        case 'data':
            console.log('Data:', data.message);
            break;
        case 'error':
            console.log('Error:', data.message);
            break;
    }
};

eventSource.onerror = function(event) {
    console.error('SSE Error:', event);
};
```

### Python (with aiohttp)
```python
import aiohttp
import json

async with aiohttp.ClientSession() as session:
    async with session.get('/sessions/session-id/stream') as response:
        async for line in response.content:
            line = line.decode('utf-8').strip()
            if line.startswith('data: '):
                data = json.loads(line[6:])
                print(f"Received: {data}")
```

### cURL
```bash
curl -N -H "Accept: text/event-stream" \
  "http://localhost:8000/sessions/session-id/stream?client_id=curl-client"
```

## Test Tools

### HTML Test Page
Open `test_sse.html` in a browser to:
- Connect to any session's SSE stream
- View real-time events with formatting
- Monitor connection metrics
- Get session streaming status

### Python Test Client
Run the test client:
```bash
python test_sse_client.py <session_id> [base_url]
```

Example:
```bash
python test_sse_client.py abc-123-def http://localhost:8000
```

## Features

### Real-time Processing Updates
- Reaction parsing progress
- SMILES/SMARTS/SMIRKS generation updates
- Validation results
- Error handling

### Multiple Client Support
- Multiple clients can connect to the same session stream
- Each client gets the same events
- Client tracking for monitoring

### Graceful Fallback
- If no active worker, streams historical content from database
- Handles worker completion and cleanup
- Error recovery

### Production Ready
- Proper SSE formatting
- CORS headers for browser compatibility
- Connection management
- Resource cleanup

## Event Flow

1. **Session Worker Starts**
   ```
   run_service_agent() yields results
   ↓
   Results saved to database (existing functionality)
   ↓
   Results added to streaming queue
   ↓
   Queued results available for SSE streaming
   ```

2. **Client Connects**
   ```
   GET /sessions/{id}/stream
   ↓
   Find active worker for session
   ↓
   Stream from worker's queue
   ↓
   Send completion signal when done
   ```

3. **Historical Content**
   ```
   No active worker found
   ↓
   Get existing content from database
   ↓
   Stream historical content
   ↓
   Send completion signal
   ```

## Configuration

### Headers
The SSE endpoint sets appropriate headers:
```
Content-Type: text/event-stream
Cache-Control: no-cache
Connection: keep-alive
X-Accel-Buffering: no  # Disable nginx buffering
```

### Queue Management
- Uses `asyncio.Queue` for thread-safe operations
- Automatic cleanup when workers complete
- Timeout handling for client connections

## Error Handling

- Network errors: Client automatically reconnects (browser default)
- Worker errors: Sent as completion events with error status
- Queue errors: Graceful degradation with error messages
- JSON parsing errors: Handled in client code

This implementation provides a robust, real-time streaming solution for session content while maintaining compatibility with existing database storage functionality.
