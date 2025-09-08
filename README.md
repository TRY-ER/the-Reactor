# Reactor - AI-Powered Chemical Reaction Analysis

![Chemical Reactor](https://img.shields.io/badge/Chemical-Reactor-blue?style=for-the-badge)
![Python](https://img.shields.io/badge/Python-3.12+-brightgreen?style=for-the-badge&logo=python)
![React](https://img.shields.io/badge/React-19.1+-blue?style=for-the-badge&logo=react)
![FastAPI](https://img.shields.io/badge/FastAPI-Latest-teal?style=for-the-badge&logo=fastapi)

A comprehensive AI-powered platform for analyzing chemical reactions, extracting molecular representations (SMILES, SMARTS, SMIRKS), and providing interactive web-based visualization of chemical structures.

## 🧪 Features

- **AI-Powered Analysis**: Specialized AI agents for chemical reaction parsing and molecular representation
- **Multiple Molecular Formats**: Generate SMILES, SMARTS, and SMIRKS representations
- **Interactive Web Interface**: React-based frontend for project and session management
- **Real-time Processing**: Asynchronous analysis with progress tracking
- **Chemical Structure Visualization**: Built-in chemical renderer component
- **Project Management**: Organize analysis work into projects and sessions
- **RESTful API**: Comprehensive FastAPI backend with full CRUD operations
- **Database Persistence**: SQLite database for storing analysis results

## 🏗️ Architecture

The project consists of several integrated components:

```
reactor/
├── frontend/          # React TypeScript application
├── backend/           # Python FastAPI services
│   ├── api_endpoint/  # Main API server
│   └── app/          # Core application logic
├── oss-server/        # Open-source LLM server integration
├── analysis/          # Jupyter notebooks and data analysis
└── docs/             # Documentation
```

## 🚀 Quick Start

### Prerequisites

- **Python 3.12+**
- **Node.js 16+**
- **npm or yarn**

### 1. Backend Setup

#### API Server
```bash
# Navigate to backend API directory
cd backend/api_endpoint

# Install Python dependencies
pip install -r requirements.txt

# Set up environment variables
cp .env.example .env
# Edit .env with your configuration

# Initialize database
python migrate_db.py

# Start the API server
uvicorn main:app --reload --host 0.0.0.0 --port 8000
```

The API server will be available at: `http://localhost:8000`
- API Documentation: `http://localhost:8000/docs`
- Interactive API: `http://localhost:8000/redoc`

#### Core Backend Services
```bash
# Navigate to backend directory
cd backend

# Install dependencies
pip install -r requirements.txt

# Run tests
python run_tests.py
```

### 2. Frontend Setup

```bash
# Navigate to frontend directory
cd frontend

# Install dependencies
npm install

# Start development server
npm start
```

The frontend will be available at: `http://localhost:3000`

### 3. OSS Server (Optional)

```bash
# Navigate to oss-server directory
cd oss-server

# Install dependencies (if different from backend)
pip install -r requirements.txt

# Configure API keys for external services
# Edit configuration files as needed

# Run server-specific tests
python test_krutrim_api.py
python test_openai_api.py
```

## 🔧 Configuration

### Environment Variables

Create `.env` files in the respective directories:

**Backend API (.env)**
```env
# Database
DATABASE_URL=sqlite:///./projects.db

# CORS Settings
ALLOWED_ORIGINS=http://localhost:3000,http://localhost:3001

# API Configuration
API_KEY=your_api_key_here
DEBUG=True

# External Services
GROQ_API_KEY=your_groq_api_key
OPENAI_API_KEY=your_openai_api_key
```

**Frontend (.env)**
```env
REACT_APP_API_BASE_URL=http://localhost:8000
REACT_APP_ENVIRONMENT=development
```

## 📖 Usage

### 1. Creating a Project

```bash
# Using curl
curl -X POST "http://localhost:8000/projects/" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Organic Synthesis Project",
    "description": "Analysis of organic synthesis reactions",
    "sessions": []
  }'
```

### 2. Creating a Session

```bash
curl -X POST "http://localhost:8000/sessions/" \
  -H "Content-Type: application/json" \
  -d '{
    "name": "Reaction Analysis Session",
    "project_id": "your-project-id"
  }'
```

### 3. Web Interface Workflow

1. **Dashboard**: View all projects and sessions
2. **Project Details**: Manage sessions within a project
3. **Session Details**: Set queries, start analysis workers, monitor progress
4. **Chemical Renderer**: Visualize SMILES and SMARTS structures

### 4. AI Analysis Agents

The system includes specialized AI agents:

- **Parser Agent**: Extracts individual reactions from text
- **SMILES Agent**: Generates SMILES representations for molecules  
- **Reaction SMILES Agent**: Creates reaction SMILES strings
- **SMARTS Agent**: Generates SMARTS patterns for reactions
- **SMIRKS Agent**: Creates SMIRKS transformation patterns

## 🧪 Example Chemical Reactions

The system can analyze various reaction types:

```
# Synthesis
8 Fe + S₈ → 8 FeS

# Decomposition  
2 H₂O₂ → 2 H₂O + O₂

# Single Replacement
Zn + CuSO₄ → ZnSO₄ + Cu

# Double Replacement
AgNO₃ + NaCl → AgCl + NaNO₃
```

## 🛠️ Development

### Running Tests

**Backend Tests**
```bash
cd backend
python run_tests.py

# Or with pytest
pytest
```

**API Tests**
```bash
cd backend/api_endpoint
pytest tests/
```

**Frontend Tests**
```bash
cd frontend
npm test
```

### Code Structure

**Backend API Architecture**
- `main.py`: FastAPI application setup
- `router.py`: API route definitions
- `models.py`: Database models
- `db.py`: Database configuration
- `session_runner.py`: Asynchronous analysis workers

**Frontend Architecture**
- `src/components/`: React components including ChemicalRenderer
- `src/services/`: API service calls
- `src/components/Dashboard/`: Project management interface

## 📊 Database Schema

The system uses SQLite with the following main entities:

- **Projects**: Organize analysis work
- **Sessions**: Individual analysis instances within projects
- **Analysis Results**: Store molecular representations and metadata

## 🤝 API Reference

### Projects API
- `GET /projects/`: List all projects
- `POST /projects/`: Create new project
- `GET /projects/{id}`: Get project by ID
- `PUT /projects/{id}`: Update project
- `DELETE /projects/{id}`: Delete project

### Sessions API
- `GET /sessions/`: List all sessions
- `POST /sessions/`: Create new session
- `GET /sessions/{id}`: Get session by ID
- `PUT /sessions/{id}`: Update session
- `DELETE /sessions/{id}`: Delete session

### Analysis Workers
- Session state management: `init` → `active` → `completed`
- Real-time progress tracking
- Asynchronous processing with worker assignment

## 🔬 Chemical Data

The system includes comprehensive chemical reaction data:
- 67+ different reaction types
- Detailed mechanism descriptions
- Example reactions with balanced equations
- Reference sources and documentation

## 🐛 Troubleshooting

### Common Issues

1. **Database Connection Issues**
   ```bash
   # Reset database
   cd backend/api_endpoint
   rm projects.db
   python migrate_db.py
   ```

2. **CORS Errors**
   - Check `ALLOWED_ORIGINS` in backend `.env`
   - Ensure frontend URL is included

3. **API Key Issues**
   - Verify API keys in environment variables
   - Check external service connectivity

4. **Port Conflicts**
   - Backend: Default port 8000
   - Frontend: Default port 3000
   - Modify ports in startup commands if needed

## 🤝 Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🙏 Acknowledgments

- Chemical reaction data sourced from various educational and research institutions
- Built with FastAPI, React, and modern web technologies
- Utilizes RDKit for chemical informatics
- Integrates with various AI/LLM services for analysis

## 📞 Support

For issues, questions, or contributions:
- Create an issue on GitHub
- Check the API documentation at `/docs` endpoint
- Review the comprehensive test suites for examples

---

**Chemical Reactor** - Transforming chemical analysis through AI-powered molecular intelligence.
