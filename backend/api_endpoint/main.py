
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from dotenv import load_dotenv
import os

# Load environment variables from .env file
load_dotenv()

# Import db module to ensure tables are created
import db

app = FastAPI(
    title="Chemical Reactor API",
    description="""
    ## Chemical Reaction Analysis API
    
    This API provides comprehensive tools for analyzing chemical reactions using AI agents.
    
    ### Key Features
    
    * **Project Management**: Organize your chemical analysis work into projects
    * **Session Management**: Create sessions within projects for individual reaction analyses
    * **AI-Powered Analysis**: Automatically extract SMILES, SMARTS, and SMIRKS representations
    * **Asynchronous Processing**: Run analysis in background with real-time progress tracking
    * **Query Management**: Store, update, and rerun analyses with different chemical reactions
    * **State Persistence**: All progress is saved to database, resilient to disconnections
    
    ### Workflow
    
    1. **Create Project** → Organize your work
    2. **Create Session** → Individual analysis workspace  
    3. **Set Query** → Chemical reaction text to analyze
    4. **Start Worker** → Begin AI analysis in background
    5. **Monitor Progress** → Track analysis status and results
    6. **Retrieve Results** → Get SMILES/SMARTS/SMIRKS data
    
    ### AI Agents
    
    The system uses specialized AI agents for:
    - **Parser Agent**: Extracts individual reactions from text
    - **SMILES Agent**: Generates SMILES representations for molecules
    - **Reaction SMILES Agent**: Creates reaction SMILES strings
    - **SMARTS Agent**: Generates SMARTS patterns for reactions
    - **SMIRKS Agent**: Creates SMIRKS transformation patterns
    
    ### Example Chemical Reactions
    
    - Synthesis: `8 Fe + S₈ → 8 FeS`
    - Decomposition: `2 H₂O₂ → 2 H₂O + O₂`
    - Single Replacement: `Zn + CuSO₄ → ZnSO₄ + Cu`
    - Double Replacement: `AgNO₃ + NaCl → AgCl + NaNO₃`
    """,
    version="1.0.0",
    contact={
        "name": "Chemical Reactor API Support",
        "email": "support@example.com",
    },
    license_info={
        "name": "MIT",
    },
)

# Get allowed origins from environment variable, default to localhost:3000
origins = os.getenv("ALLOWED_ORIGINS", "http://localhost:3000").split(",")

# Add CORS middleware
app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# Import and include the router
from router import router as main_router
app.include_router(main_router, prefix="")


# Run the app with uvicorn if executed directly
if __name__ == "__main__":
    import uvicorn
    uvicorn.run("main:app", host="0.0.0.0", port=8000, reload=True)

@app.get("/")
def read_root():
    return {"Hello": "World"}
