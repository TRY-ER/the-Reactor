import sys
import os
from pathlib import Path

# Add the benchmarking directory to Python path
benchmarking_dir = Path(__file__).parent
sys.path.insert(0, str(benchmarking_dir))

# Also add the reactor root directory
reactor_dir = benchmarking_dir.parent
sys.path.insert(0, str(reactor_dir))

# Ensure the backend directory is also accessible
backend_dir = reactor_dir / "backend"
if backend_dir.exists():
    sys.path.insert(0, str(backend_dir))
