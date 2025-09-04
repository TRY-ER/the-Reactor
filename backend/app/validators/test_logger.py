import logging
import io
from rdkit import rdBase, Chem

# --- Setup logging to capture warnings ---
# Redirect RDKit's logs to the Python logging system
rdBase.LogToPythonLogger()

# Create an in-memory stream and a handler to capture logs
stream = io.StringIO()
logger = logging.getLogger("rdkit")
handler = logging.StreamHandler(stream)
formatter = logging.Formatter("%(levelname)s: %(message)s")
handler.setFormatter(formatter)
logger.addHandler(handler)
logger.setLevel(logging.WARNING)  # Set the level to capture warnings

# --- RDKit operation that causes a warning ---
# This invalid SMILES string will cause a parsing warning
bad_smiles = "[5*]NC(N[5*].[*]Nc1ccc(NCCC[4*])cc1"
print(f"Attempting to process SMILES: {bad_smiles}")

# The warning will be captured by the stream, not printed to the console
mol = Chem.MolFromSmiles(bad_smiles)

# --- Process captured warnings ---
captured_output = stream.getvalue()

if captured_output:
    print("\n--- Captured warnings found ---")
    print(captured_output)
else:
    print("\n--- No warnings were captured ---")

# --- Cleanup ---
# It's good practice to remove the handler after you're done
logger.removeHandler(handler)
