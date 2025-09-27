from rdkit import Chem
from rdkit.Chem import rdChemReactions, Draw
import base64
from io import BytesIO


class Viewer:
    def __init__(self, input_type):
        ALLOWED_INPUT_TYPES = ["SMILES", "SMARTS"] 
        if input_type not in ALLOWED_INPUT_TYPES:
            raise ValueError(f"Invalid input_type '{input_type}'. Must be one of {ALLOWED_INPUT_TYPES}.")
        self.input_type = input_type

    def render(self, data):
        if self.input_type == "SMILES":
            return self.encode_image_to_base64(self._render_smiles(data))

        elif self.input_type == "SMARTS":
            return self.encode_image_to_base64(self._render_smarts(data))

    def encode_image_to_base64(self, image):
        buffered = BytesIO()
        image.save(buffered, format="PNG")
        img_str = base64.b64encode(buffered.getvalue()).decode("utf-8")
        return img_str

    def _render_smiles(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return "Invalid SMILES"
        img = Draw.MolToImage(mol)
        return img

    def _render_smarts(self, smarts):
        try:
            rxn = rdChemReactions.ReactionFromSmarts(smarts)
            if rxn is None:
                return "Invalid SMARTS"
            img = Draw.ReactionToImage(rxn)
            return img
        except Exception as e:
            return f"Error processing SMARTS: {e}"

if __name__ == "__main__":
    def show_image(img_str):
        from PIL import Image
        import io
        img_data = base64.b64decode(img_str)
        image = Image.open(io.BytesIO(img_data))
        image.show()

    # Example usage
    viewer_smiles = Viewer("SMILES")
    img_smiles = viewer_smiles.render("CCO")  # Ethanol
    print(f"SMILES Image (base64): {img_smiles[:50]}...")  # Print first 50 chars
    show_image(img_smiles)


    viewer_smarts = Viewer("SMARTS")
    img_smarts = viewer_smarts.render("[C:1]=[O:2]>>[C:1][OH:2]")  # Simple reaction SMARTS
    print(f"SMARTS Image (base64): {img_smarts[:50]}...")  # Print first 50 chars
    show_image(img_smarts)