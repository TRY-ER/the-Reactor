import re

def find_texts_within_square_brackets(value: str) -> list[str]:
    pattern = r"\[([^\]]+)\]"
    matches = re.findall(pattern, value)
    return matches

def replace_aux_text(val: str) -> str:
    replace_chars = ["@"]
    for char in replace_chars:
        val = val.replace(char, "")
    return val

def validate_smirk_text(candidate: str) -> None:
    parts = candidate.split(">>")
    assert len(parts) == 2, "SMIRKS must contain exactly one '>>' to separate reactants and products."
    reactant = parts[0].strip()
    product = parts[1].strip()

    reactants = find_texts_within_square_brackets(reactant)
    products = find_texts_within_square_brackets(product)

    assert all(':' in r for r in reactants), "Each reactant must contain at least one atom mapping (e.g., [C:1])."
    assert all(':' in p for p in products), "Each product must contain at least one atom mapping (e.g., [C:1])."

    assert reactants, "Reactants must not be empty."
    assert products, "Products must not be empty."

    reactants = list(map(replace_aux_text, reactants))
    products = list(map(replace_aux_text, products))

    assert sorted(reactants) == sorted(products), "Reactants and products must have the same atom mappings."

if __name__ == "__main__":
    # candidate = "[C:1](=[O:2])[Cl:3].[H:99][NH:4][C:0]>> [C:1](=[O:2])[NH:4][C:0].[Cl:3][H:99]"
    candidate = "[*:1][C@:2]([*:3])([*:4])[*:5]>>[*:1][C@@:2]([*:3])([*:4])[*:5]" # chiral changes
    validate_smirk_text(candidate)