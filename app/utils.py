import rdkit.Chem as Chem
from rdkit.Chem import AllChem, rdMolDescriptors, DataStructs
from sympy import symbols, Eq, solve
import re
from collections import Counter
from pymatgen.core import Composition, Element
from math import gcd


def calculate_properties(structure):
    mol = Chem.MolFromSmiles(structure)
    molecular_weight = rdMolDescriptors.CalcExactMolWt(mol)
    return {"molecular_weight": molecular_weight}

def find_similar_molecules(target_smiles, molecule_collection, target_id, tanimoto_weight=0.5, dice_weight=0.3, cosine_weight=0.2):  
    # Convert the target molecule to RDKit format
    target_mol = Chem.MolFromSmiles(target_smiles)
    if not target_mol:
        raise ValueError(f"Invalid SMILES string: {target_smiles}")

    # Compute Morgan fingerprint for the target molecule
    target_fp = AllChem.GetMorganFingerprintAsBitVect(target_mol, radius=2, nBits=2048)
    
    similar_molecules = []

    for molecule in molecule_collection:
        db_smiles = molecule.get("structure")
        db_id = str(molecule["_id"])  # Convert ObjectId to string

        # Skip if the molecule has no structure or is the same as the target
        if not db_smiles or db_id == target_id:
            continue

        # Convert database molecule to RDKit format
        db_mol = Chem.MolFromSmiles(db_smiles)
        if not db_mol:
            continue  # Skip invalid molecules

        # Compute fingerprints
        db_fp = AllChem.GetMorganFingerprintAsBitVect(db_mol, radius=2, nBits=2048)

        # Compute similarity scores
        tanimoto_sim = DataStructs.TanimotoSimilarity(target_fp, db_fp)
        dice_sim = DataStructs.DiceSimilarity(target_fp, db_fp)
        cosine_sim = DataStructs.FingerprintSimilarity(target_fp, db_fp, metric=DataStructs.CosineSimilarity)

        # Weighted similarity score
        weighted_similarity = (tanimoto_weight * tanimoto_sim) + (dice_weight * dice_sim) + (cosine_weight * cosine_sim)

        # Keep only molecules with similarity > 0.2 and < 0.95
        if 0.25 < weighted_similarity < 0.95:  
            similar_molecules.append({
                "_id": db_id,
                "name": molecule.get("name", "Unknown"),
                "formula": molecule.get("formula", "N/A"),
                "structure": db_smiles,
                "similarity": round(weighted_similarity, 3)  # Round for cleaner output
            })

    # Sort results by similarity score (highest first)
    return sorted(similar_molecules, key=lambda x: x["similarity"], reverse=True)

def calculate_bond_count(smiles: str) -> int:
    """
    Calculate the bond count of a molecule based on its SMILES structure.
    :param smiles: SMILES string of the molecule.
    :return: Total number of bonds in the molecule.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return mol.GetNumBonds()
        else:
            return 0
    except Exception as e:
        print(f"Error calculating bond count for SMILES '{smiles}': {e}")
        return 0
    
def get_atom_types(smiles: str) -> list:
    """
    Extract unique atom types from the molecule.
    :param smiles: SMILES string of the molecule.
    :return: List of unique atom types.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            atom_types = set(atom.GetSymbol() for atom in mol.GetAtoms())
            return sorted(atom_types)  # Return sorted list of atom types
        else:
            return []
    except Exception as e:
        print(f"Error extracting atom types for SMILES '{smiles}': {e}")
        return []

def simulate_reaction(reaction_input: str) -> str:
    """
    Simulate and balance a chemical reaction using SymPy.
    :param reaction_input: Reaction input string, e.g., "H2 + O2 -> H2O".
    :return: Balanced reaction as a string.
    """
    try:
        # Parse the reaction input
        reactants, products = reaction_input.split("->")
        reactants = [r.strip() for r in reactants.split("+")]
        products = [p.strip() for p in products.split("+")]

        # Create a unique symbol for each compound
        compounds = reactants + products
        coefficients = symbols(f"a:{len(compounds)}")

        # Generate equations based on elemental balance
        element_equations = []
        unique_elements = set()

        # Collect element counts for each compound
        compound_elements = {}
        for i, compound in enumerate(compounds):
            mol = Chem.MolFromSmiles(compound)
            if mol is None:
                raise ValueError(f"Invalid SMILES string: {compound}")
            element_count = {}
            for atom in mol.GetAtoms():
                element = atom.GetSymbol()
                element_count[element] = element_count.get(element, 0) + 1
                unique_elements.add(element)
            compound_elements[compound] = element_count

        # Balance each element
        for element in unique_elements:
            equation = Eq(
                sum(coefficients[i] * compound_elements[compounds[i]].get(element, 0) for i in range(len(reactants))),
                sum(coefficients[i + len(reactants)] * compound_elements[compounds[i + len(reactants)]].get(element, 0) for i in range(len(products)))
            )
            element_equations.append(equation)

        # Solve the system of equations
        solution = solve(element_equations, coefficients)

        # Find the smallest integer coefficients
        lcm = 1
        for coeff in solution.values():
            lcm = lcm * coeff.q // lcm.gcd(coeff.q)  # Least common multiple
        integer_coeffs = [int(solution.get(c, 0) * lcm) if c in solution else 0 for c in coefficients]

        # Format the balanced reaction
        balanced_reaction = " + ".join(
            f"{integer_coeffs[i]} {reactants[i]}" for i in range(len(reactants))
        ) + " -> " + " + ".join(
            f"{integer_coeffs[i + len(reactants)]} {products[i]}" for i in range(len(products))
        )
        return balanced_reaction
    except Exception as e:
        return f"Error balancing reaction: {e}"

def handle_small_molecule(compound):
    """Convert simple molecules (like H2, O2) into valid SMILES strings."""
    if not compound or not compound.strip():
        raise ValueError("Empty molecule input received.")

    diatomic_pattern = re.compile(r"([A-Za-z]+)(\d*)")

    match = diatomic_pattern.fullmatch(compound.strip())
    if match:
        element = match.group(1)
        count = int(match.group(2)) if match.group(2) else 2  # Default to 2 if unspecified
        smiles = f"[{element}]" * count  # Example: O2 -> [O][O]
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return Chem.MolToSmiles(mol)

    # Convert standard molecules using RDKit
    mol = Chem.MolFromSmiles(compound.strip())
    if mol:
        return Chem.MolToSmiles(mol)  # Ensure output is a string

    raise ValueError(f"Invalid molecule input: {compound}")

def parse_formula(formula):
    """Parse chemical formula into element counts."""
    element_pattern = re.compile(r"([A-Z][a-z]*)(\d*)")
    from collections import Counter
    parsed_elements = Counter()
    for match in element_pattern.finditer(formula):
        element = match.group(1)
        count = int(match.group(2)) if match.group(2) else 1
        parsed_elements[element] += count
    return parsed_elements

def parse_reaction(reaction_str):
    """
    Extract reactants from user input (left of '->'), e.g., "Na + Cl2 -> ?".
    Validate them using pymatgen's Composition.
    """
    left_side = reaction_str.split("->")[0].strip()
    reactants = [r.strip() for r in left_side.split("+") if r.strip()]
    for r in reactants:
        Composition(r)  # raises ValueError if invalid
    return reactants

def get_oxidation_state(element):
    """Retrieve the most common oxidation state of an element."""
    try:
        el = Element(element)
        ox_states = el.common_oxidation_states
        return ox_states[0] if ox_states else None
    except:
        return None

def predict_product(reactants):
    """
    Predict a product by combining the first two distinct elements from reactants.
    Example: "Na" + "Cl2" => "NaCl".
    """
    if not reactants:
        return "Unknown Product"

    # Gather unique element symbols from all reactants
    unique_elems = []
    for r in reactants:
        for el in Composition(r).elements:
            if el.symbol not in unique_elems:
                unique_elems.append(el.symbol)

    if len(unique_elems) < 2:
        return "Unknown Product"

    e1, e2 = unique_elems[0], unique_elems[1]
    ox1, ox2 = get_oxidation_state(e1), get_oxidation_state(e2)
    if ox1 is None or ox2 is None:
        return "Unknown Product"

    # Cross-multiply absolute oxidation numbers for simplest ratio
    g1, g2 = abs(ox2), abs(ox1)
    return f"{e1}{g1 if g1 > 1 else ''}{e2}{g2 if g2 > 1 else ''}"
    
    return product_formula

def extract_element_count(formula, element):
    """Extract element count safely from a chemical formula."""
    match = re.findall(rf'({element})(\d*)', formula)
    if match:
        count = match[0][1] or "1"  # If no number, assume 1
        return int(count)
    return 1

def balance_reaction(reactants, product):
    """
    Balance a reaction with one or more reactants and a single product
    using a matrix null-space approach.
    Example:
      reactants = ["Na", "Cl2"], product = "NaCl"
      => "2 Na + 1 Cl2 -> 2 NaCl"
    """
    from sympy import Matrix

    # Build a list of species: reactants + [product]
    species = reactants + [product]

    # Find all elements in those species
    all_text = "".join(species)
    elements = re.findall(r"[A-Z][a-z]?", all_text)
    unique_elems = sorted(set(elements))

    # Construct a matrix: rows = elements, columns = species
    # Reactants get positive counts, product is negative (so sum = 0).
    stoich_rows = []
    for elem in unique_elems:
        row = []
        for i, sp in enumerate(species):
            # Count how many times elem appears in sp
            c = parse_formula(sp).get(elem, 0)
            # For the product (last column), we use negative
            row.append(c if i < len(reactants) else -c)
        stoich_rows.append(row)

    # Convert to a sympy Matrix and compute the nullspace
    M = Matrix(stoich_rows)
    ns = M.nullspace()
    if not ns:
        return "Unable to balance"

    # Usually there's a single vector in ns for a simple reaction
    for vec in ns:
        # Scale vector to get integer coefficients
        denominators = [v.q for v in vec]
        lcm_denom = 1
        for d in denominators:
            lcm_denom = lcm_denom * d // gcd(lcm_denom, d)

        scaled = [int(lcm_denom * v) for v in vec]

        # Ensure product coefficient is positive
        if scaled[-1] < 0:
            scaled = [-x for x in scaled]

        # Skip trivial solutions if product <= 0 or all reactants = 0
        if scaled[-1] <= 0 or all(x == 0 for x in scaled[:-1]):
            continue

        # If any reactant coefficient <= 0, skip
        if any(x <= 0 for x in scaled[:-1]):
            continue

        # Build final string "2 Na + 1 Cl2 -> 2 NaCl"
        reactant_side = " + ".join(f"{scaled[i]} {r}" for i, r in enumerate(reactants))
        product_side = f"{scaled[-1]} {product}"
        return f"{reactant_side} -> {product_side}"

    return "Unable to balance"