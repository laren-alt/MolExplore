import rdkit.Chem as Chem
from rdkit.Chem import AllChem, rdMolDescriptors, DataStructs
from sympy import symbols, Eq, solve
import re
from collections import Counter

def calculate_properties(structure):
    mol = Chem.MolFromSmiles(structure)
    molecular_weight = rdMolDescriptors.CalcExactMolWt(mol)
    return {"molecular_weight": molecular_weight}

def find_similar_molecules(target_molecule, molecule_collection, threshold=0.7):
    # Convert the target molecule to an RDKit molecule object
    target_mol = Chem.MolFromSmiles(target_molecule)
    if target_mol is None:
        raise ValueError(f"Invalid SMILES string for target molecule: {target_molecule}")

    # Generate a fingerprint for the target molecule
    target_fp = AllChem.GetMorganFingerprintAsBitVect(target_mol, radius=2, nBits=2048)

    similar_molecules = []

    for molecule in molecule_collection:
        # Get the SMILES string of the molecule in the collection
        db_smiles = molecule.get("structure")
        db_name = molecule.get("name", "Unknown")
        db_formula = molecule.get("formula")
        if not db_smiles:
            continue

        # Convert to RDKit molecule and compute fingerprint
        db_mol = Chem.MolFromSmiles(db_smiles)
        if db_mol is None:
            continue
        db_fp = AllChem.GetMorganFingerprintAsBitVect(db_mol, radius=2, nBits=2048)

        # Calculate Tanimoto similarity
        similarity = DataStructs.TanimotoSimilarity(target_fp, db_fp)
        if similarity >= threshold:
            similar_molecules.append({"_id": str(molecule["_id"]), "name": db_name, "formula": db_formula, "structure": db_smiles, "similarity": round(similarity, 2)})

    # Sort by similarity in descending order
    similar_molecules = sorted(similar_molecules, key=lambda x: x["similarity"], reverse=True)

    return similar_molecules

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
    """Convert small molecule names (like H2, O2) into valid RDKit objects."""
    # Regular expression to detect diatomic molecules (e.g., H2, O2, N2)
    diatomic_pattern = re.compile(r"([A-Za-z]+)(\d*)")
    
    # Handle simple diatomic molecules (like H2, O2, N2)
    match = diatomic_pattern.match(compound)
    if match:
        element = match.group(1)
        count = match.group(2)
        count = int(count) if count else 2  # Default to 2 if no count is specified
        # Generate SMILES for diatomic molecules (e.g., H2 -> [H][H], O2 -> [O][O])
        smiles = f"[{element}][{element}]"
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError(f"Failed to create molecule from SMILES: {smiles}")
        return mol

    # For non-diatomic molecules, return the original compound
    mol = Chem.MolFromSmiles(compound)
    if mol is None:
        raise ValueError(f"Failed to create molecule from SMILES: {compound}")
    return mol

def parse_formula(formula):
    """Parse chemical formula into element counts."""
    element_pattern = re.compile(r"([A-Z][a-z]*)(\d*)")
    parsed_elements = Counter()
    for match in element_pattern.finditer(formula):
        element = match.group(1)
        count = int(match.group(2)) if match.group(2) else 1
        parsed_elements[element] += count
    return parsed_elements