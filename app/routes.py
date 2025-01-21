from flask import Blueprint, render_template, request, redirect, url_for, flash, session, abort
from werkzeug.security import generate_password_hash, check_password_hash
from bson.objectid import ObjectId
from .models import get_user_collection, get_molecule_collection
from .forms import MoleculeForm, LoginForm, ReactionForm
from .utils import find_similar_molecules, calculate_bond_count, get_atom_types, handle_small_molecule, parse_formula
from rdkit import Chem
import json
import re
from collections import Counter
from datetime import datetime
from sympy import symbols, Eq, solve

main = Blueprint("main", __name__)

# Seed admin user (can be called manually or on app startup)
def seed_admin_user():
    user_collection = get_user_collection()
    if not user_collection.find_one({"username": "admin"}):
        current_time = datetime.utcnow()
        admin_user = {
            "username": "admin",
            "password": generate_password_hash("admin123"),  # Default password
            "role": "admin",
            "createdAt": current_time,
            "updatedAt": current_time,
        }
        user_collection.insert_one(admin_user)
        print("Admin user created with username: 'admin' and password: 'admin123'")

# Uncomment to run seed function only once during the app setup
seed_admin_user()

@main.route("/")
def home():
    # Load JSON data
    json_file_path = "molecules.json"
    with open(json_file_path, "r") as file:
        molecule_data = json.load(file)

    # Generate 2D molecule images for each structure
    for molecule in molecule_data:
        smiles = molecule["structure"]
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            molecule["mol_block"] = Chem.MolToMolBlock(mol) 

    return render_template("home.html", molecule_data=molecule_data)

@main.route("/gallery")
def gallery():
    """Displays a paginated and filterable gallery of molecules."""
    query = {}
    search = request.args.get("search", "")
    sort_by = request.args.get("sort_by", "name")
    min_weight = request.args.get("min_weight", type=float, default=None)
    max_weight = request.args.get("max_weight", type=float, default=None)
    page = request.args.get("page", type=int, default=1)
    per_page = 8

    # Search by name
    if search:
        query["name"] = {"$regex": search, "$options": "i"}

    # Filter by molecular weight range
    if min_weight is not None or max_weight is not None:
        query["molecular_weight"] = {}
        if min_weight is not None:
            query["molecular_weight"]["$gte"] = min_weight
        if max_weight is not None:
            query["molecular_weight"]["$lte"] = max_weight

    # Get the molecule collection
    molecule_collection = get_molecule_collection()

    # Total number of molecules for pagination
    total_molecules = molecule_collection.count_documents(query)
    total_pages = (total_molecules + per_page - 1) // per_page

    # Fetch molecules with sorting and pagination
    molecules_cursor = (
        molecule_collection.find(query)
        .sort(sort_by, 1)
        .skip((page - 1) * per_page)
        .limit(per_page)
    )

    # Convert the cursor to a list
    molecules = list(molecules_cursor)
    for molecule in molecules:
        molecule["_id"] = str(molecule["_id"])

    return render_template(
        "gallery.html",
        molecules=molecules,
        search=search,
        sort_by=sort_by,
        min_weight=min_weight,
        max_weight=max_weight,
        page=page,
        total_pages=total_pages,
    )


@main.route("/details/<molecule_id>", methods=["GET", "POST"])
def details(molecule_id):
    form = ReactionForm()  # Initialize the form
    reaction_result = None  # Initialize reaction result

    try:
        # Get the molecule collection
        molecule_collection = get_molecule_collection()
        
        # Fetch the target molecule
        molecule = molecule_collection.find_one({"_id": ObjectId(molecule_id)})
        if not molecule:
            flash("Molecule not found.", "danger")
            return redirect(url_for("main.gallery"))

        # Add bond count and atom types
        molecule["bond_count"] = calculate_bond_count(molecule["structure"])
        molecule["atom_types"] = get_atom_types(molecule["structure"])

        # Fetch all molecules for similarity search
        all_molecules = list(molecule_collection.find({}, {"_id": 1,"name": 1, "structure": 1}))
        similar_molecules = find_similar_molecules(molecule["structure"], all_molecules)

        # Handle chemical reaction simulation if form is submitted
        if request.method == "POST" and form.validate_on_submit():
            reaction_input = form.reaction.data.strip()  # Remove extra spaces

            if not reaction_input:
                flash("Please enter a valid reaction.", "warning")
                return redirect(url_for("main.details", molecule_id=molecule_id))

            # Ensure the reaction input has the right format (reactants -> products)
            if "->" not in reaction_input:
                flash("Invalid reaction format. Make sure the reaction has '->'.", "warning")
                return redirect(url_for("main.details", molecule_id=molecule_id))

            # Parse the reaction input
            try:
                reactants, products = reaction_input.split("->")
                reactants = [x.strip() for x in reactants.split("+")]
                products = [x.strip() for x in products.split("+")]
            except ValueError:
                flash("Invalid reaction format. Make sure to separate reactants and products with '->'.", "warning")
                return redirect(url_for("main.details", molecule_id=molecule_id))

            # Handle small molecule names (e.g., H2 -> H2)
            reactants = [handle_small_molecule(r) for r in reactants]
            products = [handle_small_molecule(p) for p in products]

            # Create variables for each compound
            try:
                variables = {compound: symbols(compound, integer=True) for compound in reactants + products}
            except Exception as e:
                flash(f"Error in creating symbols: {str(e)}", "danger")
                return redirect(url_for("main.details", molecule_id=molecule_id))

            # Debug: Print variables to ensure correct symbol creation
            print("Variables:", variables)

            # Create equations based on conservation of elements
            element_counts = {}
            for compound in variables:
                mol = handle_small_molecule(compound)  # Use the handle_small_molecule function
                if mol is None:
                    flash(f"Invalid molecule in reaction: {compound}", "danger")
                    return redirect(url_for("main.details", molecule_id=molecule_id))

                # Count elements in the molecule
                for atom in mol.GetAtoms():
                    element = atom.GetSymbol()
                    count = atom.GetTotalNumHs() + atom.GetNumImplicitHs() + 1  # Includes implicit Hs
                    element_counts.setdefault(element, []).append((variables[compound], count))

            # Build equations
            equations = []
            for element, compounds in element_counts.items():
                equations.append(Eq(sum(c * coeff for c, coeff in compounds[:len(reactants)]),
                                    sum(c * coeff for c, coeff in compounds[len(reactants):])))
            print(equations, "--------------")
            # Solve the system of equations
            solution = solve(equations, list(variables.values()), dict=True)
            if not solution:
                flash("Unable to balance the reaction.", "danger")
                return redirect(url_for("main.details", molecule_id=molecule_id))

            # Generate balanced reaction string
            balanced_reactants = " + ".join(f"{solution[var]} {var.name}" for var in variables.values()[:len(reactants)])
            balanced_products = " + ".join(f"{solution[var]} {var.name}" for var in variables.values()[len(reactants):])
            reaction_result = f"{balanced_reactants} -> {balanced_products}"

        return render_template(
            "details.html",
            molecule=molecule,
            similar_molecules=similar_molecules,
            reaction_result=reaction_result,
            form=form  # Pass the form to the template
        )
    except Exception as e:
        flash(f"An error occurred: {e}", "danger")
        return redirect(url_for("main.gallery"))
    
@main.route("/admin/molecule/details/<molecule_id>")
def molecule_details(molecule_id):
    try:
        # Fetch the molecule from the database
        molecule = get_molecule_collection().find_one({"_id": ObjectId(molecule_id)})
        if not molecule:
            flash("Molecule not found.", "danger")
            return redirect(url_for("main.admin"))

        # Validate and fetch the structure (SMILES) for similarity search
        structure = molecule.get("structure")
        if not structure:
            flash("Molecule structure data is missing.", "warning")
            return render_template("molecule_details.html", molecule=molecule, similar_molecules=[])

        # Fetch all molecules for similarity search
        molecule_collection = list(get_molecule_collection().find({}, {"name": 1, "formula": 1, "structure": 1}))
        similar_molecules = find_similar_molecules(structure, molecule_collection)
        return render_template(
            "molecule_details.html",
            molecule=molecule,
            similar_molecules=similar_molecules
        )
    except Exception as e:
        flash("An error occurred while processing the molecule details.", "danger")
        return redirect(url_for("main.admin"))


@main.route("/stats")
def stats():
    molecule_collection = get_molecule_collection()
    molecules = list(molecule_collection.find({}, {"name": 1, "formula": 1, "molecular_weight": 1, "structure": 1}))

    # Calculate element distribution
    element_distribution = Counter()

    for molecule in molecules:
        molecule["_id"] = str(molecule["_id"])
        formula = molecule.get("formula", "")
        parsed_elements = parse_formula(formula)
        molecule["parsed_formula"] = dict(parsed_elements)  # Add parsed formula to each molecule
        element_distribution.update(parsed_elements)

    # Format element_distribution for pie chart
    element_distribution_3d = [[element, int(count)] for element, count in element_distribution.items()]

    # Molecular weights for histogram
    molecular_weights = [molecule.get("molecular_weight", 0) for molecule in molecules]
    weight_bins = {f"{i}-{i+50}": 0 for i in range(0, int(max(molecular_weights, default=0)) + 50, 50)}
    for weight in molecular_weights:
        for bin_range in weight_bins:
            low, high = map(int, bin_range.split("-"))
            if low <= weight < high:
                weight_bins[bin_range] += 1
                break

    # Scatter data: molecular weight vs bond count
    scatter_data = [
        {"x": molecule.get("molecular_weight", 0), "y": calculate_bond_count(molecule.get("structure", ""))}
        for molecule in molecules
    ]

    # Find the most common element
    most_common_element = max(element_distribution.items(), key=lambda x: x[1])[0] if element_distribution else None

    return render_template(
        "stats.html",
        molecules=molecules,
        molecular_weights=weight_bins,
        element_distribution=dict(element_distribution),
        element_distribution_3d=element_distribution_3d,
        scatter_data=scatter_data,
        most_common_element=most_common_element,
    )



@main.route("/molecule/<molecule_name>")
def molecule_detail(molecule_name):
    # Load JSON data
    json_file_path = "molecules.json"
    with open(json_file_path, "r") as file:
        molecule_data = json.load(file)

    # Find the molecule by name
    molecule = next((m for m in molecule_data if m["name"] == molecule_name), None)

    if not molecule:
        # Return a 404 error if the molecule is not found
        abort(404)

    return render_template("explore_detail.html", molecule=molecule)

@main.route("/admin")
def admin():
    """Admin Panel: Displays the list of molecules with pagination."""
    if not session.get("user"):
        flash("You need to log in to access the admin panel.", "danger")
        return redirect(url_for("main.login"))

    molecule_collection = get_molecule_collection()

    # Pagination parameters
    per_page = 8
    page = int(request.args.get("page", 1))
    total_molecules = molecule_collection.count_documents({})
    total_pages = (total_molecules + per_page - 1) // per_page

    # Fetch paginated molecules
    molecules = list(
        molecule_collection.find()
        .skip((page - 1) * per_page)
        .limit(per_page)
    )
    for molecule in molecules:
        molecule["_id"] = str(molecule["_id"])

    return render_template(
        "admin.html",
        molecules=molecules,
        page=page,
        total_pages=total_pages,
    )


@main.route("/admin/create_molecule", methods=["GET", "POST"])
def create_molecule():
    """Add a new molecule."""
    if not session.get("user"):
        flash("You need to log in to add a molecule.", "danger")
        return redirect(url_for("main.login"))
    
    form = MoleculeForm()
    molecule_collection = get_molecule_collection()
    
    if form.validate_on_submit():
        data = {
            "name": form.name.data,
            "formula": form.formula.data,
            "molecular_weight": form.molecular_weight.data,
            "structure": form.structure.data,
        }
        molecule_collection.insert_one(data)
        flash("Molecule added successfully!", "success")
        return redirect(url_for("main.admin"))
    
    return render_template("create_molecule.html", form=form)

@main.route("/admin/update_molecule/<molecule_id>", methods=["GET", "POST"])
def update_molecule(molecule_id):
    """Edit an existing molecule."""
    if not session.get("user"):
        flash("You need to log in to update a molecule.", "danger")
        return redirect(url_for("main.login"))
    
    molecule_collection = get_molecule_collection()
    molecule = molecule_collection.find_one({"_id": ObjectId(molecule_id)})
    
    if not molecule:
        flash("Molecule not found!", "danger")
        return redirect(url_for("main.admin"))
    
    form = MoleculeForm(data=molecule)
    
    if form.validate_on_submit():
        updated_data = {
            "name": form.name.data,
            "formula": form.formula.data,
            "molecular_weight": form.molecular_weight.data,
            "structure": form.structure.data,
        }
        molecule_collection.update_one({"_id": ObjectId(molecule_id)}, {"$set": updated_data})
        flash("Molecule updated successfully!", "success")
        return redirect(url_for("main.admin"))
    
    return render_template("update_molecule.html", form=form, molecule=molecule)

@main.route("/admin/delete/<molecule_id>")
def delete_molecule(molecule_id):
    if not session.get("user"):
        flash("You need to log in to perform this action.", "danger")
        return redirect(url_for("main.login"))

    try:
        molecule_collection = get_molecule_collection()
        molecule_collection.delete_one({"_id": ObjectId(molecule_id)})
        flash("Molecule deleted successfully!", "success")
    except Exception as e:
        flash("Failed to delete molecule.", "danger")
    return redirect(url_for("main.admin"))

@main.route("/login", methods=["GET", "POST"])
def login():
    form = LoginForm()
    if form.validate_on_submit():
        user_collection = get_user_collection()
        user = user_collection.find_one({"username": form.username.data})
        if user and check_password_hash(user["password"], form.password.data):
            session["user"] = user["username"]
            flash("Login successful!", "success")
            return redirect(url_for("main.admin"))
        flash("Invalid username or password.", "danger")
    return render_template("login.html", form=form)

@main.route("/logout")
def logout():
    session.pop("user", None)
    flash("Logged out successfully.", "success")
    return redirect(url_for("main.home"))
