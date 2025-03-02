from flask import Blueprint, render_template, request, redirect, url_for, flash, session, abort, jsonify
from werkzeug.security import generate_password_hash, check_password_hash
from bson.objectid import ObjectId
from rdkit import Chem
import json
import os
import random
from math import ceil
from collections import Counter
from datetime import datetime

from .models import get_user_collection, get_molecule_collection, get_comment_collection
from .forms import MoleculeForm, LoginForm, ReactionForm
from .utils import (
    find_similar_molecules,
    calculate_bond_count,
    get_atom_types,
    handle_small_molecule,
    parse_formula,
    parse_reaction,       # Matrix-based parsing
    get_oxidation_state,  # Used for product prediction
    predict_product,      # Uses oxidation states
    balance_reaction      # Matrix-based balancing
)


main = Blueprint("main", __name__)

# Seed admin user (runs once during app setup)
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

# Seed molecule data from molecules.json if the database is empty
def seed_molecule_data():
    molecule_collection = get_molecule_collection()
    
    # Check if the database is empty
    if molecule_collection.count_documents({}) == 0:
        json_file_path = "molecules.json"

        # Ensure file exists
        if os.path.exists(json_file_path):
            with open(json_file_path, "r") as file:
                molecules = json.load(file)
            
            for molecule in molecules:
                molecule["createdAt"] = datetime.utcnow()
                molecule["updatedAt"] = datetime.utcnow()
                molecule_collection.insert_one(molecule)

            print(f"Seeded {len(molecules)} molecules from molecules.json into the database.")
        else:
            print("molecules.json not found. No data seeded.")

# Call seeding functions
seed_admin_user()
seed_molecule_data()

@main.route("/")
def home():
    # Get the current page number and search query from the URL parameters
    page = request.args.get('page', 1, type=int)  # Get the current page number from the query string
    search_query = request.args.get('search', '', type=str)  # Get the search query from the query string
    per_page = 12  # Define how many molecules per page

    molecule_collection = get_molecule_collection()

    # Construct the search filter
    search_filter = {}
    if search_query:
        search_filter = {
            "$or": [
                {"name": {"$regex": search_query, "$options": "i"}},
                {"formula": {"$regex": search_query, "$options": "i"}},
                {"molecular_weight": {"$regex": search_query, "$options": "i"}}
            ]
        }

    # Fetch paginated molecules based on search filter
    molecule_data = list(molecule_collection.find(search_filter, {"_id": 0})
                         .skip((page - 1) * per_page).limit(per_page))

    # Get the total number of filtered molecules to calculate the number of pages
    total_molecules = molecule_collection.count_documents(search_filter)
    total_pages = ceil(total_molecules / per_page)

    # Generate 2D molecule images for each structure
    for molecule in molecule_data:
        smiles = molecule.get("structure", "")
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            molecule["mol_block"] = Chem.MolToMolBlock(mol)

    return render_template(
        "home.html", 
        molecule_data=molecule_data, 
        total_pages=total_pages, 
        current_page=page,
        search_query=search_query  # Pass the search query back to the template
    )


@main.route("/gallery")
def gallery():
    """Displays a paginated and filterable gallery of molecules."""
    query = {}
    search = request.args.get("search", "")
    sort_by = request.args.get("sort_by", "name")
    min_weight = request.args.get("min_weight", type=str, default=None)
    max_weight = request.args.get("max_weight", type=str, default=None)
    page = request.args.get("page", type=int, default=1)
    per_page = 8

    # Convert min_weight and max_weight safely
    try:
        min_weight = float(min_weight) if min_weight and min_weight.strip() else None
    except ValueError:
        min_weight = None
    
    try:
        max_weight = float(max_weight) if max_weight and max_weight.strip() else None
    except ValueError:
        max_weight = None

    # Validate: min_weight should not be greater than max_weight
    if min_weight is not None and max_weight is not None and min_weight > max_weight:
        flash("Minimum weight cannot be greater than maximum weight.", "danger")
        return redirect(url_for("main.gallery"))

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
    form = ReactionForm()

    try:
        molecule_collection = get_molecule_collection()
        molecule = molecule_collection.find_one({"_id": ObjectId(molecule_id)})

        if not molecule:
            flash("Molecule not found.", "danger")
            return redirect(url_for("main.gallery"))

        # Calculate any additional info you need
        molecule["bond_count"] = calculate_bond_count(molecule["structure"])
        molecule["atom_types"] = get_atom_types(molecule["structure"])

        # Find similar molecules
        all_molecules = list(molecule_collection.find({}, {"_id": 1, "name": 1, "structure": 1, "formula": 1}))
        similar_molecules = find_similar_molecules(molecule["structure"], all_molecules, str(molecule["_id"]))

        # Fetch recent comments
        comments = list(
            get_comment_collection()
            .find({"molecule_id": ObjectId(molecule_id)})
            .sort("timestamp", -1)
            .limit(10)
        )
        for comment in comments:
            comment["_id"] = str(comment["_id"])

        return render_template(
            "details.html",
            molecule=molecule,
            similar_molecules=similar_molecules,
            form=form,
            comments=comments,
            user_role=session.get("role")
        )

    except Exception as e:
        print("Error in details route:", e)
        flash(f"An error occurred: {e}", "danger")
        return redirect(url_for("main.gallery"))


@main.route('/comment/unfavorite/<comment_id>', methods=['POST'])
def unfavorite_comment(comment_id):
    comment = get_comment_collection().find_one({"_id": ObjectId(comment_id)})
    if not comment:
        return jsonify({"success": False, "message": "Comment not found"}), 404

    new_favorite_count = max(0, comment.get("favorite_count", 0) - 1)
    get_comment_collection().update_one(
        {"_id": ObjectId(comment_id)},
        {"$set": {"favorite_count": new_favorite_count}}
    )

    return jsonify({"success": True, "new_favorite_count": new_favorite_count})


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
        similar_molecules = find_similar_molecules(molecule["structure"], molecule_collection, str(molecule["_id"]))

        
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
        molecule['molecular_weight'] = float(molecule['molecular_weight'])
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
            "molecular_weight": float(form.molecular_weight.data),
            "structure": form.structure.data,
            "description": form.description.data,
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
            "molecular_weight": float(form.molecular_weight.data),
            "structure": form.structure.data,
            "description": form.description.data,
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
            session["role"] = user["role"]
            session["user_id"] = str(user["_id"])
            flash("Login successful!", "success")
            return redirect(url_for("main.admin"))
        flash("Invalid username or password.", "danger")
    return render_template("login.html", form=form)

@main.route("/logout")
def logout():
    session.pop("user", None)
    flash("Logged out successfully.", "success")
    return redirect(url_for("main.home"))

@main.route("/details/<molecule_id>/add_comment", methods=["POST"])
def add_comment(molecule_id):
    comment_text = request.json.get("comment_text")
    commenter_name = request.json.get("commenter_name", "Anonymous").strip()
    if not commenter_name:
        return jsonify({"success": False, "message": "Commenter Name cannot be empty."}), 400
    
    if not comment_text:
        return jsonify({"success": False, "message": "Comment cannot be empty."}), 400

    comment_data = {
        "molecule_id": ObjectId(molecule_id),
        "user_id": session.get("user_id") if session.get("role") == "admin" else None,
        "text": comment_text,
        "commenter_name": commenter_name,
        "timestamp": datetime.utcnow(),
    }
    comment_collection = get_comment_collection()
    inserted_comment = comment_collection.insert_one(comment_data)

    # Convert MongoDB ObjectId to string for JSON response
    comment_data["_id"] = str(inserted_comment.inserted_id)

    return jsonify({
        "success": True,
        "message": "Comment added successfully!",
        "comment": {
            "id": comment_data["_id"],
            "molecule_id": str(comment_data["molecule_id"]),
            "user_id": str(comment_data["user_id"]) if comment_data["user_id"] else None,
            "text": comment_data["text"],
            "timestamp": comment_data["timestamp"].strftime("%Y-%m-%d %H:%M:%S"),
        }
    })


@main.route("/comment/delete/<comment_id>", methods=["POST"])
def delete_comment(comment_id):
    """Delete a comment (admin or owner only) via AJAX."""
    if not session.get("user"):
        return jsonify({"success": False, "message": "You need to log in to delete a comment."}), 403

    comment = get_comment_collection().find_one({"_id": ObjectId(comment_id)})
    if not comment:
        return jsonify({"success": False, "message": "Comment not found."}), 404

    if session.get("role") == "admin" or session.get("user_id") == comment["user_id"]:
        get_comment_collection().delete_one({"_id": ObjectId(comment_id)})
        return jsonify({"success": True, "message": "Comment deleted successfully!"})
    else:
        return jsonify({"success": False, "message": "You don't have permission to delete this comment."}), 403
    
@main.route('/comment/favorite/<comment_id>', methods=['POST'])
def favorite_comment(comment_id):
    comment = get_comment_collection().find_one({"_id": ObjectId(comment_id)})

    if not comment:
        return jsonify({"success": False, "message": "Comment not found"}), 404

    # Increase favorite count by 1
    new_favorite_count = comment.get("favorite_count", 0) + 1

    get_comment_collection().update_one(
        {"_id": ObjectId(comment_id)},
        {"$set": {"favorite_count": new_favorite_count}}
    )

    return jsonify({"success": True, "new_favorite_count": new_favorite_count})


@main.route("/quiz")
def chemistry_quiz():
    """Render the quiz page with a random molecule."""
    molecule = get_molecule_collection().aggregate([{"$sample": {"size": 1}}]).next()
    # Determine quiz type (Multiple Choice or Text Input)
    quiz_type = random.choice(["multiple_choice", "text"])
    print(quiz_type, "----")
    # Generate answer choices for multiple choice
    options = []
    if quiz_type == "multiple_choice":
        options = [molecule["name"]]
        other_molecules = list(get_molecule_collection().aggregate([{"$sample": {"size": 3}}]))
        options += [m["name"] for m in other_molecules]
        random.shuffle(options)

    return render_template("quiz.html", molecule=molecule, quiz_type=quiz_type, options=options)

@main.route("/quiz/check-answer", methods=["POST"])
def check_quiz_answer():
    """Validate the user's answer."""
    data = request.json
    molecule_id = data.get("molecule_id")
    user_answer = data.get("user_answer", "").strip().lower()

    molecule = get_molecule_collection().find_one({"_id": ObjectId(molecule_id)})
    if not molecule:
        return jsonify({"success": False, "message": "Molecule not found."})

    correct_answer = molecule["name"].strip().lower()
    is_correct = user_answer == correct_answer

    return jsonify({
        "correct": is_correct,
        "message": "Correct! 🎉" if is_correct else f"Incorrect. The correct answer is {molecule['name']}."
    })

@main.route('/predict', methods=['POST'])
def predict():
    """API endpoint to predict the product of a chemical reaction."""
    data = request.get_json()
    if 'reactants' not in data:
        return jsonify({"error": "Missing reactants"}), 400

    try:
        reactants = parse_reaction(data['reactants'])
        product = predict_product(reactants)

        if product == "Unknown Product":
            return jsonify({"error": "Could not predict the reaction."}), 400

        return jsonify({"reactants": reactants, "predicted_product": product})

    except ValueError as e:
        return jsonify({"error": str(e)}), 400


@main.route('/balance', methods=['POST'])
def balance():
    """API endpoint to balance a chemical reaction."""
    data = request.get_json()
    if 'reactants' not in data or 'product' not in data:
        return jsonify({"error": "Missing reactants or product"}), 400

    try:
        reactants = parse_reaction(data['reactants'])
        product = data['product']
        balanced_equation = balance_reaction(reactants, product)

        return jsonify({"reactants": reactants, "product": product, "balanced_equation": balanced_equation})

    except ValueError as e:
        return jsonify({"error": str(e)}), 400

