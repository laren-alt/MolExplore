from flask import Blueprint, render_template, request, redirect, url_for, flash
from .models import get_molecule_collection
from .forms import MoleculeForm
from .utils import calculate_properties, find_similar_molecules
from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import os

main = Blueprint("main", __name__)

@main.route("/")
def home():
    # Generate example molecules
    smiles_list = [
        "CCO", "CCN", "CCC", "CCCl", "CCBr", "CCCC", "CCOCC", "CCOCCC", "CCOCCCC"
    ]
    labels = ["alcohol", "amine", "alkane", "alkyl halide", "alkyl halide", "alkane", "ether", "ether", "ether"]

    # Ensure the static/images directory exists
    image_dir = os.path.join("app", "static", "images", "chemical")
    os.makedirs(image_dir, exist_ok=True)

    # Compute molecular properties
    data_points = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(smi)
        if mol:
            fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=128)
            data_points.append(fp)

    # Generate molecule images and save them
    molecule_images = []
    for i, smi in enumerate(smiles_list):
        mol = Chem.MolFromSmiles(smi)
        if mol:
            img = Draw.MolToImage(mol, size=(150, 150))  # Generate molecule image
            img_path = os.path.join(image_dir, f"mol_{i}.png")
            img.save(img_path)  # Save the image using PIL
            print(img_path.replace("app", ""))
            image_dir1 = os.path.join("static", "images", "chemical", img_path.split("\\")[-1])
            molecule_images.append((smi, image_dir1))

    return render_template("home.html", molecule_images=molecule_images)

@main.route("/gallery")
def gallery():
    molecules = get_molecule_collection().find()
    return render_template("gallery.html", molecules=molecules)

@main.route("/details/<molecule_id>")
def details(molecule_id):
    molecule = get_molecule_collection().find_one({"_id": molecule_id})
    similar_molecules = find_similar_molecules(molecule)
    return render_template("details.html", molecule=molecule, similar_molecules=similar_molecules)

@main.route("/admin", methods=["GET", "POST"])
def admin():
    form = MoleculeForm()
    if form.validate_on_submit():
        data = form.data
        data["properties"] = calculate_properties(data["formula"])
        get_molecule_collection().insert_one(data)
        flash("Molecule added successfully!", "success")
        return redirect(url_for("main.admin"))
    return render_template("admin.html", form=form)
