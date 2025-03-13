# MolExplore

MolExplore is an **interactive web application** for exploring and managing a database of chemical molecules. It provides features for **searching**, **viewing**, **analyzing**, and **commenting** on molecules, as well as an **admin panel** for managing the data. Additionally, the project includes a **Quiz** for practicing chemical identification and a **Reaction Simulator** for predicting and balancing chemical equations.

## Table of Contents
1. [Key Features](#key-features)
2. [Technologies Used](#technologies-used)
3. [Project Structure](#project-structure)
4. [Setup and Installation](#setup-and-installation)
5. [Usage](#usage)

## Key Features

1. **Home Page & Gallery**  
   - Displays a list of molecules (from MongoDB) with options to search by name, filter by weight range, and sort results.
   - Each molecule has a thumbnail preview of its 2D structure.

2. **Molecule Details**  
   - Shows comprehensive data for a selected molecule: formula, molecular weight, bond count, and atom types.
   - Renders 2D and 3D views of the molecule (using RDKit and 3Dmol.js).
   - Lists similar molecules based on fingerprint similarity.

3. **Community Comments**  
   - Users can add comments/notes to each molecule and “like” these comments.
   - Only the admin can delete comments.

4. **Quiz (Multiple-Choice)**  
   - Randomly displays a molecule’s structure.
   - Users select the correct name from multiple-choice options.
   - Provides immediate feedback on correctness.

5. **Reaction Simulator**
   - **Predict Product**: Suggests a possible product based on common chemical rules and oxidation states.
   - **Balance Reaction**: Balances the equation by finding stoichiometric coefficients (uses matrix-based math from SymPy).

6. **Statistics Page**
   - Displays charts (pie, bar, scatter) of overall molecule data (element distribution, molecular weight distribution, etc.).

7. **Admin Panel**
   - Protected by a simple login system.
   - Allows CRUD (Create, Read, Update, Delete) operations on molecules.

## Technologies Used

### Back-End
- **Python**
- **Flask** (Web Framework)
- **MongoDB** (Data storage)
- **pymongo** (Python driver for MongoDB)

### Chemistry Libraries
- **RDKit**
  - For molecule parsing, fingerprint calculation, 2D structure generation.
- **SymPy**
  - For balancing chemical equations using a symbolic/matrix approach.
- **pymatgen**
  - For analyzing elements, oxidation states, and providing additional chemistry logic.

### Front-End
- **HTML & CSS** (with **Bootstrap** for responsive design)
- **JavaScript** (with **jQuery** and **AJAX** for dynamic interactions)
- **3Dmol.js** for 3D molecular visualization
- **Chart.js** / **Highcharts** (optional) for data visualization

### Other
- **Jinja2** (templating in Flask)
- **Flask-WTF** (form handling & validation)
- **Werkzeug** (security, password hashing)

## Project Structure

```
MolExplore/
├── app/
│   ├── __init__.py           # Flask app factory, initializes extensions
│   ├── routes.py             # Main routes for molecules, comments, quiz, etc.
│   ├── models.py             # MongoDB collection references
│   ├── forms.py              # Flask-WTF forms for login, molecule input, etc.
│   ├── utils.py              # Helper functions for RDKit, reaction balancing
│   ├── templates/            # Jinja2 HTML templates
│   │   ├── base.html
│   │   ├── home.html
│   │   ├── gallery.html
│   │   ├── details.html
│   │   ├── quiz.html
│   │   ├── admin.html
│   │   └── ...
│   ├── static/               # CSS, JS, images, etc.
│   │   ├── css/
│   │   ├── js/
│   │   └── images/
│   └── ...
├── config.py                 # Configuration (SECRET_KEY, MONGO_URI)
├── run.py                    # Entry point to start the Flask server
├── molecules.json            # Example data for seeding
└── README.md                 # Project overview (this file)
```

## Setup and Installation

1. **Clone the Repository**  
   ```bash
   git clone https://github.com/laren-alt/MolExplore.git
   cd MolExplore
   ```

2. **Create and Activate a Virtual Environment (Recommended)**  
   ```bash
   python -m venv venv
   source venv/bin/activate
   # or on Windows:
   # venv\Scriptsctivate
   ```

3. **Install Dependencies**  
   ```bash
   pip install -r requirements.txt
   ```
   Make sure you have **MongoDB** running locally or update `MONGO_URI` in `config.py` to point to your own server.

4. **Run the App**  
   ```bash
   python run.py
   ```
   By default, the Flask app will run at **http://127.0.0.1:5000**.

5. **Seed Admin User & Molecule Data**  
   - The project automatically seeds an **admin** user (`username: admin`, `password: admin123`) if none exists.
   - It also seeds molecules from `molecules.json` if the database is empty.

## Usage

- **Homepage**: Lists molecules (paginated), allows basic searching.
- **Gallery**: Advanced search, filter by molecular weight, sorting, etc.
- **Details Page**: For each molecule, view 2D/3D representation, comments, and a mini reaction simulator.
- **Admin Panel**: Go to `/admin` (login required) to add, edit, or delete molecules.
- **Comments**: All users can post comments or like/unlike them; only admin can delete comments.
- **Quiz**: Access `/quiz` for a multiple-choice test on random molecules.
- **Reaction Simulator**: Input reactants, predict product, and then balance the resulting chemical equation.
