from flask_pymongo import PyMongo
from flask import current_app

def get_molecule_collection():
    """
    Returns the MongoDB molecule collection.
    Assumes Flask app is configured with MongoDB URI.
    """
    mongo = PyMongo(current_app)
    return mongo.db.molecules
