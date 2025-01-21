from flask_pymongo import PyMongo
from flask import current_app

def get_user_collection():
    """
    Returns the MongoDB users collection.
    """
    mongo = PyMongo(current_app)
    return mongo.db.users

def get_molecule_collection():
    """
    Returns the MongoDB molecules collection.
    """
    mongo = PyMongo(current_app)
    return mongo.db.molecules
