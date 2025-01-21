from flask import Flask
from flask_pymongo import PyMongo
from flask_bootstrap import Bootstrap
from flask_wtf.csrf import CSRFProtect

mongo = PyMongo()
csrf = CSRFProtect()

def create_app():
    app = Flask(__name__)
    app.config.from_object("config.Config")
    
    # Initialize extensions
    mongo.init_app(app)
    Bootstrap(app)
    csrf.init_app(app)
    
    # Set up application context explicitly for any components that need it
    with app.app_context():
        # Register blueprints/routes within app context
        from .routes import main
        app.register_blueprint(main)
    
    return app
