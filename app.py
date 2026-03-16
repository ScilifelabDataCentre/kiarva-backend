# Setup script for Flask app

from flask import Flask
from flask_smorest import Api
from flask_migrate import Migrate
from flask_cors import CORS
from sqlalchemy import inspect
from sqlalchemy.exc import SQLAlchemyError

from db import db

from resources.immunediscoverdata import blp as ImmuneDiscoverDataBlueprint
from loaders import *

def create_app(config_class=None):
    app = Flask(__name__)
    # CORS(app, origins=[os.getenv("FRONTEND_URL")])
    CORS(app, origins=["*"])

    app.config.from_object(config_class or 'config.Config')

    db.init_app(app)
    if not app.config.get("TESTING"):
        migrate = Migrate(app, db)
    api = Api(app)

    api.register_blueprint(ImmuneDiscoverDataBlueprint)

    # Load tsv data to db as well as pre-load plots if in prod
    with app.app_context():
        try:
            inspector = inspect(db.engine)
            if "immunediscoverdata" in inspector.get_table_names():
                load_tsv_to_db()
                if not app.debug and not app.config.get("TESTING"):
                    load_plot_data_to_dict()
            else:
                print("Table 'immunediscoverdata' not found (if first-time setup, retry after DB upgrade).")
        except SQLAlchemyError as e:
            print("SQLAlchemy error during data preload:", e)

    return app