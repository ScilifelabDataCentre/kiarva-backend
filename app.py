# Setup script for Flask app

from flask import Flask
from flask_smorest import Api
from flask_migrate import Migrate
from flask_cors import CORS
from sqlalchemy import inspect

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

    # Load tsv data to db as well as pre-load plots if in prod.
    #
    # Deliberately not wrapped in a try/except. A database error here means there is no
    # data to serve, and catching it used to log the error and return a fully routed app
    # anyway - so the pod passed its readiness check and answered requests against an
    # empty database. Letting it raise means the worker never boots, the pod never
    # becomes ready, and the previous one goes on serving.
    with app.app_context():
        inspector = inspect(db.engine)
        if "immunediscoverdata" in inspector.get_table_names():
            load_tsv_to_db()
            if not app.debug and not app.config.get("TESTING"):
                load_plot_data_to_dict()
        else:
            # Not an error: the table does not exist yet on a first-time setup, before
            # 'flask db upgrade' has run.
            print("Table 'immunediscoverdata' not found (if first-time setup, retry after DB upgrade).")

    return app