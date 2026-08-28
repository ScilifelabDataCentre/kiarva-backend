# Setup script for Flask app
#
# The service is stateless by design. There is no PVC: DATABASE_URL defaults to a SQLite
# file inside the container, and every pod rebuilds it from the TSVs in data/ on startup.
# Nothing a request does writes to it - the app only ever reads, transforms and serves - so
# a pod can be restarted at any time to get a clean database, and replacing the source data
# is a matter of shipping a new image rather than migrating anything.
#
# Two consequences that are easy to misread as bugs:
#   - Recording loaded files in loaded_from_tsv only skips work within one container's
#     lifetime. A restart starts from an empty database and re-reads everything.
#   - 'flask db upgrade' in docker/entrypoint.sh runs against a database that has no tables
#     yet, which is why the loader below has a branch for that and does not treat it as an
#     error. It is needed on every fresh SQLite file and implies nothing about persistence.

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
            load_plot_data_to_dict()
        else:
            # Not an error: the table does not exist yet on a first-time setup, before
            # 'flask db upgrade' has run.
            print("Table 'immunediscoverdata' not found (if first-time setup, retry after DB upgrade).")

    return app