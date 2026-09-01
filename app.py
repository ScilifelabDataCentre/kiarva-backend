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
#     yet, and has to build an app through this factory to do it. That step sets
#     SKIP_DATA_LOAD to say so. It is needed on every fresh SQLite file and implies nothing
#     about persistence.

import copy
import os

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

    # from_object copies the class attribute by reference, so every app built from a given
    # config would share one API_SPEC_OPTIONS dict - and apispec's to_dict() ends by
    # deep-merging the spec it just built into it. Left shared, the first app's schemas
    # accumulate there and then win for every app built after it in the process: under
    # pytest that leaks specs between tests, and under a threaded worker two concurrent
    # /openapi.json requests merge into one object. Each app gets its own copy.
    app.config["API_SPEC_OPTIONS"] = copy.deepcopy(app.config.get("API_SPEC_OPTIONS", {}))

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
    # SKIP_DATA_LOAD is how a process states that it is not the server: docker/entrypoint.sh
    # sets it for 'flask db upgrade', which cannot have the tables it is about to create. It
    # is checked before the table is looked at, because the whole point is that there is
    # nothing to look at yet.
    if os.environ.get("SKIP_DATA_LOAD"):
        print("SKIP_DATA_LOAD is set: no data loaded and nothing validated.", flush = True)
        return app

    # Anywhere else, a missing table means the migration did not run - and answering that by
    # printing a note and returning a routed app was the same silent degradation the removed
    # try/except used to cause. It is the shape that is hardest to notice, because /health
    # has no database in it: readiness passes and every data endpoint raises "no such table".
    with app.app_context():
        if "immunediscoverdata" not in inspect(db.engine).get_table_names():
            raise SourceDataError(
                "Table 'immunediscoverdata' does not exist, so there is no data to serve. "
                "Run 'SKIP_DATA_LOAD=1 flask db upgrade' first. The flag is part of the "
                "command because the flask CLI builds this same app to run a migration, and "
                "at that point the tables it is about to create do not exist yet - so "
                "without it the migration fails on this very error. docker/entrypoint.sh "
                "runs it exactly that way on every container start.")

        load_tsv_to_db()
        if not app.debug and not app.config.get("TESTING"):
            load_plot_data_to_dict()

    return app