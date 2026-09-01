#!/usr/bin/env bash
# set -e so a failed migration stops the container rather than being followed by the server.
# Without it, 'flask db upgrade' could fail and the next line still exec'd gunicorn, which
# then passed its readiness check while every data endpoint raised "no such table" as a 500.
set -e

echo "Starting setup..."
# SKIP_DATA_LOAD tells create_app() that this process is the migration and not the server.
# The flask CLI has to build an app through that factory to run a migration, and at that
# point the tables the migration is about to create do not exist yet - so the load has to be
# skipped here, and only here. Scoped to this one command, so the server below still loads
# the TSVs and validates them, and a missing table there is the error it should be.
SKIP_DATA_LOAD=1 flask db upgrade
echo "Done with DB upgrade. Starting Gunicorn..."
exec gunicorn --bind 0.0.0.0:5000 'app:create_app()' --access-logfile - --access-logformat '%(h)s - - [%(t)s] %(r)s %(s)s %(b)s'
