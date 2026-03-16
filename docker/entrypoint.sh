#!/usr/bin/env bash
echo "Starting setup..."
flask db upgrade
echo "Done with DB upgrade. Starting Gunicorn..."
exec gunicorn --bind 0.0.0.0:5000 'app:create_app()' --access-logfile - --access-logformat '%(h)s - - [%(t)s] %(r)s %(s)s %(b)s'
