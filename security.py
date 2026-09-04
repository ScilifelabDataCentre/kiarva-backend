import os

import functools
from flask import request
from flask_smorest import abort

# Query parameter validation used to live here as validated_name(). It is now declared by
# the marshmallow schemas in schemas.py, which flask_smorest applies before a view runs -
# see NAME_PATTERN there for the character set these names are held to.

def is_valid(api_key):
    return os.getenv("API_KEY") == api_key

def api_key_required(func):
    @functools.wraps(func)
    def decorator(*args, **kwargs):
        api_key = request.headers.get("X-api-key")
        if not api_key:
            abort(
                400, 
                message="Please provide an API key."
            )
        # Check if API key is correct and valid
        if is_valid(api_key):
            return func(*args, **kwargs)
        else:
            abort(
                403, 
                message="The provided API key is not valid."
            )
    return decorator