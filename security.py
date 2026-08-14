import os
import re

import functools
from flask import request
from flask_smorest import abort

# Allele, gene and plot-selection names in the source data are built from a narrow
# character set - letters, digits, and * / , _ - - and the longest value in the dataset
# is 28 characters. Validating a query parameter against that before anything else
# touches it keeps untrusted text out of database queries, response bodies and
# filenames, and it is what lets the endpoints below report an error without echoing
# the request back.
NAME_PARAMETER_PATTERN = re.compile(r"\A[A-Za-z0-9*/,_-]{1,64}\Z")

def validated_name(value, param_name):
    """Return an allele/gene/selection parameter, or abort with 400 if it is unusable.

    param_name is a literal supplied by the endpoint, not request data, so it is safe to
    name in the response. The offending value deliberately is not: reflecting
    unvalidated request input back to the caller is what the security scanners flag, and
    the caller already knows what it sent.
    """
    if not value:
        abort(400, message="Missing required query parameter '" + param_name + "'.")

    if not NAME_PARAMETER_PATTERN.match(value):
        abort(400, message="Malformed value for query parameter '" + param_name + "'.")

    return value

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