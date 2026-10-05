import re

# We want to match against a particular string called "string1",
# we want the match to be exact and only found once.
# The string in which we search may look like
# "string1" or "string1,string2,string1-2,..." or "string1-2,string1,...",
# so we need an expression that matches ONLY "string1" in all of these situations.
# To do this we use ',{0,1}string1(,|$)':
# First match the character ',' 0 or 1 times, then match the desired searched after string
# then match either the character ',' or end of line (nothing).

# re.escape rather than escaping a couple of characters by hand: the name comes from a
# request, so any metacharacter left unescaped changes what the pattern matches, and an
# unbalanced one - a bare '(' or '[' - made re.search raise inside the user-defined
# function SQLite calls for REGEXP, surfacing as a 500. SQLAlchemy implements
# regexp_match on SQLite with Python's re.search, so re.escape is the matching escape.
def plot_options_regex(string_name):
    return ',{0,1}' + re.escape(string_name) + '(,|$)'