
# We want to match against a particular string called "string1",
# we want the match to be exact and only found once.
# The string in which we search may look like
# "string1" or "string1,string2,string1-2,..." or "string1-2,string1,...",
# so we need an expression that matches ONLY "string1" in all of these situations.
# To do this we use ',{0,1}string1(,|$)':
# First match the character ',' 0 or 1 times, then match the desired searched after string
# then match either the character ',' or end of line (nothing).
def plot_options_regex(string_name):
    string_name_regex_safe = string_name.replace("*", "\\*").replace("/", "\\/")
    return ',{0,1}' + string_name_regex_safe + '(,|$)'