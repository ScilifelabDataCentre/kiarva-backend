# In the frontend, allele selection is split into 4 parts:
# Loci, Gene type, Gene, Allele
# This script fetches the next selectable part based on currently 
# selected options
# Example:
# Existing data in db - [IGHV1-2*01, IGHV1-3*02]
# Function input - [IGHV], output - [1-2, 1-3]
# Function input - [IGHV1-2*], output [01]

from models.immunediscoverdata import ImmuneDiscoverDataModel
from repositories.filters import plot_selection_criteria
from utils.regex import plot_options_regex

def get_plot_options(selection):
    data_out = []

    if not selection:
        return data_out
    
    # if '*' exists in function input, that means we are expecting an
    # 'allele' output, which means the last part of the selection
    select_allele = ('*' in selection)
    
    if (not select_allele):
        data = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.gene
            # startswith(autoescape=True) rather than like(selection + '%'): the selection
            # comes from a request, and LIKE reads '_' as a single-character wildcard. A
            # selection of "_" matched every gene. '%' is already rejected by the schema,
            # but '_' cannot be - it is a legitimate character in allele names.
            ).distinct().filter(ImmuneDiscoverDataModel.gene.startswith(selection, autoescape=True)).filter(*plot_selection_criteria()).all()
    else:
        data = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.allele,
            ImmuneDiscoverDataModel.gene
            # selection[:-1] to omit the '*' from query
            # i.e. input "IGHV1-2*"" -> match "IGHV1-2"
            ).distinct().filter(ImmuneDiscoverDataModel.gene.regexp_match(plot_options_regex(selection[:-1]))).filter(*plot_selection_criteria()).all()

    try:
        for row in data:
            if (not select_allele):
                if "," in row[0]:
                    next_selection = [item[len(selection):] for item in row[0].split(",")]
                else:
                    next_selection = [row[0][len(selection):]]
            else:
                next_selection = [row[0]]
            for item in next_selection:
                data_out.append(item)
    except IndexError as e:
        print(e)

    data_out = list(set(data_out))
    data_out.sort()
    return data_out