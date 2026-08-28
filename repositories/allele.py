# Scripts relating to fetching allele data from db

from db import db
from models.immunediscoverdata import ImmuneDiscoverDataModel
from repositories.filters import allele_column, plot_selection_criteria
from utils.regex import plot_options_regex

# Whether an allele exists and is one the plots cover. Asked of the data rather than of the
# dictionaries pre-calculated at startup, because those are only populated in prod - reading
# plottability off them made the same request 404 in prod and return an all-zero plot under
# debug and pytest, which left the 404 impossible to cover with a test.
def is_plottable_allele(allele_name, plot_type):
    column = allele_column(plot_type)

    return db.session.query(
        ImmuneDiscoverDataModel.query.with_entities(ImmuneDiscoverDataModel.id
            ).filter(column == allele_name).filter(*plot_selection_criteria()).exists()
    ).scalar()

def get_allele_sequence(allele_name):
    allele_data = ImmuneDiscoverDataModel.query.with_entities(
    ImmuneDiscoverDataModel.db_name,
    ImmuneDiscoverDataModel.sequence
    ).where(ImmuneDiscoverDataModel.db_name == allele_name).distinct().all()

    if len(allele_data) < 1:
        return {}
    else:
        return {'allele': allele_name, 'sequence': allele_data[0][1]}
    
# There is a mismatch between actual allele names and their plot options
# in some cases. We therefore translate from the plot options to the
# actual name before fetching data.
# plot_selection_criteria() is what makes taking data[0] correct rather than
# arbitrary: restricted to the plottable rows, a (gene, allele) pair maps to
# exactly one db_name. Without the '_F' exclusion it maps to up to six, since
# a flanking variant shares its parent's 'allele' value, and which one came
# back first was left to the query plan.
def get_db_name_from_options(gene, allele):
    data = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.db_name,
        ImmuneDiscoverDataModel.gene,
        ImmuneDiscoverDataModel.allele
        ).where(ImmuneDiscoverDataModel.gene.regexp_match(plot_options_regex(gene)), ImmuneDiscoverDataModel.allele == allele).filter(*plot_selection_criteria()).distinct().all()

    if data:
        return data[0][0]
    else:
        return "Not found"
        