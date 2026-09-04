# Scripts relating to fetching allele data from db

from models.immunediscoverdata import ImmuneDiscoverDataModel
from repositories.filters import plot_selection_criteria
from utils.regex import plot_options_regex

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
        