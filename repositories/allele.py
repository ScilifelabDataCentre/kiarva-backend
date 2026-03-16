# Scripts relating to fetching allele data from db

from models.immunediscoverdata import ImmuneDiscoverDataModel
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
    # actual name before fetching data
def get_db_name_from_options(gene, allele):
    data = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.db_name,
        ImmuneDiscoverDataModel.gene,
        ImmuneDiscoverDataModel.allele
        ).where(ImmuneDiscoverDataModel.gene.regexp_match(plot_options_regex(gene)), ImmuneDiscoverDataModel.allele == allele).distinct().all()

    if data:
        return data[0][0]
    else:
        return "Not found"
        