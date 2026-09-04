# Scripts relating to fetching amino acid data from db

from models.immunediscoverdata import ImmuneDiscoverDataModel
from utils.regex import plot_options_regex

def get_aminoacid_sequence(aa_allele_name):
    allele_data = ImmuneDiscoverDataModel.query.with_entities(
    ImmuneDiscoverDataModel.db_name_AA,
    ImmuneDiscoverDataModel.sequence_AA
    ).where(ImmuneDiscoverDataModel.db_name == aa_allele_name).distinct().all()

    if len(allele_data) < 1:
        return {}
    else:
        return allele_data[0][1]
    
# Several different alleles can translate to the same amino acid, the
# amino acid lists contain the names of such alleles under one master amino acid
# allele.
# This function returns this list.
# For example: allele1, allele5, allele6 all translate to aminoacid1, then
# db_name_AA: allele1
# db_name_AA_list: [allele1, allele5, allele6]
# return: db_name_AA_list
def get_aminoacid_allele_list(aa_allele_name):
    aa_allele_data = ImmuneDiscoverDataModel.query.with_entities(
    ImmuneDiscoverDataModel.db_name_AA,
    ImmuneDiscoverDataModel.db_name_AA_list
    ).where(ImmuneDiscoverDataModel.db_name_AA == aa_allele_name).distinct().all()

    if len(aa_allele_data) == 0:
        return {'aa_allele_list': None}
    
    # Testing the value rather than its length: db_name_AA_list is nullable, and len(None)
    # raises where this function already has an answer for "no list". validate_loaded_data()
    # requires the column wherever db_name_AA is set, so this is the second line of defence
    # rather than the only one - but it is the same shape as the fix in repositories/igsnper.py,
    # and a function that already reports absence should not raise on it.
    if not aa_allele_data[0][1]:
        return {'aa_allele_list': None}

    return {'aa_allele_list': aa_allele_data[0][1].split(',')}

# Several different alleles can translate to the same amino acid, the
# amino acid lists contain the names of such alleles under one master amino acid
# allele.
# This function returns the "master" amino acid allele of the list.
# For example: allele1, allele5, allele6 all translate to aminoacid1, then
# db_name_AA: allele1
# db_name_AA_list: [allele1, allele5, allele6]
# return: db_name_AA
def get_aminoacid_top_allele(aa_allele_name):
    allele_data = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.db_name_AA,
        ImmuneDiscoverDataModel.gene,
        ImmuneDiscoverDataModel.allele,
        ImmuneDiscoverDataModel.db_name_AA_list
        ).filter(ImmuneDiscoverDataModel.db_name_AA_list.regexp_match(plot_options_regex(aa_allele_name))).distinct().all()

    if len(allele_data) < 1:
        return {}
    else:
        return {'allele': aa_allele_name, 'allele_aa': allele_data[0][0]}