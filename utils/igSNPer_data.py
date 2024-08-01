
from models.immunediscoverdata import ImmuneDiscoverDataModel


def get_igSNPer_data(allele_name):
    igSNPer_data = ImmuneDiscoverDataModel.query.with_entities(
    ImmuneDiscoverDataModel.IgSNPer_uncommon,
    ImmuneDiscoverDataModel.IgSNPer_SNPs,
    ImmuneDiscoverDataModel.db_name
    ).where(ImmuneDiscoverDataModel.db_name == allele_name).distinct().all()

    igSNPer_score = igSNPer_data[0][0]
    
    if len(igSNPer_data[0][1]) > 0:
        igSNPer_SNPs = igSNPer_data[0][1].split(';')
        if len(igSNPer_SNPs) > 1:
            igSNPer_SNPs = igSNPer_SNPs[:-1]
    else:
        igSNPer_SNPs = []

    return {'igSNPer_score': igSNPer_score, 'igSNPer_SNPs': igSNPer_SNPs}