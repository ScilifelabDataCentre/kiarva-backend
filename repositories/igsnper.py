# Scripts relating to fetching igSNPer data from db
# for more info on igSNPer check https://kiarva.scilifelab.se/methodology

from models.immunediscoverdata import ImmuneDiscoverDataModel

def get_igSNPer_data(allele_name):
    igSNPer_data = ImmuneDiscoverDataModel.query.with_entities(
    ImmuneDiscoverDataModel.IgSNPer_uncommon,
    ImmuneDiscoverDataModel.IgSNPer_SNPs,
    ImmuneDiscoverDataModel.db_name
    ).where(ImmuneDiscoverDataModel.db_name == allele_name).distinct().all()

    # Result from db should be on the form [(score, SNPs, allele_name)].
    # If current allele has no associated igSNPer data, the result should be [(None, None, allele_name)]
    # As far as I understand, it should not be possible to have a score without SNPs and vice versa,
    # but we check here if both are None before sending an "empty" response.
    if igSNPer_data[0][0] is None and igSNPer_data[0][1] is None:
        return {'igSNPer_score': None, 'igSNPer_SNPs': []}
    
    # Should contain floating point value
    igSNPer_score = igSNPer_data[0][0]
    
    # igSNPer_data[0][1] should be a string of SNPs separated by semicolons, which we split by to
    # respond with a list of SNPs.
    if len(igSNPer_data[0][1]) > 0:
        igSNPer_SNPs = igSNPer_data[0][1].split(';')
        # remove empty strings from list
        igSNPer_SNPs = [x for x in igSNPer_SNPs if x.strip()]
    else:
        igSNPer_SNPs = []

    return {'igSNPer_score': igSNPer_score, 'igSNPer_SNPs': igSNPer_SNPs}