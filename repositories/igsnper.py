# Scripts relating to fetching igSNPer data from db
# for more info on igSNPer check https://kiarva.scilifelab.se/methodology

from models.immunediscoverdata import ImmuneDiscoverDataModel

def get_igSNPer_data(allele_name):
    igSNPer_data = ImmuneDiscoverDataModel.query.with_entities(
    ImmuneDiscoverDataModel.IgSNPer_uncommon,
    ImmuneDiscoverDataModel.IgSNPer_SNPs,
    ImmuneDiscoverDataModel.db_name
    ).where(ImmuneDiscoverDataModel.db_name == allele_name).distinct().all()

    # An allele that is not in the data at all yields no rows, which is different from an
    # allele that is present but has no IgSNPer columns (handled below). Returning {} here
    # follows the other repository functions and lets the resource answer with a 404; the
    # index below used to raise IndexError and surface as a 500.
    if not igSNPer_data:
        return {}

    # Result from db should be on the form [(score, SNPs, allele_name)].
    # If current allele has no associated igSNPer data, the result should be [(None, None, allele_name)]
    if igSNPer_data[0][0] is None and igSNPer_data[0][1] is None:
        return {'igSNPer_score': None, 'igSNPer_SNPs': []}

    # Should contain floating point value
    igSNPer_score = igSNPer_data[0][0]

    # igSNPer_data[0][1] should be a string of SNPs separated by semicolons, which we split by to
    # respond with a list of SNPs.
    #
    # A score with no SNPs is normal rather than a contradiction: the score counts uncommon
    # SNPs, so a score of 0.0 means there were none to list and the column is null. That is
    # the majority shape in the data - 209,867 rows against 216,496 with both set - and no
    # row anywhere has a score above 0.0 with no SNPs. Testing the value rather than its
    # length covers the null and the empty string in one, and stops 29 of the 732 plottable
    # alleles answering 500.
    if igSNPer_data[0][1]:
        igSNPer_SNPs = igSNPer_data[0][1].split(';')
        # remove empty strings from list
        igSNPer_SNPs = [x for x in igSNPer_SNPs if x.strip()]
    else:
        igSNPer_SNPs = []

    return {'igSNPer_score': igSNPer_score, 'igSNPer_SNPs': igSNPer_SNPs}