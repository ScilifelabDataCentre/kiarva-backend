# Script that fetches any matching alleles' sequence to an input gene sequence
# segment.
# Example:
# Data in db - 
# [allele: allele1, sequence: EXAMPLESEQ]
# [allele: allele2, sequence: SEQEXAMPLE]
# [allele: allele3, sequence: DIFFERENTSEQ]
# Input: 
# "EXAMPLE"
# Output: 
# [allele: allele1, sequence: EXAMPLESEQ, position: 0]
# [allele: allele2, sequence: SEQEXAMPLE, position: 3]

import re

from models.immunediscoverdata import ImmuneDiscoverDataModel

def sequence_search(sequence_str):
    sequence_str = sequence_str.upper()
    # autoescape=True because the search term is user input: LIKE treats '%' and '_' as
    # wildcards, so without it a search for a single '%' matched - and returned - every
    # sequence in the database.
    sequence_data = ImmuneDiscoverDataModel.query.with_entities(
    ImmuneDiscoverDataModel.db_name,
    ImmuneDiscoverDataModel.sequence
    ).filter(ImmuneDiscoverDataModel.sequence.contains(sequence_str, autoescape=True)).filter(~ImmuneDiscoverDataModel.db_name.contains('_F', autoescape=True)).distinct().all()

    data_out = []
    if len(sequence_data) < 1:
        data_out = [{'allele': '', 'sequence': '', 'positions': []}]
        return data_out
    else:
        for row in sequence_data:
            # re.escape because the term is searched for literally. Unescaped, its
            # metacharacters would be interpreted as a pattern, and an unbalanced one
            # would raise.
            positions = [m.start() for m in re.finditer(re.escape(sequence_str), row[1])]
            data_out.append({'allele': row[0], 'sequence': row[1], 'positions': positions})

    return sorted(data_out, key=lambda d: d['allele'])