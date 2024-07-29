from collections import Counter, OrderedDict
from models.immunediscoverdata import ImmuneDiscoverDataModel

# calculate the frequency that an allele appears in a population, alt a superpopulation
def calculate_allele_frequencies(allele_name, population_type):
    # Hedestam group requested specific order that follows their 
    # research paper, sort in this order
    superpopulation_order = [
        "AFR",
        "EUR",
        "EAS",
        "SAS",
        "AMR",
    ]
    subpopulation_order = [
        'ACB',
        'ASW',
        'ESN',
        'GWD',
        'LWK',
        'MSL',
        'YRI',
        'FIN',
        'GBR',
        'IBS',
        'TSI',
        'CDX',
        'CHB',
        'CHS',
        'JPT',
        'KHV',
        'BEB',
        'GIH',
        'ITU',
        'PJL',
        'STU',
        'CLM',
        'MXL',
        'PEL',
        'PUR'
    ]

    if population_type == "superpopulation":
        pop_order = superpopulation_order
    else:
        pop_order = subpopulation_order

    cases = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.case,
        getattr(ImmuneDiscoverDataModel, population_type)
        ).distinct().all()
    cases_with_allele = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.case,
        getattr(ImmuneDiscoverDataModel, population_type),
        ImmuneDiscoverDataModel.db_name
        ).where(ImmuneDiscoverDataModel.db_name == allele_name).distinct().all()

    populations = [col[1] for col in cases]
    populations_with_allele = [col[1] for col in cases_with_allele]
    pop_count = Counter(populations)
    pop_with_allele_count = Counter(populations_with_allele)

    pop_count = OrderedDict(sorted(pop_count.items()))

    data_out = []
    for pop_placement in pop_order:
        for pop in pop_count:
            if pop_placement == pop:
                data_out.append({
                    'population': pop,
                    'n': pop_with_allele_count[pop],
                    'frequency': pop_with_allele_count[pop]/pop_count[pop]
                })
    
    return data_out