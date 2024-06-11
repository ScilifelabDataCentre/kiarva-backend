# calculate the frequency that an allele appears in a population, alt a superpopulation
from collections import Counter, OrderedDict
from models.immunediscoverdata import ImmuneDiscoverDataModel


def calculate_allele_frequencies(allele_name, calc_for_superpopulation = False):
        cases = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.case
            ).distinct().all()
        cases_with_allele = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.case,
            ImmuneDiscoverDataModel.db_name
            ).where(ImmuneDiscoverDataModel.db_name == allele_name).distinct().all()

        split_index = 2 if calc_for_superpopulation else 1

        populations = [col[0].split('_')[split_index] for col in cases]
        populations_with_allele = [col[0].split('_')[split_index] for col in cases_with_allele]
        pop_count = Counter(populations)
        pop_with_allele_count = Counter(populations_with_allele)

        pop_count = OrderedDict(sorted(pop_count.items()))

        data_out = []
        for pop in pop_count:
            data_out.append({
                'population': pop,
                'n': pop_with_allele_count[pop],
                'frequency': pop_with_allele_count[pop]/pop_count[pop]
            })
        
        return data_out