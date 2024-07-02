from sqlalchemy import func
from sqlalchemy.orm import aliased
from models.immunediscoverdata import ImmuneDiscoverDataModel

# calculate the frequency that an allele appears in a population, alt a superpopulation
def calculate_allele_frequencies(allele_name, population_type):
        distinct_cases = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.case,
            getattr(ImmuneDiscoverDataModel, population_type),
            ).distinct().subquery()
        
        aliased_distinct_cases = aliased(ImmuneDiscoverDataModel, distinct_cases)

        pop_count = aliased_distinct_cases.query.with_entities(
            getattr(aliased_distinct_cases, population_type),
            aliased_distinct_cases.case,
            func.count(getattr(aliased_distinct_cases, population_type))
            ).group_by(getattr(aliased_distinct_cases, population_type)).all()
        
        pop_count_with_allele = aliased_distinct_cases.query.with_entities(
            getattr(aliased_distinct_cases, population_type),
            aliased_distinct_cases.case,
            ImmuneDiscoverDataModel.db_name,
            func.count(getattr(aliased_distinct_cases, population_type))
            ).join(ImmuneDiscoverDataModel, aliased_distinct_cases.case == ImmuneDiscoverDataModel.case).where(ImmuneDiscoverDataModel.db_name == allele_name).group_by(getattr(aliased_distinct_cases, population_type)).all()

        data_out = []
        for i in range(len(pop_count)):
            data_out.append({
                'population': pop_count[i][0],
                'n': pop_count_with_allele[i][3],
                'frequency': pop_count_with_allele[i][3]/pop_count[i][2]
            })
        
        return data_out