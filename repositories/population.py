# Script to fetch population data (name of superpopulations and subpopulations)
# from db

from models.immunediscoverdata import ImmuneDiscoverDataModel

def get_populations():
    data_from_db = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.superpopulation,
            ImmuneDiscoverDataModel.population,
            ).distinct().all()
    
    data = []

    for entry in data_from_db:
        data.append({
            "population": entry[1],
            "superpopulation": entry[0]
        })
    # "ALL" population for an aggregated column of all population data
    data.append({
        "population": "ALL",
        "superpopulation": "ALL"
        })
    return data