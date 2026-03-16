# Scripts for calculating frequencies of alleles and amino acids in different
# populations

from collections import Counter

from flask import current_app
from models.immunediscoverdata import ImmuneDiscoverDataModel
from repositories.aminoacid import get_aminoacid_allele_list, get_aminoacid_top_allele
from repositories.population import get_populations

from constants import allele_superpopulation_frequencies, allele_population_frequencies, aminoacid_allele_superpopulation_frequencies, aminoacid_allele_population_frequencies
from utils.regex import plot_options_regex

# Hedestam group requested specific order that follows their 
# research paper, sort in this order
superpopulation_order = [
    "AFR",
    "EUR",
    "EAS",
    "SAS",
    "AMR",
    "ALL"
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
    'PUR',
    'ALL'
]

# calculate the frequency that an allele or aminoacid appears in a population, alt a superpopulation
def calculate_frequencies(allele_name, population_type, plot_type):
    if plot_type == "genomic":
        db_name = "db_name"
    elif plot_type == "aminoacid":
        db_name = "db_name_AA"

    if population_type == "superpopulation":
        pop_order = superpopulation_order
    elif population_type == "population":
        pop_order = subpopulation_order

    cases = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.case,
        getattr(ImmuneDiscoverDataModel, population_type)
        ).distinct().all()
    cases_with_allele = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.case,
        getattr(ImmuneDiscoverDataModel, population_type),
        getattr(ImmuneDiscoverDataModel, db_name)
        ).where(getattr(ImmuneDiscoverDataModel, db_name) == allele_name).distinct().all()

    populations = [col[1] for col in cases]
    populations_with_allele = [col[1] for col in cases_with_allele]
    pop_count = Counter(populations)
    pop_with_allele_count = Counter(populations_with_allele)

    pop_count["ALL"] = sum(pop_count.values())
    pop_with_allele_count["ALL"] = sum(pop_with_allele_count.values())

    data_out = []
    for pop_placement in pop_order:
        for pop in pop_count:
            if pop_placement == pop:
                data_out.append({
                    'population': pop,
                    'n': pop_with_allele_count[pop],
                    'frequency': round(pop_with_allele_count[pop]/pop_count[pop], 5)
                })
    
    return data_out

# create a .tsv formated table with frequency data for the requested allele/gene and type (genomic or amino acid),
# which can then be downloaded by a user.
def create_frequencies_table(allele_or_gene, plot_type, full_gene = False):

    alleles = []

    # if full gene is requested, create a query to fetch 
    # all alleles of that gene
    if full_gene:
        if plot_type == "genomic":
            allele_data = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.db_name,
            ImmuneDiscoverDataModel.gene,
            ).filter(ImmuneDiscoverDataModel.gene.regexp_match(plot_options_regex(allele_or_gene))).filter(ImmuneDiscoverDataModel.db_name.notlike('%_F%')).distinct().all()

            for item in allele_data:
                alleles.append({'allele': item[0]})
    
        elif plot_type == "aminoacid":
            allele_data = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.db_name_AA,
            ImmuneDiscoverDataModel.db_name_AA_list,
            ImmuneDiscoverDataModel.gene,
            ).filter(ImmuneDiscoverDataModel.gene.regexp_match(plot_options_regex(allele_or_gene))).filter(ImmuneDiscoverDataModel.db_name_AA != None).distinct(ImmuneDiscoverDataModel.db_name_AA).all()

            for item in allele_data:
                # split on second item to get aa_list in form ['aa1','aa2] instead of 'aa1,aa2'
                alleles.append({'allele': item[0], 'aa_list': item[1].split(",")})

    # if single allele, just use requested allele name
    else:
        alleles = [{'allele': allele_or_gene}]

    plot_data_all_alleles = []

    # fetch population names and set up a dictionary that gives us the desired
    # superpopulation, provided a population as key
    populations = get_populations()
    pop_dict = {}
    superpops = ["AFR",
                "EUR",
                "EAS",
                "SAS",
                "AMR"]
    for pop in populations + [{'population': pop, 'superpopulation': pop} for pop in superpops]:
        pop_dict[pop['population']] = pop['superpopulation']
    
    # loop through all alleles. If single allele it's a list of length 1, if full gene
    # a list of length >= 1
    for allele_data in alleles:
        allele_name = allele_data['allele']
        if plot_type == "aminoacid":
            if full_gene:
                aa_list = allele_data['aa_list']
            else:
                # if single allele, we did not fetch db_name_AA and db_name_AA_list in the query above,
                # fetch them using our repository functions
                allele_name = get_aminoacid_top_allele(allele_name)['allele_aa']
                aa_list = get_aminoacid_allele_list(allele_name)['aa_allele_list']

        plot_data_subpops = {}
        plot_data_superpops = {}
        # If running on prod we have population frequencies pre-calculated in dictionaries,
        # use them directly
        if not current_app.debug and not current_app.config.get("TESTING"):
            if plot_type == "genomic":
                plot_data_superpops = allele_superpopulation_frequencies[allele_name]
                plot_data_subpops = allele_population_frequencies[allele_name]
            elif plot_type == "aminoacid":
                plot_data_superpops = aminoacid_allele_superpopulation_frequencies[allele_name]
                plot_data_subpops = aminoacid_allele_population_frequencies[allele_name]
        else:
            plot_data_superpops = calculate_frequencies(allele_name, "superpopulation", plot_type)
            plot_data_subpops = calculate_frequencies(allele_name, "population", plot_type)

        # assuming for now that we do not keep "ALL". If it later turns out that we need to show the aggregated
        # "ALL" data, we need to rename them to show for each one if it's referring to aggregated 
        # subpop or superpop data.
        plot_data_superpops = [item for item in plot_data_superpops if item['population'] != "ALL"]
        plot_data_subpops = [item for item in plot_data_subpops if item['population'] != "ALL"]
        populations = [item for item in populations if item['population'] != "ALL"]

        plot_data_combined = []

        # Add the plot data together and add the new dict keys 'allele' and 'superpopulation'.
        # For superpopulation plot data, 'population'=='superpopulation'.
        # For amino acids, add new key 'collapsed_translated_sequence' which contains db_name_AA_list.
        for item in plot_data_superpops + plot_data_subpops:
            item['allele'] = allele_name
            if plot_type == "aminoacid":
                item['collapsed_translated_sequence'] = aa_list
            item['superpopulation'] = pop_dict[item['population']]
            
            plot_data_combined.append(item)

        plot_data_all_alleles += plot_data_combined

    # Use a hard coded column order to put the keys in the desired order of the data of the final tsv file
    col_order = []
    if plot_type == "aminoacid":
        col_order += ['collapsed_translated_sequence']
    col_order += ['allele',
                'population',
                'superpopulation',
                'frequency',
                'n']

    # convert plot_data_all_alleles into a .tsv formated string
    plot_data_tsv_string = '\t'.join(col_order)
    for item in plot_data_all_alleles:
        item_ordered = {}
        for col in col_order:
            if col == 'collapsed_translated_sequence' and plot_type == 'genomic':
                continue
            item_ordered[col] = item[col]
        plot_data_tsv_string += '\n' + '\t'.join([str(item_ordered[k]) for k in item_ordered])
    
    return plot_data_tsv_string


    