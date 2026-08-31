# Scripts for calculating frequencies of alleles and amino acids in different
# populations

from collections import Counter

from flask import current_app
from sqlalchemy import func

from db import db
from models.immunediscoverdata import ImmuneDiscoverDataModel
from repositories.aminoacid import get_aminoacid_allele_list, get_aminoacid_top_allele
from repositories.filters import allele_column, plot_selection_criteria
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

def population_display_order(population_type):
    """The population order a frequency response is reported in.

    Raises for the same reason as allele_column: chosen with an if/elif and no else, an
    unrecognised population_type left the variable unassigned and failed one line later.
    """
    if population_type == "superpopulation":
        return superpopulation_order
    if population_type == "population":
        return subpopulation_order
    raise ValueError("unknown population_type: " + repr(population_type))

# Key the per-population case totals are cached under. They are stored on the Flask app
# rather than at module level because the pytest fixtures build one app per test, each
# with its own in-memory database, and a module-level cache would let one test's totals
# be used for another's data.
POPULATION_TOTALS_KEY = "kiarva_population_totals"

def population_totals(population_type):
    """Number of distinct cases per population, plus the aggregated "ALL" total.

    This is the denominator of every frequency: how many individuals a population
    contains, which has nothing to do with which allele is being asked about. It costs a
    full scan of the table, and used to be re-run inside every calculate_frequencies()
    call - thousands of times during the startup pre-load, for one unchanging result.

    The data is read-only once loaded, so it is calculated once per app and cached.
    """
    cache = current_app.extensions.setdefault(POPULATION_TOTALS_KEY, {})

    if population_type not in cache:
        cases = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.case,
            getattr(ImmuneDiscoverDataModel, population_type)
            ).distinct().all()
        totals = Counter(case[1] for case in cases)
        totals["ALL"] = sum(totals.values())
        cache[population_type] = totals

    # Copied so that a caller cannot mutate the cached totals
    return Counter(cache[population_type])

def frequency_entries(pop_totals, cases_with_allele, pop_order):
    """Assemble the per-population entries for one allele, in the requested order.

    cases_with_allele is a Counter, so a population in which no case carries the allele
    contributes an n=0 entry rather than being left out - "we looked and found none" is
    a different statement from "no data", and the plots show it as a zero bar. A
    population present in the data but missing from pop_order is not reported at all,
    which is how the ordering has always behaved.
    """
    return [
        {
            'population': pop,
            'n': cases_with_allele[pop],
            'frequency': round(cases_with_allele[pop]/pop_totals[pop], 5)
        }
        for pop in pop_order if pop in pop_totals
    ]

# calculate the frequency that an allele or aminoacid appears in a population, alt a superpopulation
def calculate_frequencies(allele_name, population_type, plot_type):
    db_name_column = allele_column(plot_type)
    pop_order = population_display_order(population_type)

    pop_count = population_totals(population_type)

    cases_with_allele = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.case,
        getattr(ImmuneDiscoverDataModel, population_type),
        db_name_column
        ).where(db_name_column == allele_name).distinct().all()

    pop_with_allele_count = Counter(col[1] for col in cases_with_allele)
    pop_with_allele_count["ALL"] = sum(pop_with_allele_count.values())

    return frequency_entries(pop_count, pop_with_allele_count, pop_order)

def calculate_all_frequencies(population_type, plot_type):
    """Frequencies for every plottable allele at once, as {allele_name: [entries]}.

    Produces exactly what calling calculate_frequencies() once per allele produces, but
    in two queries rather than two per allele. That per-allele loop is what made the
    startup pre-load take tens of minutes: it re-read the whole table for each allele,
    and re-derived the same population totals every time.
    """
    db_name_column = allele_column(plot_type)
    pop_order = population_display_order(population_type)

    # The denominator is deliberately not restricted to the rows selected below. It means
    # "how many individuals are in this population", so it covers the whole cohort - not
    # only the individuals that happen to carry an amino acid allele or a plottable one.
    pop_count = population_totals(population_type)

    # One row per (allele, case, population), so that counting rows per
    # (allele, population) counts distinct cases - what the DISTINCT in
    # calculate_frequencies() does for a single allele.
    distinct_cases = ImmuneDiscoverDataModel.query.with_entities(
        db_name_column.label("allele"),
        ImmuneDiscoverDataModel.case,
        getattr(ImmuneDiscoverDataModel, population_type).label("population"),
        ).filter(db_name_column != None).filter(*plot_selection_criteria()).distinct().subquery()

    cases_per_allele = db.session.query(
        distinct_cases.c.allele,
        distinct_cases.c.population,
        func.count(),
        ).group_by(distinct_cases.c.allele, distinct_cases.c.population).all()

    counts_by_allele = {}
    for allele_name, population, n in cases_per_allele:
        counts_by_allele.setdefault(allele_name, Counter())[population] = n

    data_out = {}
    for allele_name, cases_with_allele in counts_by_allele.items():
        cases_with_allele["ALL"] = sum(cases_with_allele.values())
        data_out[allele_name] = frequency_entries(pop_count, cases_with_allele, pop_order)

    return data_out

# create a .tsv formated table with frequency data for the requested allele/gene and type (genomic or amino acid),
# which can then be downloaded by a user.
def create_frequencies_table(allele_or_gene, plot_type, full_gene = False):
    """Build the downloadable tsv, or return None when there is nothing to report.

    None means the caller asked for an amino acid table by a name that resolves to no
    master amino acid allele, which the resource turns into a 404.
    """
    # Validates plot_type up front. The branches below choose on it with if/elif and if/else,
    # so an unrecognised value silently produced a header-only tsv for a full gene and was
    # treated as amino acid everywhere else.
    allele_column(plot_type)

    alleles = []

    # if full gene is requested, create a query to fetch
    # all alleles of that gene
    # The gene column is filtered on but deliberately not selected: one allele can be
    # stored under more than one gene value (a name like "IGHV3-30,IGHV3-30-5" and a
    # plain "IGHV3-30"), and selecting it would make DISTINCT treat those as different
    # alleles and repeat the allele in the table.
    if full_gene:
        if plot_type == "genomic":
            allele_data = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.db_name,
            ).filter(ImmuneDiscoverDataModel.gene.regexp_match(plot_options_regex(allele_or_gene))).filter(~ImmuneDiscoverDataModel.db_name.contains('_F', autoescape=True)).distinct().all()

            for item in allele_data:
                alleles.append({'allele': item[0]})

        elif plot_type == "aminoacid":
            # db_name_AA_list is a property of db_name_AA, so a plain DISTINCT over the
            # two of them returns one row per amino acid allele. DISTINCT ON is not an
            # option here: it is PostgreSQL only, and is silently dropped on SQLite.
            allele_data = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.db_name_AA,
            ImmuneDiscoverDataModel.db_name_AA_list,
            ).filter(ImmuneDiscoverDataModel.gene.regexp_match(plot_options_regex(allele_or_gene))).filter(ImmuneDiscoverDataModel.db_name_AA != None).distinct().all()

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
                # get_aminoacid_top_allele returns {} when the name is in no row's
                # db_name_AA_list, so subscripting it was a KeyError and a 500 for any
                # schema-valid name that is not a genomic allele with a translation - a
                # flanking variant, an IGHD allele, a typo. The sibling of the {} the
                # IgSNPer endpoint now turns into a 404, in the one path this function
                # still had left.
                top_allele = get_aminoacid_top_allele(allele_name)
                if not top_allele:
                    return None
                allele_name = top_allele['allele_aa']
                aa_list = get_aminoacid_allele_list(allele_name)['aa_allele_list']

        plot_data_subpops = {}
        plot_data_superpops = {}
        if plot_type == "genomic":
            superpop_cache = allele_superpopulation_frequencies
            subpop_cache = allele_population_frequencies
        else:
            superpop_cache = aminoacid_allele_superpopulation_frequencies
            subpop_cache = aminoacid_allele_population_frequencies

        # If running on prod we have population frequencies pre-calculated in dictionaries,
        # use them directly. Only the plottable alleles are pre-calculated though, and a
        # download can legitimately ask for one that is not - a flanking ('_F') variant, or
        # an IGHD/IGHJ allele - so fall back to calculating those here. Under debug and
        # testing nothing is pre-loaded and everything is calculated on demand.
        if not current_app.debug and not current_app.config.get("TESTING") and allele_name in superpop_cache:
            # Copied, not referenced: the loop below writes 'allele', 'superpopulation' and
            # 'collapsed_translated_sequence' into each entry, which would otherwise mutate
            # the cached dictionaries in place and leak those keys into every later response
            # from the frequency plot endpoints.
            plot_data_superpops = [dict(entry) for entry in superpop_cache[allele_name]]
            plot_data_subpops = [dict(entry) for entry in subpop_cache[allele_name]]
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


    