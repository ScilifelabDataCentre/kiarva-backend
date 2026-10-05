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
    full scan of the table, and used to be re-run once per allele during the startup
    pre-load - thousands of times, for one unchanging result.

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

def calculate_all_frequencies(population_type, plot_type):
    """Frequencies for every plottable allele at once, as {allele_name: [entries]}.

    Two queries for the whole dataset. This replaced a loop that called a single-allele
    version once per allele, which re-read the whole table each time and re-derived the
    same population totals with it - tens of minutes of startup for one unchanging result.

    Two things here are easy to "tidy" into being wrong, and both are silent:

    - The denominator must stay unfiltered. Restricting it to the rows selected below -
      to individuals who carry an amino acid allele, say - inflates every amino acid
      frequency. It reads like an oversight because the filter is right there on the
      numerator, and on the real dataset the two happen to agree, because every one of the
      2486 cases carries amino acid data. They do not agree in general: an individual who
      carries only a deletion belongs in the denominator uncounted.
    - A population no case carries the allele in must still be reported, as n=0. GROUP BY
      returns no row for those at all, and 1094 of the 1450 genomic alleles are absent from
      at least one superpopulation. In a plot a missing bar and a zero bar say different
      things.

    test_aminoacid_populationfrequencies and test_populationfrequencies are what hold both
    of these: their exact ALL and MSL figures are 25/26 and 0/1 precisely because one case
    in the mock data carries only a deletion.
    """
    db_name_column = allele_column(plot_type)
    pop_order = population_display_order(population_type)

    # The denominator is deliberately not restricted to the rows selected below. It means
    # "how many individuals are in this population", so it covers the whole cohort - not
    # only the individuals that happen to carry an amino acid allele or a plottable one.
    pop_count = population_totals(population_type)

    # One row per (allele, case, population), so that counting rows per
    # (allele, population) counts distinct cases rather than rows: an individual appears
    # once per gene and flank position, not once overall.
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
    """Build the downloadable tsv, or return None if there is nothing to report.

    Everything reported comes from the dictionaries pre-calculated at startup, so the
    alleles this can describe are exactly the alleles that can be plotted. None means the
    caller asked for something outside that - an allele or gene that is not plottable, or
    not in the data at all - and the resource turns it into a 404.
    """
    # Validated up front. The branches below choose on plot_type with if/elif and if/else,
    # so an unrecognised value silently produced a header-only tsv for a full gene and was
    # treated as amino acid everywhere else.
    allele_column(plot_type)

    if plot_type == "genomic":
        superpop_cache = allele_superpopulation_frequencies
        subpop_cache = allele_population_frequencies
    else:
        superpop_cache = aminoacid_allele_superpopulation_frequencies
        subpop_cache = aminoacid_allele_population_frequencies
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
            ).filter(ImmuneDiscoverDataModel.gene.regexp_match(plot_options_regex(allele_or_gene))).filter(*plot_selection_criteria()).distinct().all()

            for item in allele_data:
                alleles.append({'allele': item[0]})

        elif plot_type == "aminoacid":
            # db_name_AA_list is a property of db_name_AA, so a plain DISTINCT over the
            # two of them returns one row per amino acid allele. DISTINCT ON is not an
            # option here: it is PostgreSQL only, and is silently dropped on SQLite.
            allele_data = ImmuneDiscoverDataModel.query.with_entities(
            ImmuneDiscoverDataModel.db_name_AA,
            ImmuneDiscoverDataModel.db_name_AA_list,
            ).filter(ImmuneDiscoverDataModel.gene.regexp_match(plot_options_regex(allele_or_gene))).filter(ImmuneDiscoverDataModel.db_name_AA != None).filter(*plot_selection_criteria()).distinct().all()

            for item in allele_data:
                # split on second item to get aa_list in form ['aa1','aa2] instead of 'aa1,aa2'
                alleles.append({'allele': item[0], 'aa_list': item[1].split(",")})

    # if single allele, just use requested allele name
    elif plot_type == "genomic":
        alleles = [{'allele': allele_or_gene}]

    else:
        # The amino acid download is asked for by genomic allele name - that is what the
        # frontend has resolved from the plot selection - so it has to be translated to the
        # master amino acid allele it collapses into. 243 of the 732 plottable alleles are
        # not their own master, so resolving this after the guard below rather than before it
        # made every one of them a 404.
        top_allele = get_aminoacid_top_allele(allele_or_gene)
        # .get('allele_aa') rather than the dict itself, which is truthy even when the master
        # is null - the shape of a row carrying db_name_AA_list with no db_name_AA, which
        # validate_loaded_data() rejects at startup. The cache guard below already answers
        # None for that, since None is not a key of the pre-calculated dictionaries; said
        # here as well so this branch states its own precondition instead of resting on a
        # later test of something else.
        if not top_allele.get('allele_aa'):
            return None
        alleles = [{'allele': top_allele['allele_aa'],
                    'aa_list': get_aminoacid_allele_list(top_allele['allele_aa'])['aa_allele_list']}]

    # A gene with no plottable alleles, or an allele that is not one, has no table. Checked
    # against the resolved names, which is what the caches are keyed on.
    if not alleles or any(entry['allele'] not in superpop_cache for entry in alleles):
        return None

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
            aa_list = allele_data['aa_list']

        # Copied, not referenced: the loop below writes 'allele', 'superpopulation' and
        # 'collapsed_translated_sequence' into each entry, which would otherwise mutate the
        # pre-calculated dictionaries in place and leak those keys into every later response
        # from the frequency plot endpoints.
        plot_data_superpops = [dict(entry) for entry in superpop_cache[allele_name]]
        plot_data_subpops = [dict(entry) for entry in subpop_cache[allele_name]]

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


    