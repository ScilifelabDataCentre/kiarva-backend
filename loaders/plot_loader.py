# set up pre-calculated allele population frequencies dict
# for faster loading, only done in prod because it takes a long time
# to load on startup. Running "flask run --debug" allows running without
# pre-loading plots

from models.immunediscoverdata import ImmuneDiscoverDataModel
from constants import allele_superpopulation_frequencies, allele_population_frequencies, aminoacid_allele_superpopulation_frequencies, aminoacid_allele_population_frequencies
from services import calculate_frequencies
from utils import print_progress_bar
from datetime import datetime

def load_plot_data_to_dict():
    alleles_in_db = ImmuneDiscoverDataModel.query.with_entities(ImmuneDiscoverDataModel.db_name).distinct().all()
    alleles = [allele[0] for allele in alleles_in_db]

    aminoacid_alleles_in_db = ImmuneDiscoverDataModel.query.with_entities(ImmuneDiscoverDataModel.db_name_AA).distinct().all()
    aminoacid_alleles = [aminoacid_allele[0] for aminoacid_allele in aminoacid_alleles_in_db]

    nr_of_alleles = len(alleles_in_db)
    nr_of_aa_alleles = len(aminoacid_alleles_in_db)
    # calculate frequencies and add results to dictionaries for instant fetching in prod
    print("Loading genomic alleles...")
    start_time = datetime.now()
    print("Started at " + str(start_time))
    allele_counter_increment = 0
    print_progress_bar(allele_counter_increment, nr_of_alleles, prefix = 'Progress:', suffix = 'Complete', length = 50)
    progress_bar_last_printed = datetime.now()
    for allele in alleles:
        allele_counter_increment += 1

        progress_timer_current = datetime.now()

        if (progress_timer_current-progress_bar_last_printed).total_seconds() > 10:
            print("")
            print_progress_bar(allele_counter_increment, nr_of_alleles, prefix = 'Progress:', suffix = 'Complete', length = 50)
            progress_bar_last_printed = progress_timer_current

        allele_superpopulation_frequencies[allele] = calculate_frequencies(allele, "superpopulation", "genomic")
        allele_population_frequencies[allele] = calculate_frequencies(allele, "population", "genomic")

    finish_time = datetime.now()
    print("Loading genomic alleles complete at " + str(finish_time))
    print("Time elapsed: " + str(finish_time-start_time))
    
    print("Loading amino acid alleles...")
    start_time = datetime.now()
    print("Started at " + str(start_time))
    aa_allele_counter_increment = 0
    print_progress_bar(aa_allele_counter_increment, nr_of_aa_alleles, prefix = 'Progress:', suffix = 'Complete', length = 50)
    progress_bar_last_printed = datetime.now()
    for aminoacid_allele in aminoacid_alleles:
        aa_allele_counter_increment += 1

        progress_timer_current = datetime.now()

        if (progress_timer_current-progress_bar_last_printed).total_seconds() > 10:
            print("")
            print_progress_bar(aa_allele_counter_increment, nr_of_aa_alleles, prefix = 'Progress:', suffix = 'Complete', length = 50)
            progress_bar_last_printed = progress_timer_current
        aminoacid_allele_superpopulation_frequencies[aminoacid_allele] = calculate_frequencies(aminoacid_allele, "superpopulation", "aminoacid")
        aminoacid_allele_population_frequencies[aminoacid_allele] = calculate_frequencies(aminoacid_allele, "population", "aminoacid")

    finish_time = datetime.now()
    print("Loading amino acid alleles complete at " + str(finish_time))
    print("Time elapsed: " + str(finish_time-start_time))