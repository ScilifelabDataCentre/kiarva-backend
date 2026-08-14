# set up pre-calculated allele population frequencies dict
# for faster loading, only done in prod because it takes a long time
# to load on startup. Running "flask run --debug" allows running without
# pre-loading plots

from constants import allele_superpopulation_frequencies, allele_population_frequencies, aminoacid_allele_superpopulation_frequencies, aminoacid_allele_population_frequencies
from services import calculate_all_frequencies
from datetime import datetime

# The four dictionaries are filled with update() rather than reassigned. constants.py
# exposes them as module-level dictionaries which frequencies.py and the resources
# module import by name, so rebinding the names here would leave those modules holding
# the original empty dictionaries and every plot request would fail.
def load_plot_data_to_dict():
    for description, target, population_type, plot_type in (
        ("genomic alleles by superpopulation", allele_superpopulation_frequencies, "superpopulation", "genomic"),
        ("genomic alleles by population", allele_population_frequencies, "population", "genomic"),
        ("amino acid alleles by superpopulation", aminoacid_allele_superpopulation_frequencies, "superpopulation", "aminoacid"),
        ("amino acid alleles by population", aminoacid_allele_population_frequencies, "population", "aminoacid"),
    ):
        start_time = datetime.now()
        target.update(calculate_all_frequencies(population_type, plot_type))
        elapsed = datetime.now() - start_time
        print("Loaded " + str(len(target)) + " " + description + " in " + str(elapsed))
