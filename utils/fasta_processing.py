# Utility scripts for handling fasta conversions.

import os

# Converts a dict with {allele, sequence} into a string in fasta format
def dict_to_fasta_str(sequence_data):
    fasta_out = ''
    for gene in sequence_data:
        fasta_out += '>' + gene + '\n' + sequence_data[gene] + '\n'
    return fasta_out

# Reads a fasta file and converts the content to a dict with {allele, sequence} format
def fasta_to_dict(fasta_filename):
    fasta_path = os.getcwd() + "/tmp/" + fasta_filename
    sequence_data = {}
    fasta_strings = []
    try:
        with open(fasta_path, 'r') as file:
            for line in file:
                fasta_strings.append(line.strip())
    except FileNotFoundError:
        print("File not found at " + fasta_path)

    current_gene = ""
    current_sequence = ""

    for line in fasta_strings:
        if line.startswith(">"):
            if current_gene:
                sequence_data[current_gene] = current_sequence
            current_gene = line.replace(">", "")
            current_sequence = ""
        else:
            current_sequence += line

    # add last sequence to data
    sequence_data[current_gene] = current_sequence

    # print(sequence_data)
    return sequence_data