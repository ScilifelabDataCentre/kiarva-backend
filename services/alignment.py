# Scripts related to sequence alignment using MAFFT (Multiple Alignment Fast
# Fourier Transforms). Hedestam group uses MAFFT to do sequence alignment, so 
# we needed to use MAFFT to reflect the same results. Due to problems with how
# MAFFT handles amino acid sequences however, we first use MAFFT to align
# genetic sequences, then a custom script to transcribe to amino acid sequences
# with gaps.

import os
import shutil
import subprocess
import tempfile

from models.immunediscoverdata import ImmuneDiscoverDataModel
from utils import dict_to_fasta_str, fasta_to_dict
from utils.regex import plot_options_regex
from utils.translation import translate

# Script to apply mafft to a dict of sequences. MAFFT is run from command
# line, expects a .fasta file input and gives .fasta file output, so we employ 
# scripts to convert from dict to fasta, write to file, run MAFFT commands, read file
# and convert back to dict.
def align_with_mafft(seq_dict):
    fasta_str = dict_to_fasta_str(seq_dict)

    # A private directory per call. These were fixed filenames under os.getcwd() + "/tmp/",
    # so two alignments running at once - two gunicorn workers, or two threads in one - wrote
    # the same input file and read the same output file. The loser silently returned a blank
    # record or, worse, the other request's alignment, and both went out as 200.
    tmp_path = tempfile.mkdtemp(prefix="mafft-")
    input_file = os.path.join(tmp_path, "unaligned.fasta")
    output_file = os.path.join(tmp_path, "aligned.fasta")
    try:
        with open(input_file, "w") as f:
            f.write(fasta_str)
        # check=True and no except: a failed alignment used to be printed and then read back
        # out of the output file anyway, which is empty - so MAFFT failing looked exactly
        # like a gene with no sequences and was served as a 200. Letting it raise makes it
        # a 500, which is what a failure to produce the answer is.
        with open(output_file, "w") as out_f:
            subprocess.run(["mafft", "--auto", "--quiet", input_file],
                           stdout=out_f,
                           check=True)
        result = fasta_to_dict(output_file)
    finally:
        shutil.rmtree(tmp_path, ignore_errors=True)

    return dict(sorted(result.items()))

# Function to translate an aligned nucleotide sequence into an amino acid sequence, 
# while looking for frameshifts. 
# Handles frameshifts by putting a cutoff point where they appear and filling rest 
# of string with 'X'.

# Frameshifts are places in a gene sequence where the frame of translation has been 
# shifted enough that it completely changes the rest of the translation. 
def translate_nt_to_aa(sequence):
    aa_full_seq_len = len(sequence) // 3
    seq_to_be_translated = ""
    consecutive_gaps = ""
    tmp_triplet = ""

    # go through the nucleotide sequence as a series of triplets
    for i in range(0, len(sequence), 3):
        if len(sequence)-i >= 3:
            current_triplet = sequence[i:i+3]
            if (i > 2):
                prev_triplet = sequence[i-3:i]
            
            # if the current triplet contains only gaps ('-')
            if (current_triplet.count('-') == 3):
                # if we are still counting gaps from previous triplets, add them to consecutive gaps
                if (len(consecutive_gaps) == 0):
                    seq_to_be_translated += current_triplet
                # else, add them to the string to be translated
                else:
                    consecutive_gaps += current_triplet
            # if the current triplet contains both nucleotides and gaps, separate and save them
            # in tmp_triplet and consecutive_gaps for later
            elif ('-' in current_triplet):
                tmp_triplet += current_triplet.strip('-')
                consecutive_gaps += '-'*current_triplet.count('-')
            # last else is if we find a triplet with no gaps
            else:
                # if consecutive_gaps and/or tmp_triplet are not empty, then we are currently processing
                # previous triplets containing gaps
                if (consecutive_gaps or tmp_triplet):
                    # if the series of gaps we have counted is divisible by 3, then the
                    # alignment is mostly preserved and we can continue
                    if (len(consecutive_gaps) % 3 == 0):
                        # if the previous triplet only had one gap, then we add the gaps first
                        # then the nucleotides, otherwise the other way around
                        if (prev_triplet.count('-') == 1):
                            seq_to_be_translated += consecutive_gaps + tmp_triplet
                        else:
                            seq_to_be_translated += tmp_triplet + consecutive_gaps
                        seq_to_be_translated += current_triplet
                        tmp_triplet = ""
                        consecutive_gaps = ""
                    # if the series of consecutive gaps is not divisible by 3, we've found a frameshift
                    # and break out of the process
                    else:
                        break
                # if the current triplet only contains nucleotides and we are not processing previous
                # triplets, just add them to the string to be translated
                else:
                    seq_to_be_translated += current_triplet

    
    # translate the processed nucleotide sequence using the standard genetic code
    translated_seq = translate(seq_to_be_translated)

    # if we found a frameshift and broke out of the processing above
    # the string will be shortened and the rest of it will be filled with X
    translated_seq += 'X'*(aa_full_seq_len-len(translated_seq))

    return translated_seq
        
def align_sequences(gene):
    sequence_data = ImmuneDiscoverDataModel.query.with_entities(
    ImmuneDiscoverDataModel.db_name,
    ImmuneDiscoverDataModel.sequence,
    ImmuneDiscoverDataModel.gene
    ).filter(ImmuneDiscoverDataModel.gene.regexp_match(plot_options_regex(gene))).filter(~ImmuneDiscoverDataModel.db_name.contains('_F', autoescape=True)).filter(ImmuneDiscoverDataModel.db_name.notlike('%*DEL')).distinct().all()

    nt_seq_dict = {}
    for row in sequence_data:
        nt_seq_dict[row[0]] = row[1]

    aligned_nt = align_with_mafft(nt_seq_dict)

    formated_result = []
    for key in aligned_nt:
        sequence_aa = translate_nt_to_aa(aligned_nt[key].upper())
        formated_result.append({"allele": key, "sequence_nt": aligned_nt[key].upper(), "sequence_aa": sequence_aa})
    return formated_result