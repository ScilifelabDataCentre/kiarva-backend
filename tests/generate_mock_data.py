# Generates the mock dataset used by the tests.
#
# Rules for this dataset:
#
#   * It uses the real naming format (db_name "IGHV1-8*01", gene "IGHV1-8", cases
#     ending in population_superpopulation, ...) and mirrors the column layout and
#     value formats of the real 1KGP tsv in data/compressed/. Several loader quirks
#     are keyed on real gene names and can only be reached with real names, and the
#     frontend E2E fixtures (kiarva-frontend/next-app/tests/fixtures, recorded with
#     next-app/scripts/record-fixtures.sh) query real names too.
#   * It is deliberately minimal: every row exists for a stated reason, written in
#     the comment above it. Do not add rows "for coverage" - if a row is not needed
#     by a test in tests/test_routes.py or by an edge case named below, it does not
#     belong here. Sequences are short fragments rather than full alleles so that
#     alignment and fasta output can be asserted in full.
#
# Both the .tsv and the .zip in tests/mock_data/ are committed, and the .zip is the
# one that counts: load_tsv_to_db() unpacks it over tests/mock_data/in/ on every
# test run, the same way the real data is unpacked in production. Run this module
# (python -m tests.generate_mock_data) after editing it to rewrite both.

import csv
import os
import zipfile

from constants import ROOT_DIR

MOCK_DATA_DIR = ROOT_DIR + "/tests/mock_data/"
TSV_NAME = "mock_allele_data.tsv"

HEADERS = [
    "cohort", "case", "db_name", "gene", "allele", "sequence", "prefix", "suffix",
    "flank_index", "count", "full_count", "IgSNPer_uncommon", "IgSNPer_common",
    "IgSNPer_uncommon_str", "IgSNPer_common_str", "IgSNPer_SNPs", "db_name_AA",
    "db_name_AA_list", "sequence_AA", "file",
]

# The 25 1KGP populations and the superpopulation each belongs to. /data/populationregions
# returns this mapping, and the frequency endpoints need one case per population to
# report a frequency for all of them, so this list sets the size of the dataset.
POPULATIONS = [
    ("ACB", "AFR"), ("ASW", "AFR"), ("ESN", "AFR"), ("GWD", "AFR"), ("LWK", "AFR"),
    ("MSL", "AFR"), ("YRI", "AFR"),
    ("FIN", "EUR"), ("GBR", "EUR"), ("IBS", "EUR"), ("TSI", "EUR"),
    ("CDX", "EAS"), ("CHB", "EAS"), ("CHS", "EAS"), ("JPT", "EAS"), ("KHV", "EAS"),
    ("BEB", "SAS"), ("GIH", "SAS"), ("ITU", "SAS"), ("PJL", "SAS"), ("STU", "SAS"),
    ("CLM", "AMR"), ("MXL", "AMR"), ("PEL", "AMR"), ("PUR", "AMR"),
]

SOURCE_FILE = "1KGP_mock_genotypes.tsv.gz"

# Flanking regions of IGHV1-8, kept in the prefix/suffix columns of every IGHV1-8 row
# and spliced into the sequence of the flanking (_F) row further down.
IGHV1_8_PREFIX = "AGGTGCCCACTCC"
IGHV1_8_SUFFIX = "CACAGTGTGAAA"

# The three IGHV1-8 alleles, chosen to cover the three translation cases of
# services/alignment.py once MAFFT has aligned them against each other.
IGHV1_8_SEQUENCES = {
    # Reference: no gaps, translates straight through.
    "01": "CTGGATTCACCTTTACTAGCTCTGCTATGCAGTGGGTGCGACAGGCTCGTGGACAACGCC",
    # 3 nt shorter, so the alignment gap is a whole codon: the reading frame survives
    # and the translation keeps going past the gap.
    "02": "CTGGATTCACCTTTACGTATGCTATGCAGTGGGTGCGACAGGCTCGTGGACAACGCC",
    # 1 nt shorter, so the alignment gap is a frameshift: translate_nt_to_aa() has to
    # stop at the gap and pad the rest of the amino acid sequence with X.
    "04": "CTGGATTCACCTTTACTAGCTCTGCTATGCAGTGGGTGCGAAGGCTCGTGGACAACGCC",
}

# All three IGHV1-8 alleles collapse into one amino acid allele. That is what makes
# the translated endpoints differ from the genomic ones: an allele missing from a
# population genomically can still be present in translation. (Nucleotide deletions
# would not collapse like this in real data - the deletions above are there for the
# alignment cases, and the collapsing for the amino acid ones.)
IGHV1_8_AA_ALLELE = "IGHV1-8*01"
IGHV1_8_AA_LIST = "IGHV1-8*01,IGHV1-8*02,IGHV1-8*04"
IGHV1_8_AA_SEQUENCE = "LDSPLLALLCSGCDRLVDNA"

# Real IgSNPer values are a float score plus a ';'-separated list of SNPs, where each
# SNP is an rs id followed by its base, count and position.
IGHV1_8_SCORE = "2.0"
IGHV1_8_SNPS = "rs12345678(A:87,106539163);rs87654321(G:204,106539280);"
# An allele with no rs ids of its own still gets one "(:0,<position>)" entry rather
# than an empty list.
NO_RS_SNPS = "(:0,106791243);"

# Not a nucleotide sequence: it is here to be searched for. /data/sequences is the one
# endpoint a user types input into, so the mock hit is recognisable in test output and
# in the frontend fixture. Contains "ESEARCHTES" once (at 14) and "SE" twice (7, 15).
SEARCH_SEQUENCE = "THISISASEQUENCESEARCHTEST123"

# Every population carries IGHV1-8*01 except these three. Without the deviations every
# genomic frequency would come out 1.0 and every superpopulation would look the same,
# which hides mistakes in the frequency maths.
POPULATION_DEVIATIONS = {
    # A different allele of the same gene: IGHV1-8*01 is absent here genomically but
    # present in translation, since both alleles collapse to IGHV1-8*01.
    "TSI": "02",
    "CDX": "04",
    # A homozygous deletion. Only db_name and the columns before it are filled in real
    # *DEL rows, which is where the loader's empty-value handling gets exercised
    # (sequence -> None, flank_index -> -1, IgSNPer/amino acid columns -> None), and
    # where /data/igsnperdata has to answer with an empty score and SNP list. Also the
    # only population without any amino acid data.
    "MSL": "DEL",
}


def row(
    case,
    db_name,
    gene,
    allele,
    sequence="",
    prefix="",
    suffix="",
    flank_index="1.0",
    score="",
    snps="",
    aa_allele="",
    aa_list="",
    aa_sequence="",
    source_file=SOURCE_FILE,
):
    # count/full_count and the IgSNPer_*_str columns are part of the real file but are
    # not read by the loader, so they are kept constant here.
    return {
        "cohort": "1KGP",
        "case": case,
        "db_name": db_name,
        "gene": gene,
        "allele": allele,
        "sequence": sequence,
        "prefix": prefix,
        "suffix": suffix,
        "flank_index": flank_index,
        "count": "200.0" if sequence else "",
        "full_count": "190.0" if sequence else "",
        "IgSNPer_uncommon": score,
        "IgSNPer_common": "1.0" if score else "",
        "IgSNPer_uncommon_str": "",
        "IgSNPer_common_str": snps,
        "IgSNPer_SNPs": snps,
        "db_name_AA": aa_allele,
        "db_name_AA_list": aa_list,
        "sequence_AA": aa_sequence,
        "file": source_file,
    }


def ighv1_8_row(case, allele):
    """One IGHV1-8 row, either a real allele or the *DEL placeholder."""
    if allele == "DEL":
        return row(case, "IGHV1-8*DEL", "IGHV1-8", "DEL", flank_index="", source_file="")
    return row(
        case,
        "IGHV1-8*" + allele,
        "IGHV1-8",
        allele,
        sequence=IGHV1_8_SEQUENCES[allele],
        prefix=IGHV1_8_PREFIX,
        suffix=IGHV1_8_SUFFIX,
        score=IGHV1_8_SCORE,
        snps=IGHV1_8_SNPS,
        aa_allele=IGHV1_8_AA_ALLELE,
        aa_list=IGHV1_8_AA_LIST,
        aa_sequence=IGHV1_8_AA_SEQUENCE,
    )


def mock_rows():
    rows = []

    # One case per population, so that the frequency endpoints cover all 25
    # populations and all 5 superpopulations.
    for population, superpopulation in POPULATIONS:
        case = f"case_{population}_{superpopulation}"
        rows.append(ighv1_8_row(case, POPULATION_DEVIATIONS.get(population, "01")))

    # A second case in TSI. Frequencies count distinct cases, not rows, so this is the
    # only population where the denominator is not 1 and a frequency comes out
    # fractional (IGHV1-8*01 in TSI: 1 of 2 cases).
    rows.append(ighv1_8_row("case2_TSI_EUR", "01"))

    # The flanking sequence of IGHV1-8*01, as its own row. In real IGHV rows only
    # db_name carries the _F suffix - the allele column repeats "01" - and no amino
    # acid data is attached. Every query except /fasta/genomic_fl filters these rows
    # out (db_name NOT LIKE '%_F%'), including /data/db_name, which has to keep
    # answering IGHV1-8*01 for gene+allele "IGHV1-8,01".
    rows.append(
        row(
            "case_ACB_AFR",
            "IGHV1-8*01_F1",
            "IGHV1-8",
            "01",
            sequence=IGHV1_8_PREFIX + IGHV1_8_SEQUENCES["01"] + IGHV1_8_SUFFIX,
            prefix=IGHV1_8_PREFIX,
            suffix=IGHV1_8_SUFFIX,
        )
    )

    # The loader rewrites a handful of real gene names, so each rewrite needs a row.
    # gene "IGHV1-69/1-69D" is stored as "IGHV1-69". The allele also carries a strain
    # suffix, which has to survive into the /data/plotoptions allele list.
    rows.append(
        row(
            "case_CHB_EAS",
            "IGHV1-69*04_S7754",
            "IGHV1-69/1-69D",
            "04_S7754",
            sequence="CAGGTCCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGTCCTCGGTGAAGGTC",
            prefix="AGGTGTCCAGTCC",
            suffix="CACAGTGTGAAA",
            score="4.0",
            snps="rs11417200(C:77,107170005);",
            aa_allele="IGHV1-69*04_S7754",
            aa_list="IGHV1-69*04_S7754",
            aa_sequence="QVQLVQSGAEVKKPGSSVKV",
        )
    )
    # ... and for that gene alone, *DEL rows are dropped by the loader instead of being
    # stored (genes_without_del), so this row must not reach the database at all.
    rows.append(
        row("case_CHB_EAS", "IGHV1-69*DEL", "IGHV1-69/1-69D", "DEL", flank_index="", source_file="")
    )

    # gene "IGHV3-30+" is stored as "IGHV3-30".
    rows.append(
        row(
            "case_JPT_EAS",
            "IGHV3-30*01",
            "IGHV3-30+",
            "01",
            sequence="CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGGGGTCCCTGAGACTC",
            prefix="AGGTGTCCAGTGT",
            suffix="CACAGTGAGGGG",
            score="0.0",
            snps=NO_RS_SNPS,
            aa_allele="IGHV3-30*01",
            aa_list="IGHV3-30*01",
            aa_sequence="QVQLVESGGGVVQPGGSLRL",
        )
    )
    # One db_name can belong to two genes. The loader turns this one into gene
    # "IGHV3-30,IGHV3-30-5" and copies db_name into the allele column, which
    # /data/plotoptions has to split back into two gene options, and which /data/db_name
    # is queried with as an allele.
    rows.append(
        row(
            "case_JPT_EAS",
            "IGHV3-30*02/IGHV3-30-5*02",
            "IGHV3-30+",
            "02",
            sequence="CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGGGGTCCCTGAGACTG",
            prefix="AGGTGTCCAGTGT",
            suffix="CACAGTGAGGGG",
            score="7.0",
            snps="rs10137773(C:251,106791253);",
            aa_allele="IGHV3-30*02/IGHV3-30-5*02",
            # A two-name amino acid list on a db_name that itself contains '*' and '/',
            # which the plot option regex has to escape.
            aa_list="IGHV3-30*02/IGHV3-30-5*02,IGHV3-30*02_S9143",
            aa_sequence="QVQLVESGGGVVQPGGSLRL",
        )
    )
    # A db_name that is rewritten to a gene the gene column never mentions
    # ("IGHV3-30+" -> "IGHV3-30-5"), which is what makes IGHV3-30-5 selectable.
    rows.append(
        row(
            "case_KHV_EAS",
            "IGHV3-30-5*03_S1123",
            "IGHV3-30+",
            "03_S1123",
            sequence="CAGGTGCAGCTGGTGGAGTCTGGGGGAGGCGTGGTCCAGCCTGGGGGGTCCCTGAGACTA",
            prefix="AGGTGTCCAGTGT",
            suffix="CACAGTGAGGGG",
            score="0.0",
            snps=NO_RS_SNPS,
            aa_allele="IGHV3-30-5*03_S1123",
            aa_list="IGHV3-30-5*03_S1123",
            aa_sequence="QVQLVESGGGVVQPGGSLRL",
        )
    )
    # gene "IGHV3-23/3-23D" is stored as "IGHV3-23". This row also carries the sequence
    # searched for by /data/sequences.
    rows.append(
        row(
            "case_GIH_SAS",
            "IGHV3-23*01",
            "IGHV3-23/3-23D",
            "01",
            sequence=SEARCH_SEQUENCE,
            prefix="AGGTGTCCAGTGT",
            suffix="CACAGTGAGGGG",
            score="1.0",
            snps="rs56069819(G:31,106848172);",
            aa_allele="IGHV3-23*01",
            aa_list="IGHV3-23*01",
            aa_sequence="EVQLLESGGGLVQPGGSLRL",
        )
    )

    # D and J genes: a third and fourth gene type for /data/plotoptions and the download
    # page, and rows with no amino acid data at all - only V genes are translated.
    rows.append(
        row(
            "case_PUR_AMR",
            "IGHD1-1*01",
            "IGHD1-1",
            "01",
            sequence="GGTACAACTGGAACGAC",
            score="0.0",
            snps="(:0,106385361);",
        )
    )
    # D and J flanking rows suffix the allele column as well as db_name, unlike the V
    # gene flanking row above, so "01_F1" shows up as an allele of IGHD1-1.
    rows.append(
        row(
            "case_PUR_AMR",
            "IGHD1-1*01_F1",
            "IGHD1-1",
            "01_F1",
            sequence="AACAGCCCCGAGTCACGGTGGGTACAACTGGAACGACCACCGTGAGAAAAACTGTGT",
        )
    )
    rows.append(
        row(
            "case_PEL_AMR",
            "IGHJ6*02",
            "IGHJ6",
            "02",
            sequence="ATTACTACTACTACTACGGTATGGACGTCTGGGGCCAAGGGACCACGGTCACCGTCTCCTCAG",
            score="0.0",
            snps="rs1950943(G:28,106329435);rs2338628(A:44,106329451);",
        )
    )

    return rows


def generate_mock_data():
    rows = mock_rows()

    tsv_path = MOCK_DATA_DIR + "in/" + TSV_NAME
    os.makedirs(os.path.dirname(tsv_path), exist_ok=True)
    with open(tsv_path, "w", newline="") as tsv_file:
        writer = csv.DictWriter(tsv_file, fieldnames=HEADERS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    # The zip is what the test setup actually loads, so keep it in step with the .tsv.
    # Written with a fixed timestamp so that regenerating unchanged data does not
    # produce a different file.
    zip_path = MOCK_DATA_DIR + "compressed/" + TSV_NAME + ".zip"
    os.makedirs(os.path.dirname(zip_path), exist_ok=True)
    with open(tsv_path, "rb") as tsv_file:
        tsv_bytes = tsv_file.read()
    info = zipfile.ZipInfo(TSV_NAME, date_time=(1980, 1, 1, 0, 0, 0))
    info.compress_type = zipfile.ZIP_DEFLATED
    info.external_attr = 0o644 << 16
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr(info, tsv_bytes)

    print(f"Wrote {len(rows)} rows to {tsv_path} and {zip_path}")


if __name__ == "__main__":
    generate_mock_data()
