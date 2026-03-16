# Generates mock data for used in tests.

import random
import csv
from constants import ROOT_DIR

def generate_mock_data():
    random.seed(42)

    case_values = [
        "case_ACB_AFR","case_ASW_AFR","case_PJL_SAS","case_STU_SAS","case_ITU_SAS",
        "case_BEB_SAS","case_GIH_SAS","case_CHS_EAS","case_CDX_EAS","case_KHV_EAS",
        "case_CHB_EAS","case_JPT_EAS","case_CLM_AMR","case_GWD_AFR","case_ESN_AFR",
        "case_MSL_AFR","case_YRI_AFR","case_FIN_EUR","case_GBR_EUR","case_IBS_EUR",
        "case_LWK_AFR","case_TSI_EUR","case_MXL_AMR","case_PEL_AMR","case_PUR_AMR"
    ]

    aa_groups = {
        "TEST1-8*01": [
            ("TEST1-8*01", "TEST1-8", "01"),
            ("TEST1-8*02", "TEST1-8", "02"),
            ("TEST1-8*04", "TEST1-8", "04"),
        ],
        "TEST1-69*01_S01": [
            ("TEST1-69*01_S01", "TEST1-69/1-69D", "01_S01")
        ],
        "TEST1-30*01_S01": [
            ("TEST1-30*01_S01", "TEST1-30+", "01_S01")
        ],
        "TEST2-8*01": [
            ("TEST2-8*01", "TEST2-8", "01")
        ],
        "TEST3-4*01": [
            ("TEST3-4*01", "TEST3-4", "01")
        ],
        "SEQTEST1-2*01": [
            ("SEQTEST1-2*01", "SEQTEST1-2", "01")
        ],
        "ALIGNMENTTEST1-2*01": [
            ("ALIGNMENTTEST1-2*01", "ALIGNMENTTEST1-2", "01"),
            ("ALIGNMENTTEST1-2*02", "ALIGNMENTTEST1-2", "02"),
            ("ALIGNMENTTEST1-2*03", "ALIGNMENTTEST1-2", "03"),
        ]
    }

    def aa_list(master):
        return ",".join(v[0] for v in aa_groups[master])

    headers = [
        "cohort","case","db_name","gene","allele","sequence","prefix","suffix","flank_index",
        "count","full_count","IgSNPer_uncommon","IgSNPer_common","IgSNPer_uncommon_str",
        "IgSNPer_common_str","IgSNPer_SNPs","db_name_AA","db_name_AA_list","sequence_AA","file"
    ]

    with open(ROOT_DIR + "/tests/mock_data/in/mock_allele_data.tsv", "w", newline="") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(headers)

        unique = set()

        # 1. Ensure each case value has one row with db_name = TEST1-8*01
        for case in case_values:
            master = "TEST1-8*01"
            db_name, gene, allele = aa_groups[master][0]
            flank = random.choice([0.0, 1.0, 2.0, 3.0, 4.0])
            key = (case, db_name, flank)
            unique.add(key)

            writer.writerow([
                "1KGP", case, db_name, gene, allele, "TESTSEQUENCE01", "PREFIX", "SUFFIX", flank,
                200.0, 200.0, round(random.uniform(0, 1), 2), round(random.uniform(0, 10), 2), "",
                f"(:0,{random.randint(10000000,99999999)});",
                f"rs{random.randint(10000000,99999999)};rs{random.randint(10000000,99999999)};",
                master, aa_list(master), "TESTAMINOACIDSEQUENCE", "filename.tsv.gz"
            ])

        # 2. Add more unique rows to reach at least 60
        while len(unique) < 60:
            case = random.choice(case_values)
            master = random.choice(list(aa_groups.keys()))
            db_name, gene, allele = random.choice(aa_groups[master])
            flank = random.choice([0.0, 1.0, 2.0, 3.0, 4.0])
            key = (case, db_name, flank)
            if key in unique:
                continue
            unique.add(key)

            writer.writerow([
                "1KGP", case, db_name, gene, allele, "TESTSEQUENCE01", "PREFIX", "SUFFIX", flank,
                200.0, 200.0, round(random.uniform(0, 1), 2), round(random.uniform(0, 10), 2), "",
                f"(:0,{random.randint(10000000,99999999)});",
                f"rs{random.randint(10000000,99999999)};rs{random.randint(10000000,99999999)};",
                master, aa_list(master), "TESTAMINOACIDSEQUENCE", "filename.tsv.gz"
            ])

        # 3. add row to test sequence search
        master = "SEQTEST1-2*01"
        db_name, gene, allele = aa_groups[master][0]
        flank = random.choice([0.0, 1.0, 2.0, 3.0, 4.0])
        key = (case, db_name, flank)
        unique.add(key)

        writer.writerow([
            "1KGP", case, db_name, gene, allele, "THISISASEQUENCESEARCHTEST123", "PREFIX", "SUFFIX", flank,
            200.0, 200.0, round(random.uniform(0, 1), 2), round(random.uniform(0, 10), 2), "",
            f"(:0,{random.randint(10000000,99999999)});",
            f"rs{random.randint(10000000,99999999)};rs{random.randint(10000000,99999999)};",
            master, aa_list(master), "TESTAMINOACIDSEQUENCE", "filename.tsv.gz"
        ])

        # 4. add test sequences for alignments
        # one "base" sequence, one sequence with a symmetric number of deletions (3)
        # one sequence with a non-symmetric number of deletions (not divisble by 3)
        alignment_sequences = [
            "CTGGATTCACCTTTACTAGCTCTGCTATGCAGTGGGTGCGACAGGCTCGTGGACAACGCC",
            "CTGGATTCACCTTTACGTATGCTATGCAGTGGGTGCGACAGGCTCGTGGACAACGCC",
            "CTGGATTCACCTTTACTAGCTCTGCTATGCAGTGGGTGCGAAGGCTCGTGGACAACGCC",
        ]
        for i in range(3):
            master = "ALIGNMENTTEST1-2*0" + str(i+1)
            db_name, gene, allele = master, "ALIGNMENTTEST1-2", "0" + str(i+1)
            flank = random.choice([0.0, 1.0, 2.0, 3.0, 4.0])
            key = (case, db_name, flank)
            unique.add(key)

            writer.writerow([
                "1KGP", case, db_name, gene, allele, alignment_sequences[i], "PREFIX", "SUFFIX", flank,
                200.0, 200.0, round(random.uniform(0, 1), 2), round(random.uniform(0, 10), 2), "",
                f"(:0,{random.randint(10000000,99999999)});",
                f"rs{random.randint(10000000,99999999)};rs{random.randint(10000000,99999999)};",
                "ALIGNMENTTEST1-2*01", aa_list("ALIGNMENTTEST1-2*01"), "TESTAMINOACIDSEQUENCE", "filename.tsv.gz"
            ])

    print(f"Generated mock_allele_data.tsv with {len(unique)} unique rows, including one for every case with TEST1-8*01.")
