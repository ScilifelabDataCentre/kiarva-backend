# Shared query criteria for the plot and MSA data paths.
#
# Two categories of row are loaded into the db but must never be offered as a
# plot or alignment selection:
#
# - Flanking-region variants, which carry a '_F' suffix on db_name. The research
#   group has specified these are not to be shown in plots or alignments. They
#   stay in the db because the "Flanking genomic" FASTA download is built from
#   exactly these rows (services/fasta_generation.py, type='genomic_fl').
# - IGHD and IGHJ genes. Out of the IGH loci only IGHV is plotted, per the
#   research group; TRGV is the other plotted locus. IGHD/IGHJ rows are still
#   used by the FASTA downloads.
#
# Excluding '_F' has a second effect worth knowing about: it makes (gene, allele)
# resolve to exactly one db_name. With flanking variants included, 553 of 732
# IGHV/TRGV (gene, allele) pairs are ambiguous, because a flanking variant shares
# the 'allele' value of its parent for these loci - IGHV1-2*02_F1 has allele '02',
# the same as IGHV1-2*02. See get_db_name_from_options in repositories/allele.py.

from sqlalchemy import or_

from models.immunediscoverdata import ImmuneDiscoverDataModel

# Loci offered as plot / MSA selections
PLOT_LOCI = ("IGHV", "TRGV")

def plot_selection_criteria():
    """Criteria restricting a query to the rows eligible for plots and MSA.

    Splat into an existing filter chain: .filter(*plot_selection_criteria())
    """
    return (
        ~ImmuneDiscoverDataModel.db_name.contains('_F', autoescape=True),
        or_(*(ImmuneDiscoverDataModel.gene.like(locus + '%') for locus in PLOT_LOCI)),
    )
