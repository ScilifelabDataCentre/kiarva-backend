# Validation of the source data once it has been loaded into the db.
#
# These checks raise rather than warn. The research group asked to be told about data
# problems immediately: a wrong frequency in a published resource is worse than a service
# that will not start, and a warning in a log nobody reads is not a report. They run in
# every mode - pytest, local dev and production - so a broken TSV is caught while it is
# still being worked on rather than after it ships.
#
# Raising during app creation means the Gunicorn worker never boots and the pod never
# passes its readiness check. With a Deployment and maxUnavailable: 0 the previous pod goes
# on serving, so bad data cannot reach users and CrashLoopBackOff is the alarm. The data is
# baked into the image, so it can only change when a new image ships.

from sqlalchemy import func, or_

from db import db
from models.immunediscoverdata import ImmuneDiscoverDataModel
from repositories.filters import in_plot_loci, is_deletion, locus_of, plot_selection_criteria
from services.frequencies import subpopulation_order, superpopulation_order

class SourceDataError(Exception):
    """Raised when loaded data breaks an assumption the rest of the app relies on."""

# A row carries no amino acid data exactly when it is one of the kinds that is never
# translated: a flanking region variant, a homozygous deletion, or a gene outside the V loci.
#
# The last of those is derived from PLOT_LOCI rather than listing IGHD and IGHJ, which is
# what the D and J genes in the current data happen to be. Spelling them out made the two
# definitions hand-maintained complements of each other: a TRGJ row appearing in future
# source data would have been outside PLOT_LOCI, so never plotted, and simultaneously
# outside the IGHD/IGHJ list, so expected to be translated - and the service would have
# refused to boot reporting "should be translated" about a gene that never is.
def never_translated():
    return or_(
        ImmuneDiscoverDataModel.db_name.contains('_F', autoescape=True),
        is_deletion(),
        ~in_plot_loci(),
    )

def sample_names(rows, limit = 10):
    names = sorted(row[0] for row in rows)
    if len(names) <= limit:
        return ", ".join(names)
    return ", ".join(names[:limit]) + ", ... (" + str(len(names) - limit) + " more)"

def amino_acid_coverage_problems(unknown_loci = ()):
    """db_name_AA must be null on exactly the rows that are never translated.

    Both directions matter. A translated allele missing db_name_AA drops out of the amino
    acid plots and the translated FASTA download without any error; a flanking or *DEL row
    that has db_name_AA set would appear in them, and would also break the translated
    FASTA branch in services/fasta_generation.py, which has no exclusion of its own and
    relies on db_name_AA being null to leave those rows out.

    Rows of an unrecognised locus are excluded, because locus_problems() already reports
    them and this check would say something true but unhelpful about the same rows. Only
    those rows: suppressing the whole check whenever any unknown locus exists would hide a
    genuine missing db_name_AA on a known allele until the unrelated locus was dealt with,
    which is the extra round trip reporting everything at once exists to avoid.
    """
    problems = []
    known_locus_rows = (ImmuneDiscoverDataModel.gene.notin_(unknown_loci)
                        if unknown_loci else True)

    missing = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.db_name
        ).filter(ImmuneDiscoverDataModel.db_name_AA == None
        ).filter(known_locus_rows
        ).filter(~never_translated()).distinct().all()
    if missing:
        problems.append(
            str(len(missing)) + " allele(s) have no db_name_AA but should be translated: "
            + sample_names(missing))

    unexpected = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.db_name
        ).filter(ImmuneDiscoverDataModel.db_name_AA != None
        ).filter(known_locus_rows
        ).filter(never_translated()).distinct().all()
    if unexpected:
        problems.append(
            str(len(unexpected)) + " allele(s) have db_name_AA set but are never translated"
            " (a flanking variant, a *DEL, or a gene outside the plotted V loci): "
            + sample_names(unexpected))

    return problems

def unrecognised_loci():
    """Gene names in the data whose locus is not one the app knows how to present."""
    genes = {row[0] for row in ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.gene).distinct().all()}
    return sorted(gene for gene in genes if locus_of(gene) is None)

def locus_problems(unknown):
    """Every locus in the data must be one the app knows how to present.

    Deriving "never translated" from PLOT_LOCI means an unrecognised locus is treated as
    untranslated and left out of the plots without complaint, which is consistent but
    silent. A locus nobody has decided about is a change the research group should have
    announced, so it is reported here - with the right diagnosis, rather than as a
    translation problem.
    """
    if not unknown:
        return []

    return ["Gene(s) belonging to no known locus: " + sample_names([(g,) for g in unknown])
            + ". Add the locus to KNOWN_LOCI in repositories/filters.py once it is decided "
            "whether it should be plotted, which also decides whether its alleles are "
            "expected to be translated"]

def igsnper_consistency_problems():
    """An allele's IgSNPer values must not differ between the rows carrying it.

    The researchers look these up per allele name, so cohort cannot change them - which is
    why get_igSNPer_data() can answer from the first row it gets back. It queries with
    .distinct() and no ORDER BY, so if one allele ever carried two different values the row
    that won would be whichever the query plan returned first, and a real score could
    disappear silently. Zero alleles diverge today; this is what keeps that true.
    """
    def distinct_values(column):
        return func.count(func.distinct(func.coalesce(column, "~")))

    divergent = db.session.query(ImmuneDiscoverDataModel.db_name).group_by(
        ImmuneDiscoverDataModel.db_name
        ).having(or_(distinct_values(ImmuneDiscoverDataModel.IgSNPer_uncommon) > 1,
                     distinct_values(ImmuneDiscoverDataModel.IgSNPer_SNPs) > 1)).all()

    if not divergent:
        return []

    return [str(len(divergent)) + " allele(s) carry more than one set of IgSNPer values, so "
            "which one is reported depends on the query plan: " + sample_names(divergent)]

def population_problems():
    """Every population in the data must appear in the hardcoded display order.

    The frequency endpoints report populations in the order the research group asked for,
    and a value missing from those lists is silently left out of every plot rather than
    reported - so an unrecognised population is most likely a typo in the source data.
    """
    problems = []

    for column, order, label in (
        (ImmuneDiscoverDataModel.superpopulation, superpopulation_order, "superpopulation"),
        (ImmuneDiscoverDataModel.population, subpopulation_order, "population"),
    ):
        in_data = {row[0] for row in
                   ImmuneDiscoverDataModel.query.with_entities(column).distinct().all()}
        unknown = sorted(in_data - set(order))
        if unknown:
            problems.append(
                "Unknown " + label + "(s) " + sample_names([(u,) for u in unknown])
                + ", which are not in the "
                "display order in services/frequencies.py and would be left out of every plot")

    return problems

def allele_resolution_problems():
    """(gene, allele) must identify exactly one db_name once flanking rows are excluded.

    get_db_name_from_options resolves a plot selection by taking the first row it gets
    back, which is only correct because this holds. It does not hold with flanking variants
    included - they share their parent's allele value - which is why they are excluded from
    the plot and MSA paths.
    """
    distinct_db_names = func.count(func.distinct(ImmuneDiscoverDataModel.db_name))

    # The same criteria get_db_name_from_options resolves through, locus restriction
    # included. Excluding only the flanking rows made this stricter than the function it
    # protects: an ambiguous pair in IGHD or IGHJ, which the resolver can never reach and no
    # plot ever shows, would have stopped the service booting.
    ambiguous = db.session.query(
        ImmuneDiscoverDataModel.gene,
        ImmuneDiscoverDataModel.allele,
        ).filter(*plot_selection_criteria()
        ).group_by(ImmuneDiscoverDataModel.gene, ImmuneDiscoverDataModel.allele
        ).having(distinct_db_names > 1).all()

    if not ambiguous:
        return []

    pairs = [(gene + "," + allele,) for gene, allele in ambiguous]
    return [str(len(pairs)) + " gene/allele selection(s) resolve to more than one allele "
            "name, so which one a plot shows is arbitrary: " + sample_names(pairs)]

def validate_loaded_data():
    """Check the loaded data, raising SourceDataError if anything is wrong.

    Every check runs before anything is raised, so one crash reports everything that needs
    fixing rather than one problem per restart - but each problem is reported once. An
    unrecognised locus makes its rows look untranslated too, so reporting the locus and
    stopping there says one true thing rather than two overlapping ones about the same row.
    """
    unknown_loci = unrecognised_loci()

    problems = (locus_problems(unknown_loci)
                + amino_acid_coverage_problems(unknown_loci)
                + igsnper_consistency_problems()
                + population_problems()
                + allele_resolution_problems())

    if problems:
        raise SourceDataError(
            "Loaded data breaks " + str(len(problems)) + " assumption(s) this app relies "
            "on. These need reporting to the research group so the source data can be "
            "corrected:\n  - " + "\n  - ".join(problems))

    print("Source data validated: " + str(
        ImmuneDiscoverDataModel.query.with_entities(ImmuneDiscoverDataModel.id).count()
    ) + " rows.", flush = True)
