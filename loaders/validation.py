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
from repositories.filters import is_deletion, plot_selection_criteria
from services.frequencies import subpopulation_order, superpopulation_order

class SourceDataError(Exception):
    """Raised when loaded data breaks an assumption the rest of the app relies on."""

# A row carries no amino acid data exactly when it is one of the kinds that is never
# translated: a flanking region variant, a homozygous deletion, or a D or J gene.
def never_translated():
    return or_(
        ImmuneDiscoverDataModel.db_name.contains('_F', autoescape=True),
        is_deletion(),
        ImmuneDiscoverDataModel.gene.like('IGHD%'),
        ImmuneDiscoverDataModel.gene.like('IGHJ%'),
    )

def sample_names(rows, limit = 10):
    names = sorted(row[0] for row in rows)
    if len(names) <= limit:
        return ", ".join(names)
    return ", ".join(names[:limit]) + ", ... (" + str(len(names) - limit) + " more)"

def amino_acid_coverage_problems():
    """db_name_AA must be null on exactly the rows that are never translated.

    Both directions matter. A translated allele missing db_name_AA drops out of the amino
    acid plots and the translated FASTA download without any error; a flanking or *DEL row
    that has db_name_AA set would appear in them, and would also break the translated
    FASTA branch in services/fasta_generation.py, which has no exclusion of its own and
    relies on db_name_AA being null to leave those rows out.
    """
    problems = []

    missing = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.db_name
        ).filter(ImmuneDiscoverDataModel.db_name_AA == None
        ).filter(~never_translated()).distinct().all()
    if missing:
        problems.append(
            str(len(missing)) + " allele(s) have no db_name_AA but should be translated: "
            + sample_names(missing))

    unexpected = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.db_name
        ).filter(ImmuneDiscoverDataModel.db_name_AA != None
        ).filter(never_translated()).distinct().all()
    if unexpected:
        problems.append(
            str(len(unexpected)) + " allele(s) have db_name_AA set but are never translated"
            " (flanking, *DEL, IGHD or IGHJ): " + sample_names(unexpected))

    return problems

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

    Grouped over the comma-separated components of the gene column rather than its literal
    value, because that is how the resolver matches: plot_options_regex treats a composite
    value like "IGHV3-30,IGHV3-30-5" as naming both genes, and /data/plotoptions offers them
    as two separate options. So the rows one selection can reach span every literal value
    that names that gene, and grouping literally split them into separate groups that each
    looked unambiguous - a new composite row colliding with an existing plain one would have
    passed while making the resolver's answer arbitrary.

    Splitting on commas covers what the resolver matches for whole components, but not quite
    everything it matches: plot_options_regex makes the leading comma optional, so a
    selection that is only the tail of a component - "30-5" against "IGHV3-30-5" - resolves
    too, and no group here holds it. That is left alone deliberately. Every gene name in the
    data is prefixed with its locus, so none is a suffix of another (93 components, 0 such
    pairs), and /data/plotoptions only ever offers whole components - so a selection of that
    shape cannot come from the frontend and cannot collide with a real gene.
    """
    rows = ImmuneDiscoverDataModel.query.with_entities(
        ImmuneDiscoverDataModel.gene,
        ImmuneDiscoverDataModel.allele,
        ImmuneDiscoverDataModel.db_name,
        # The same criteria get_db_name_from_options resolves through, locus restriction
        # included. Excluding only the flanking rows made this stricter than the function it
        # protects: an ambiguous pair in IGHD or IGHJ, which the resolver can never reach and
        # no plot ever shows, would have stopped the service booting.
        ).filter(*plot_selection_criteria()).distinct().all()

    by_selection = {}
    for gene, allele, db_name in rows:
        for component in gene.split(","):
            by_selection.setdefault((component, allele), set()).add(db_name)

    ambiguous = sorted(pair for pair, names in by_selection.items() if len(names) > 1)

    if not ambiguous:
        return []

    return [str(len(ambiguous)) + " gene/allele selection(s) resolve to more than one allele "
            "name, so which one a plot shows is arbitrary: "
            + sample_names([(g + "," + a,) for g, a in ambiguous])]

def validate_loaded_data():
    """Check the loaded data, raising SourceDataError if anything is wrong.

    Every check runs before anything is raised, so one crash reports everything that needs
    fixing rather than one problem per restart.
    """
    problems = (amino_acid_coverage_problems()
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
