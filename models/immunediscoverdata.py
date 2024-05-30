from db import db

class ImmuneDiscoverDataModel(db.Model):
    """
    Data output from ImmuneDiscover tool
    """
    __tablename__ = 'immunediscoverdata'

    id = db.Column(db.Integer, primary_key=True)
    well = db.Column(db.Integer, nullable = True)
    case = db.Column(db.String(30), nullable = False)
    db_name = db.Column(db.String(20), nullable = False)
    sequence = db.Column(db.Text, nullable = False)
    heptamer = db.Column(db.String(10), nullable = True)
    prefix = db.Column(db.String(30), nullable = True)
    full_count = db.Column(db.Integer, nullable = True)
    count = db.Column(db.Integer, nullable = False)
    gene = db.Column(db.String(20), nullable = False)
    full_frequency = db.Column(db.Float, nullable = True)
    frequency = db.Column(db.Float, nullable = True)
    flank_index = db.Column(db.Integer, nullable = True)
    log2_count = db.Column(db.Float, nullable = True)
    file = db.Column(db.String(80), nullable = True)
    load_tsv_metadata = db.Column(db.String(80), nullable = False)