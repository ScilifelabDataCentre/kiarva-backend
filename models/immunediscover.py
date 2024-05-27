from db import db

class ImmuneDiscoverModel(db.Model):
    """
    Data output from ImmuneDiscover tool
    """
    __tablename__ = 'immunediscover'

    id = db.Column(db.Integer, primary_key=True)
    well = db.Column(db.Integer, nullable = False)
    case = db.Column(db.String(30), nullable = False)
    db_name = db.Column(db.String(20), nullable = False)
    sequence = db.Column(db.String(500), nullable = False)
    heptamer = db.Column(db.String(10), nullable = False)
    prefix = db.Column(db.String(30), nullable = False)
    full_count = db.Column(db.Integer, nullable = False)
    count = db.Column(db.Integer, nullable = False)
    gene = db.Column(db.String(20), nullable = False)
    full_frequency = db.Column(db.Float, nullable = False)
    frequency = db.Column(db.Float, nullable = False)
    flank_index = db.Column(db.Integer, nullable = False)
    log2_count = db.Column(db.Float, nullable = False)
    file = db.Column(db.String(80), nullable = False)