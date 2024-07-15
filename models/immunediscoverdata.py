from sqlalchemy import UniqueConstraint
from db import db

class ImmuneDiscoverDataModel(db.Model):
    """
    Data output from ImmuneDiscover tool
    """
    __tablename__ = 'immunediscoverdata'

    id = db.Column(db.Integer, primary_key=True)
    cohort = db.Column(db.String(30), nullable = False)
    case = db.Column(db.String(30), nullable = False)
    db_name = db.Column(db.String(20), nullable = False)
    gene = db.Column(db.String(20), nullable = False)
    sequence = db.Column(db.Text, nullable = False)
    prefix = db.Column(db.String(30), nullable = True)
    heptamer = db.Column(db.String(10), nullable = True)
    pre_heptamer = db.Column(db.String(10), nullable = True)
    post_heptamer = db.Column(db.String(10), nullable = True)
    suffix = db.Column(db.String(30), nullable = True)
    flank_index = db.Column(db.Integer, nullable = True)
    count = db.Column(db.Integer, nullable = False)
    full_count = db.Column(db.Integer, nullable = True)
    superpopulation = db.Column(db.String(20), nullable = False)
    population = db.Column(db.String(20), nullable = False)
    loaded_from_tsv = db.Column(db.String(80), nullable = False)
    loaded_at = db.Column(db.String(80), nullable = False)
    __table_args__ = (UniqueConstraint('case', 'db_name','flank_index'),)