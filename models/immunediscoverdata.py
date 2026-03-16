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
    allele = db.Column(db.String(20), nullable = False)
    sequence = db.Column(db.Text)
    flank_index = db.Column(db.Integer, nullable = False)
    IgSNPer_uncommon = db.Column(db.Float)
    IgSNPer_SNPs = db.Column(db.String(80))
    db_name_AA = db.Column(db.String(20))
    db_name_AA_list = db.Column(db.String(80))
    sequence_AA = db.Column(db.String(80))
    superpopulation = db.Column(db.String(20), nullable = False)
    population = db.Column(db.String(20), nullable = False)
    loaded_from_tsv = db.Column(db.String(80), nullable = False)
    loaded_at = db.Column(db.String(80), nullable = False)
    __table_args__ = (UniqueConstraint('case', 'db_name', 'gene','flank_index'),)

    def to_dict(self):
        return {column.name: getattr(self, column.name) for column in self.__table__.columns}