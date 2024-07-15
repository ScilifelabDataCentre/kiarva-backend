from marshmallow import Schema, fields

class ImmuneDiscoverDataGetAllSchema(Schema):
    id = fields.Str(dump_only=True)
    well = fields.Str()
    case = fields.Str(required = True)
    db_name = fields.Str(required = True)
    sequence = fields.Str(required=True)
    heptamer = fields.Str(required = False)
    prefix = fields.Str(required = False)
    full_count = fields.Str()
    count = fields.Str(required=True)
    gene = fields.Str(required = True)
    full_frequency = fields.Str()
    frequency = fields.Str()
    flank_index = fields.Str()
    log2_count = fields.Str()
    file = fields.Str()
    superpopulation = fields.Str(required = True)
    population = fields.Str(required = True)
    loaded_from_tsv = fields.Str(required = True)
    loaded_at = fields.Str(required = True)

class ImmuneDiscoverDataFrequencySchema(Schema):
    population = fields.Str()
    frequency = fields.Float()
    n = fields.Int()

class ImmuneDiscoverPopulationRegionSchema(Schema):
    superpopulation = fields.Str()
    population = fields.Str()