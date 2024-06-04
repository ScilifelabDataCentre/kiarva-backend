from marshmallow import Schema, fields

class ImmuneDiscoverDataGetSchema(Schema):
    case = fields.Str(required = True)
    db_name = fields.Str(required = True)
    sequence = fields.Str(required=True)

class ImmuneDiscoverDataUploadSchema(Schema):
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