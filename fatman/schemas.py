
from flask_marshmallow import Marshmallow
from webargs import fields

from . import app
from .models import (
    Calculation,
    CalculationCollection,
    Structure,
    StructureSet,
    BasisSet,
    Code,
    Task2,
    Artifact,
    Command,
    TestResult2,
    )


ma = Marshmallow(app)


class ArtifactSchema(ma.Schema):
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('artifactresource', aid='<id>'),
        'collection': ma.AbsoluteURLFor('artifactlistresource'),
        'download': ma.AbsoluteURLFor('artifactdownloadresource', aid='<id>'),
        })

    class Meta:
        model = Artifact
        fields = ('id', 'name', '_links')


class BasisSetSchema(ma.ModelSchema):
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('basissetresource', bid='<id>'),
        'collection': ma.AbsoluteURLFor('basissetlistresource'),
        })

    family = fields.Str(attribute='family.name')

    class Meta:
        model = BasisSet
        exclude = ('calculations', 'calculation_associations', )


class CalculationCollectionSchema(ma.ModelSchema):
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('calculationcollectionresource', ccid='<id>'),
        'collection': ma.AbsoluteURLFor('calculationcollectionlistresource'),
        })

    class Meta:
        model = CalculationCollection
        strict = True
        fields = ('_links', 'id', 'name', 'desc', )


class CalculationBasisSetAssociationSchema(ma.ModelSchema):
    basis_set = fields.Nested(
        BasisSetSchema(exclude=('basis', ))
        )
    type = fields.Str(attribute='btype')


class CalculationListSchema(ma.ModelSchema):
    id = fields.UUID()
    collection = fields.Str(attribute='collection.name')
    code = fields.Str(attribute='code.name')
    structure = fields.Str(attribute='structure.name')
    test = fields.Str(attribute='test.name')
    results_available = fields.Bool()
    current_task = fields.Nested('Task2ListSchema')
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('calculationresource', cid='<id>'),
        'collection': ma.AbsoluteURLFor('calculationlistresource'),
        'tasks': ma.AbsoluteURLFor('calculationtask2listresource', cid='<id>'),
        })


class TestResultSchema(ma.ModelSchema):
    _links = ma.Hyperlinks({
        'collection': ma.AbsoluteURLFor('testresultlistresource'),
        'self': ma.AbsoluteURLFor('testresultresource', trid='<id>'),
        })

    test = fields.Str(attribute='test.name')
    calculations = fields.Nested(CalculationListSchema, many=True)
    collections = fields.Str(attribute='collections.name')

    class Meta:
        model = TestResult2


class CalculationSchema(CalculationListSchema):
    id = fields.UUID()
    collection = fields.Str(attribute='collection.name')
    code = fields.Str(attribute='code.name')
    structure = fields.Nested('StructureSchema', only=('id', 'name', '_links', ))
    test = fields.Str(attribute='test.name')
    tasks = fields.Nested('Task2ListSchema', many=True)
    testresults = fields.Nested('TestResultSchema', many=True, exclude=('calculations',))

    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('calculationresource', cid='<id>'),
        'collection': ma.AbsoluteURLFor('calculationlistresource'),
        'tasks': ma.AbsoluteURLFor('calculationtask2listresource', cid='<id>'),
        })

    basis_sets = fields.Nested(
        CalculationBasisSetAssociationSchema,
        attribute='basis_set_associations',
        many=True)

    class Meta:
        model = Calculation
        exclude = (
            'basis_set_associations',
            'tasks_query',
            'testresults_query',
            )


class CalculationListActionSchema(ma.Schema):
    generateResults = fields.Nested(
        {'update': fields.Boolean(missing=False)},
        strict=True)

    class Meta:
        strict = True


class Task2ListSchema(ma.ModelSchema):
    id = fields.UUID()
    status = fields.Str(attribute='status.name')
    ctime = fields.DateTime()
    mtime = fields.DateTime()
    machine = fields.Str(attribute='machine.name')
    priority = fields.Int()

    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('task2resource', tid='<id>'),
        'collection': ma.AbsoluteURLFor('task2listresource'),
        'uploads': ma.AbsoluteURLFor('task2uploadresource', tid='<id>'),
        })


class Task2Schema(Task2ListSchema):
    calculation = fields.Nested(CalculationListSchema)
    infiles = fields.Nested(ArtifactSchema, many=True)
    outfiles = fields.Nested(ArtifactSchema, many=True)
    data = fields.Dict()

    class Meta:
        model = Task2
        exclude = ('artifacts', )


class StructureSchema(ma.ModelSchema):
    calculations = fields.Nested(CalculationListSchema, many=True,
                                 exclude=('structure', ))

    sets = fields.DelimitedList(fields.Str())
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('structureresource_v2', sid='<id>'),
        'collection': ma.AbsoluteURLFor('structurelistresource_v2'),
        })

    replaced_by = fields.Nested(
        'StructureSchema',
        exclude=('ase_structure', 'name', ))

    class Meta:
        model = Structure
        # tests is an old table
        exclude = ('tests', 'replaced_by_id', 'calculations', 'replaced', )


class StructureSetSchema(ma.ModelSchema):
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('structuresetresource', name='<name>'),
        'collection': ma.AbsoluteURLFor('structuresetlistresource'),
        'calculations': ma.AbsoluteURLFor('structuresetcalculationslistresource', name='<name>'),
        })

    class Meta:
        model = StructureSet
        exclude = ('id', 'structures', )


class TestResultListActionSchema(ma.Schema):
    generate = fields.Nested(
        {'update': fields.Boolean(missing=False)},
        strict=True)

    class Meta:
        strict = True


class CodeListSchema(ma.ModelSchema):
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('coderesource', cid='<id>'),
        'collection': ma.AbsoluteURLFor('codelistresource'),
        })

    id = fields.UUID()
    name = fields.Str()


class CodeCommandListSchema(ma.ModelSchema):
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('codecommandresource', cid='<code_id>', mid='<machine_id>'),
        'collection': ma.AbsoluteURLFor('codecommandlistresource', cid='<code_id>'),
        })

    machine_id = fields.UUID()
    machine_name = fields.Str(attribute='machine.name')
    machine_shortname = fields.Str(attribute='machine.shortname')


class CodeCommandSchema(ma.ModelSchema):
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('codecommandresource', cid='<code_id>', mid='<machine_id>'),
        'collection': ma.AbsoluteURLFor('codecommandlistresource', cid='<code_id>'),
        })

    code = fields.Nested(CodeListSchema)

    class Meta:
        model = Command
        strict = True


class CodeSchema(ma.ModelSchema):
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('coderesource', cid='<id>'),
        'collection': ma.AbsoluteURLFor('codelistresource'),
        })

    commands = fields.Nested(CodeCommandListSchema, many=True)

    class Meta:
        model = Code

        exclude = ('calculations', 'default_settings', 'task_runtime_settings', )
