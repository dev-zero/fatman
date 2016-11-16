
import bz2
from urllib.parse import urlsplit
from os.path import basename
from io import TextIOWrapper, BytesIO

from flask import make_response
from flask_marshmallow import Marshmallow
from flask_restful import Api, Resource
from webargs import fields, ValidationError
from webargs.flaskparser import (
    use_kwargs,
    parser,
    abort,
    )
from sqlalchemy.orm import contains_eager

from ase import io as ase_io
import numpy as np

from . import app, db, resultfiles
from .models import (
    Calculation,
    CalculationCollection,
    CalculationDefaultSettings,
    Structure,
    StructureSet,
    BasisSet,
    BasisSetFamily,
    CalculationBasisSet,
    Pseudopotential,
    PseudopotentialFamily,
    Test,
    Code,
    Task2,
    TaskStatus,
    Task2Artifact,
    Artifact,
    Machine,
    Command,
    )
from .tools import Json2Atoms, Atoms2Json
from .tools.cp2k import dict2cp2k, mergedicts

ma = Marshmallow(app)


def must_exist_in_db(model, field='id'):
    def check(model_id):
        if not model.query.filter_by(**{field: model_id}).first():
            raise ValidationError('The specified {} does not exist'.format(
                model.__name__))

    return check


class ArtifactSchema(ma.Schema):
    _links = ma.Hyperlinks({
        'self': ma.URLFor('artifactresource', aid='<id>'),
        'collection': ma.URLFor('artifactlistresource'),
        'download': ma.URLFor('artifactdownloadresource', aid='<id>'),
        })

    class Meta:
        model = Artifact
        fields = ('id', 'name', '_links')


class ArtifactListResource(Resource):
    def get(self):
        schema = ArtifactSchema(many=True)
        return schema.jsonify(Artifact.query.all())


class ArtifactResource(Resource):
    def get(self, aid):
        schema = ArtifactSchema()
        return schema.jsonify(Artifact.query.get_or_404(aid))


class ArtifactDownloadResource(Resource):
    def get(self, aid):
        artifact = Artifact.query.get_or_404(aid)

        scheme, nwloc, path, _, _ = urlsplit(artifact.path)

        if scheme == 'fkup' and nwloc == 'results':
            # Artifact.name contains a full path, including possible subdirs
            filename = basename(artifact.name)

            compressed = artifact.mdata.get('compressed', None)
            if compressed == 'bz2':
                with bz2.BZ2File(resultfiles.path(path[1:])) as infile:
                    response = make_response(infile.read())
                    response.headers["Content-Disposition"] = \
                        "attachment; filename={:}".format(filename)
                    response.headers["Content-Type"] = "text/plain"
                    return response

            elif compressed is None:
                with open(resultfiles.path(path[1:]), 'rb') as infile:
                    response = make_response(infile.read())
                    response.headers["Content-Disposition"] = \
                        "attachment; filename={:}".format(filename)
                    response.headers["Content-Type"] = "text/plain"
                    return response
            else:
                app.logger.error("unknown compression scheme found: %s",
                                 compressed)
                abort(500)

        app.logger.warning("{} contains unknown path '{path}'".format(
            artifact, **artifact.__dict__))
        abort(500)


class CalculationCollectionSchema(ma.ModelSchema):
    _links = ma.Hyperlinks({
        'self': ma.URLFor('calculationcollectionresource', ccid='<id>'),
        'collection': ma.URLFor('calculationcollectionlistresource'),
        })

    class Meta:
        model = CalculationCollection
        strict = True
        fields = ('_links', 'id', 'name', 'desc', )


class CalculationCollectionListResource(Resource):
    def get(self):
        schema = CalculationCollectionSchema(many=True)
        return schema.jsonify(CalculationCollection.query.all())


class CalculationCollectionResource(Resource):
    def get(self, ccid):
        schema = CalculationCollectionSchema()
        return schema.jsonify((CalculationCollection.query.get_or_404(ccid)))


class CalculationListSchema(ma.ModelSchema):
    id = fields.UUID()
    collection = fields.Str(attribute='collection.name')
    code = fields.Str(attribute='code.name')
    structure = fields.Str(attribute='structure.name')
    test = fields.Str(attribute='test.name')
    _links = ma.Hyperlinks({
        'self': ma.URLFor('calculationresource', cid='<id>'),
        'collection': ma.URLFor('calculationlistresource'),
        'tasks': ma.URLFor('calculationtask2listresource', cid='<id>'),
        })


class BasisSetSchema(ma.ModelSchema):
    family = fields.Str(attribute='family.name')

    class Meta:
        model = BasisSet


class CalculationBasisSetAssociationSchema(ma.ModelSchema):
    basis_set = fields.Nested(
        BasisSetSchema(only=('id', 'family', 'element', ))
        )
    type = fields.Str(attribute='btype')


class CalculationSchema(CalculationListSchema):
    basis_sets = fields.Nested(
        CalculationBasisSetAssociationSchema,
        attribute='basis_set_associations',
        many=True)

    class Meta:
        model = Calculation
        exclude = ('basis_set_associations', )


class CalculationListResource(Resource):
    def get(self):
        schema = CalculationListSchema(many=True)
        return schema.jsonify(Calculation.query.all())

    calculation_args = {
        'collection': fields.Str(
            required=True,
            validate=must_exist_in_db(CalculationCollection, 'name')),
        'test': fields.String(
            required=False, missing=None,
            validate=must_exist_in_db(Test, 'name')),
        'structure': fields.Str(
            required=True, validate=must_exist_in_db(Structure, 'name')),
        'code': fields.String(
            required=True, validate=must_exist_in_db(Code, 'name')),
        'pseudo_family': fields.Str(
            required=True,
            validate=must_exist_in_db(PseudopotentialFamily, 'name')),
        'basis_set_family': fields.Nested({
            'default': fields.Str(required=True,
                                  validate=must_exist_in_db(BasisSetFamily,
                                                            'name')),
            'ri_aux': fields.Str(validate=must_exist_in_db(BasisSetFamily,
                                                           'name')),
            }),
        'restrictions': fields.Dict(missing=None),
        }

    @use_kwargs(calculation_args)
    def post(self, collection, test, structure, code,
             pseudo_family, basis_set_family, restrictions):

        # replace IDs by ORM objects where necessary:

        code = Code.query.filter_by(name=code).one()

        if test:
            test = Test.query.filter_by(name=test).one()

        collection = (CalculationCollection.query
                      .filter_by(name=collection).one())

        structure = (Structure.query
                     .filter_by(name=structure).one())
        atoms = Json2Atoms(structure.ase_structure)
        kinds = set(atoms.get_chemical_symbols())

        pseudos = (PseudopotentialFamily.query
                   .filter_by(name=pseudo_family)
                   .one().pseudos
                   .filter(Pseudopotential.element.in_(kinds))
                   .filter(Pseudopotential.format == code.pseudo_format)
                   .all())

        if len(pseudos) != len(kinds):
            raise ValidationError(("Pseudo Family {}"
                                   " does not cover all kinds required"
                                   ).format(pseudos.family.name))

        all_basis_sets = {}
        for btype, family_name in basis_set_family.items():
            # for each family spec get the corresponding basis_sets
            all_basis_sets[btype] = (BasisSetFamily.query
                                     .filter_by(name=family_name)
                                     .one().basissets
                                     .filter(BasisSet.element.in_(kinds))
                                     .all())

            # check that we got a basis set for each kind
            if len(all_basis_sets[btype]) != len(kinds):
                raise ValidationError(("Basis Set Family {}"
                                       " does not cover all kinds required"
                                       ).format(family_name))

        default_settings = None

        if test and collection:
            default_settings = (db.session
                                .query(CalculationDefaultSettings.settings)
                                .filter_by(code=code, test=test)
                                .scalar())

        if default_settings:
            app.logger.info(("found default settings for code '%s'"
                             " and test '%s'"), code.name, test.name)

        calculation = Calculation(collection=collection, test=test,
                                  structure=structure, code=code,
                                  settings=default_settings,
                                  restrictions=restrictions,
                                  pseudos=pseudos)

        for btype, basis_sets in all_basis_sets.items():
            # create the associations between Calculation and Basis Set
            for basis_set in basis_sets:
                assoc = CalculationBasisSet(btype=btype, basis_set=basis_set)
                calculation.basis_set_associations.append(assoc)

        db.session.add(calculation)
        db.session.commit()
        schema = CalculationSchema()
        return schema.jsonify(calculation)


class CalculationResource(Resource):
    def get(self, cid):
        schema = CalculationSchema()
        return schema.jsonify((Calculation.query.get_or_404(cid)))


class Task2ListSchema(ma.ModelSchema):
    id = fields.UUID()
    status = fields.Str(attribute='status.name')
    ctime = fields.DateTime()
    mtime = fields.DateTime()
    machine = fields.Str(attribute='machine.name')
    priority = fields.Int()

    _links = ma.Hyperlinks({
        'self': ma.URLFor('task2resource', tid='<id>'),
        'collection': ma.URLFor('task2listresource'),
        'uploads': ma.URLFor('task2uploadresource', tid='<id>'),
        })


class Task2Schema(Task2ListSchema):
    calculation = fields.Nested(CalculationListSchema)
    infiles = fields.Nested(ArtifactSchema, many=True)
    outfiles = fields.Nested(ArtifactSchema, many=True)
    data = fields.Dict()

    class Meta:
        model = Task2
        exclude = ('artifacts', )


class Task2ListResource(Resource):
    # calculation is set when called as subcollection of a calculation
    @use_kwargs({
        'limit': fields.Int(missing=None),
        'status': fields.Str(validate=must_exist_in_db(TaskStatus, 'name'),
                             missing=None)
        })
    def get(self, limit, status, calculation=None):
        schema = Task2ListSchema(many=True)
        query = (Task2.query
                 .join(TaskStatus)
                 .options(contains_eager('status')))

        if calculation is not None:
            query = query.filter(Task2.calculation == calculation)

        if status is not None:
            query = query.filter(TaskStatus.name == status)

        if limit is not None:
            query = query.limit(limit)

        return schema.jsonify(query.all())

    @use_kwargs({
        'calculation': fields.UUID(required=True,
                                   validate=must_exist_in_db(Calculation)),
        })
    def post(self, calculation):
        task = Task2(calculation)
        db.session.add(task)
        db.session.commit()
        return Task2Schema().jsonify(task)


class CalculationTask2ListResource(Task2ListResource):
    def post(self, cid):
        task = Task2(cid)
        db.session.add(task)
        db.session.commit()
        return Task2Schema().jsonify(task)


# dict of possible task state transitions:
TASK_STATES = {
    'new': ['pending', 'cancelled', 'deferred'],
    'pending': ['running', 'error'],
    'cancelled': [],
    'deferred': ['new', 'cancelled'],
    'running': ['done', 'error'],
    'done': [],
    'error': [],
    }


def task_args_validate(input_dict):
    # data can only be set when ending a task (in whatever way)
    if ('data' in input_dict.keys() and
            input_dict['status'] not in ['done', 'error']):
        return False

    # machine can only be set when going into pending
    if ('machine' in input_dict.keys() and
            input_dict['status'] not in ['pending']):
        return False

    # when setting a task to pending we need a valid machine
    # to generate the correct runscript
    if input_dict['status'] == 'pending':
        return ('machine' in input_dict.keys() and
                (Machine.query
                 .filter_by(shortname=input_dict['machine'])
                 .one_or_none() is not None))

    return True


class Task2Resource(Resource):
    def get(self, tid):
        schema = Task2Schema()
        return schema.jsonify((Task2.query.get_or_404(tid)))

    @use_kwargs({
        'status': fields.Str(required=True,
                             validate=lambda k: k in TASK_STATES.keys()),
        'data': fields.Dict(),
        'machine': fields.Str(),
        }, validate=task_args_validate)
    def patch(self, tid, status, data, machine):
        task = Task2.query.get_or_404(tid)

        if status not in TASK_STATES[task.status.name]:
            raise ValidationError("can not set task to {} from {}".format(
                status, task.status.name))

        calc = task.calculation

        # setting a task to pending will generate the input files
        # since this is the point where the machine is finally known
        if status == 'pending':
            task.machine = Machine.query.filter_by(shortname=machine).one()

            # see that we have a command for this code and machine
            command = (Command.query
                       .filter_by(machine=task.machine, code=calc.code)
                       .one())

            inputs = {}
            # we put all files in one task-subfolder
            # to avoid hitting a max-file-per-dir limit
            basepath = "fkup://results/{t.id}/".format(t=task)

            if calc.code.name == 'cp2k':

                generated_input = {
                    'global': {
                        'project': "fatman.calc",
                        },
                    'force_eval': {
                        'dft': {
                            'basis_set_file_name': "./BASIS_SETS",
                            'potential_file_name': "./POTENTIALS",
                            'poisson': {
                                'periodic': None,
                                },
                            },
                        'subsys': {
                            'cell': {  # filled out below
                                'a': None,  # filled out below
                                'b': None,
                                'c': None,
                                'periodic': None,
                                },
                            'topology': {
                                'coord_file': "./struct.xyz",
                                'coord_file_format': 'XYZ',
                                },
                            'kind': [],  # filled out below
                            },
                        },
                    }

                inputs['struct'] = Artifact(name="struct.xyz",
                                            path=basepath+"{id}")

                # TODO: simplify the following:
                struct = Json2Atoms(calc.structure.ase_structure)

                if any(struct.get_pbc()) and not all(struct.get_pbc()):
                    app.logger.error("mixed-periodic not (yet) supported")
                    abort(500)

                periodic = "XYZ" if any(struct.get_pbc()) else "NONE"

                generated_input['force_eval']['dft']['poisson']['periodic'] \
                    = periodic

                cell = {
                    'a': (('[angstrom]',) + tuple(struct.get_cell()[:, 0])),
                    'b': (('[angstrom]',) + tuple(struct.get_cell()[:, 1])),
                    'c': (('[angstrom]',) + tuple(struct.get_cell()[:, 2])),
                    'periodic': periodic,
                    }
                generated_input['force_eval']['subsys']['cell'] = cell

                if 'key_value_pairs' in struct.info:
                    # it seems Python ASE is unable to handle nested dicts in
                    # the Atoms.info attribute when writing XYZ, even though it
                    # creates it in the first place
                    # see https://gitlab.com/ase/ase/issues/60
                    struct.info = dict(mergedicts(
                        {k: v for k, v in struct.info.items()
                         if k != 'key_value_pairs'},
                        struct.info['key_value_pairs']))

                bytebuf = BytesIO()
                stringbuf = TextIOWrapper(bytebuf, write_through=True)
                ase_io.write(stringbuf, struct, format='xyz')
                bytebuf.seek(0)
                inputs['struct'].save(bytebuf)

                inputs['basis'] = Artifact(name="BASIS_SETS",
                                           path=basepath+"{id}")
                bytebuf = BytesIO()
                bytebuf.write(("# {a.name}:"
                               " AUTOGENERATED by FATMAN for Task {t.id}\n")
                              .format(a=inputs['basis'], t=task)
                              .encode('utf-8'))

                # for the basis sets we have to be able to
                # lookup the entry by element
                kind = {s: {'_': s,
                            'basis_set': [],
                            'potential': None}
                        for s in struct.get_chemical_symbols()}

                # the total number of MOs to calculate the required MOs for smearing
                n_mos = 0

                for basis_set_assoc in calc.basis_set_associations:
                    element = basis_set_assoc.basis_set.element

                    kind[element]['basis_set'].append((
                        'ORB' if basis_set_assoc.btype == 'default'
                        else basis_set_assoc.btype.upper(),
                        basis_set_assoc.basis_set.family.name))

                # TODO: in case the same basis set is used for different types
                #       we write it twice
                for basis_set in set(calc.basis_sets):
                    bytebuf.write(("# Basis Set ID {b.id}\n"
                                   "{b.element} {b.family.name}\n"
                                   "{b.basis}\n")
                                  .format(b=basis_set)
                                  .encode('utf-8'))

                    # the number of MOs depends on the basis set
                    econfig_string = basis_set.basis.split('\n')[1]
                    econfig = [int(n) for n in econfig_string.split()]
                    # sum over (the number of m's per l quantum number times
                    # the number of functions per m):
                    n_mos += np.dot(
                        [2*l+1 for l in range(econfig[1], econfig[2]+1)],
                        econfig[4:])

                bytebuf.seek(0)
                inputs['basis'].save(bytebuf)

                inputs['pseudos'] = Artifact(name="POTENTIALS",
                                             path=basepath+"{id}")
                bytebuf = BytesIO()
                bytebuf.write(("# {a.name}:"
                               " AUTOGENERATED by FATMAN for Task {t.id}\n")
                              .format(a=inputs['pseudos'], t=task)
                              .encode('utf-8'))
                for pseudo in calc.pseudos:
                    kind[pseudo.element]['potential'] = (
                        "{p.family.name}-q{p.core_electrons}"
                        .format(p=pseudo))
                    # the format is checked when creating the Calculation
                    bytebuf.write(("# Pseudopotential ID {p.id}\n"
                                   "{p.element}"
                                   " {p.family.name}-q{p.core_electrons}"
                                   " {p.family.name}\n"
                                   "{p.pseudo}\n")
                                  .format(p=pseudo)
                                  .encode('utf-8'))
                bytebuf.seek(0)
                inputs['pseudos'].save(bytebuf)

                # in the CP2K Python dict struct the kinds are stored as list
                generated_input['force_eval']['subsys']['kind'] = list(
                        kind.values())

                # we _want_ to bail out here if there is no input
                user_input = task.calculation.settings['input']
                combined_input = dict(mergedicts(user_input, generated_input))

                try:
                    # if scf is itself a dict we get a reference here
                    scf = combined_input['force_eval']['dft']['scf']
                    if 'smear' in scf.keys() and 'added_mos' not in scf.keys():
                        scf['added_mos'] = int(0.3*n_mos)
                except KeyError:
                    pass

                inputs['calc'] = Artifact(name="calc.inp",
                                          path=basepath+"{id}")
                bytebuf = BytesIO()
                bytebuf.write(("# {a.name}:"
                               " AUTOGENERATED by FATMAN for task {t.id}\n")
                              .format(a=inputs['calc'], t=task)
                              .encode('utf-8'))

                dict2cp2k(combined_input, bytebuf, parameters=struct.info)
                bytebuf.seek(0)
                inputs['calc'].save(bytebuf)

            else:
                app.logger.error("code {} not (yet) supported", calc.code.name)
                abort(500)

            for artifact in inputs.values():
                db.session.add(Task2Artifact(artifact=artifact, task=task,
                                             linktype="input"))

            # intentionally merge the task settings over the other settings
            # to give the user the possibility to specify settings at task
            # creation which take precedence over any other settings

            settings = {
                'machine': task.machine.settings,
                'environment': (
                    command.environment if command.environment else {}),
                'commands': command.commands,
                }

            task.settings = dict(mergedicts(
                settings,
                task.settings if task.settings else {}))

        elif status in ['error', 'new']:
            task.data = data

        task.status = TaskStatus.query.filter_by(name=status).one()
        db.session.commit()

        schema = Task2Schema()
        return schema.jsonify(task)


class Task2UploadResource(Resource):
    upload_args = {
        'name': fields.Str(required=True),
        }
    file_args = {
        'data': fields.Field(required=True)
        }

    @use_kwargs(upload_args)
    @use_kwargs(file_args, locations=('files', ))
    def post(self, tid, name, data):
        task = Task2.query.get_or_404(tid)

        basepath = "fkup://results/{t.id}/".format(t=task)
        artifact = Artifact(name=name, path=basepath+"{id}")
        artifact.save(data)

        db.session.add(Task2Artifact(artifact=artifact, task=task,
                                     linktype="output"))
        db.session.commit()
        schema = ArtifactSchema()
        return schema.jsonify(artifact)


class StructureSchema(ma.ModelSchema):
    calculations = fields.Nested(CalculationListSchema, many=True,
                                 exclude=('structure', ))

    sets = fields.DelimitedList(fields.Str())
    _links = ma.Hyperlinks({
        'self': ma.URLFor('structureresource_v2', sid='<id>'),
        'collection': ma.URLFor('structurelistresource_v2'),
        })

    class Meta:
        model = Structure
        # tests is an old table
        exclude = ('tests', )


class StructureListResource_v2(Resource):
    def get(self):
        schema = StructureSchema(many=True,
                                 exclude=('calculations', 'ase_structure',))
        return schema.jsonify(Structure.query.all())

    structure_args = {
        'name': fields.Str(required=True),
        'sets': fields.DelimitedList(fields.Str(
            required=True,
            validate=must_exist_in_db(StructureSet, 'name'))),
        'pbc': fields.Boolean(required=False, default=True),
        'gformat': fields.Str(required=True, load_from='format')
        }
    file_args = {
        'geometry': fields.Field(required=True)
        }

    @use_kwargs(structure_args)
    @use_kwargs(file_args, locations=('files', ))
    def post(self, name, sets, pbc, gformat, geometry):
        sets = (StructureSet.query
                .filter(StructureSet.name.in_(sets))
                .all())

        struct = ase_io.read(TextIOWrapper(geometry),
                             format=gformat)

        struct.set_pbc(pbc)
        pos = struct.get_positions()

        # if the structure is passed in XYZ-format, we generate
        # a cell which is 2x the size of the molecule (3x if non-periodic)
        # but at least 5 Å (important for 2D structures)
        scale = 2. if pbc else 3.
        min_d = 5.
        cell = [max(min_d, scale*(max(pos[:, 0]) - min(pos[:, 0]))),
                max(min_d, scale*(max(pos[:, 1]) - min(pos[:, 1]))),
                max(min_d, scale*(max(pos[:, 2]) - min(pos[:, 2])))]
        # and round up
        cell = np.ceil(cell)

        struct.set_cell(cell, scale_atoms=False)
        ase_structure = Atoms2Json(struct)

        structure = Structure(name=name, sets=sets,
                              ase_structure=ase_structure)
        db.session.add(structure)
        db.session.commit()
        schema = StructureSchema()
        return schema.jsonify(structure)


class StructureResource_v2(Resource):
    def get(self, sid):
        schema = StructureSchema()
        return schema.jsonify((Structure.query.get_or_404(sid)))


# This error handler is necessary for usage with Flask-RESTful
@parser.error_handler
def handle_request_parsing_error(err):
    """webargs error handler that uses Flask-RESTful's abort function to return
    a JSON error response to the client.
    """
    abort(422, errors=err.messages)


api = Api(app, prefix='/api/v2')

api.add_resource(CalculationCollectionListResource, '/calculationcollections')
api.add_resource(CalculationCollectionResource,
                 '/calculationcollections/<uuid:ccid>')
api.add_resource(CalculationListResource, '/calculations')
api.add_resource(CalculationResource, '/calculations/<uuid:cid>')
api.add_resource(CalculationTask2ListResource,
                 '/calculations/<uuid:cid>/tasks')
api.add_resource(Task2ListResource, '/tasks')
api.add_resource(Task2Resource, '/tasks/<uuid:tid>')
api.add_resource(Task2UploadResource, '/tasks/<uuid:tid>/uploads')
api.add_resource(ArtifactListResource, '/artifacts')
api.add_resource(ArtifactResource, '/artifacts/<uuid:aid>')
api.add_resource(ArtifactDownloadResource, '/artifacts/<uuid:aid>/download')
api.add_resource(StructureListResource_v2, '/structures')
api.add_resource(StructureResource_v2, '/structures/<uuid:sid>')
