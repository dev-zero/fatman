
import bz2
from urllib.parse import urlsplit
from os.path import basename
from io import TextIOWrapper, BytesIO
import copy

import flask
from flask import make_response
from flask_marshmallow import Marshmallow
from flask_restful import Api, Resource
from webargs import fields, ValidationError
from webargs.flaskparser import (
    use_kwargs,
    parser,
    abort,
    )
from werkzeug.exceptions import HTTPException
from sqlalchemy.orm import contains_eager, joinedload

from ase import io as ase_io
import numpy as np

from . import app, db, resultfiles, apiauth
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
from .tools.slurm import generate_slurm_batch_script
from .tasks import (
    generate_calculation_results,
    generate_test_result,
    )

ma = Marshmallow(app)


def must_exist_in_db(model, field='id'):
    def check(model_id):
        if not model.query.filter_by(**{field: model_id}).first():
            raise ValidationError('The specified {} does not exist'.format(
                model.__name__))

    return check


class ArtifactSchema(ma.Schema):
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('artifactresource', aid='<id>'),
        'collection': ma.AbsoluteURLFor('artifactlistresource'),
        'download': ma.AbsoluteURLFor('artifactdownloadresource', aid='<id>'),
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


class BasisSetSchema(ma.ModelSchema):
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('basissetresource', bid='<id>'),
        'collection': ma.AbsoluteURLFor('basissetlistresource'),
        })

    family = fields.Str(attribute='family.name')

    class Meta:
        model = BasisSet
        exclude = ('calculations', 'calculation_associations', )


class BasisSetListResource(Resource):
    def get(self):
        schema = BasisSetSchema(exclude=('basis', ), many=True)
        return schema.jsonify(BasisSet.query.all())

    basisset_args = {
        'element': fields.String(required=True),
        'family': fields.Str(
            required=True,
            validate=must_exist_in_db(BasisSetFamily, 'name')),
        }
    file_args = {
        'basis': fields.Field(required=True)
        }

    @apiauth.login_required
    @use_kwargs(basisset_args)
    @use_kwargs(file_args, locations=('files', ))
    def post(self, element, family, basis):
        family = (BasisSetFamily.query
                  .filter_by(name=family)
                  .one())

        basisset = BasisSet(
            element=element,
            family=family,
            basis=TextIOWrapper(basis, encoding='utf-8').read())

        db.session.add(basisset)
        db.session.commit()
        schema = BasisSetSchema()
        return schema.jsonify(basisset)


class BasisSetResource(Resource):
    def get(self, bid):
        schema = BasisSetSchema()
        return schema.jsonify((BasisSet.query.get_or_404(bid)))


class CalculationCollectionSchema(ma.ModelSchema):
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('calculationcollectionresource', ccid='<id>'),
        'collection': ma.AbsoluteURLFor('calculationcollectionlistresource'),
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
        'self': ma.AbsoluteURLFor('calculationresource', cid='<id>'),
        'collection': ma.AbsoluteURLFor('calculationlistresource'),
        'tasks': ma.AbsoluteURLFor('calculationtask2listresource', cid='<id>'),
        })


class CalculationBasisSetAssociationSchema(ma.ModelSchema):
    basis_set = fields.Nested(
        BasisSetSchema(exclude=('basis', ))
        )
    type = fields.Str(attribute='btype')


class CalculationSchema(CalculationListSchema):
    basis_sets = fields.Nested(
        CalculationBasisSetAssociationSchema,
        attribute='basis_set_associations',
        many=True)

    tasks = fields.Nested('Task2ListSchema', many=True)

    class Meta:
        model = Calculation
        exclude = ('basis_set_associations', 'tasks_query', 'testresults_query', )


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
        'basis_set_family_fallback': fields.Nested({
            'default': fields.Str(validate=must_exist_in_db(BasisSetFamily,
                                                            'name')),
            'ri_aux': fields.Str(validate=must_exist_in_db(BasisSetFamily,
                                                           'name')),
            },
            required=False),
        'restrictions': fields.Dict(missing=None),
        }

    @staticmethod
    def new_calculation(collection, test, structure, code,
                        pseudo_family, basis_set_family, basis_set_family_fallback,
                        restrictions):

        # replace IDs by ORM objects where necessary:

        code = Code.query.filter_by(name=code).one()

        if test:
            test = Test.query.filter_by(name=test).one()

        collection = (CalculationCollection.query
                      .filter_by(name=collection)
                      .one())

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
                contained_kinds = set(b.element for b in all_basis_sets[btype])
                missing_kinds = kinds - contained_kinds

                # if there is no fallback at all
                # or no fallback for the same kind of basis set types (ORB, RI, etc.)
                # we are done here:
                if not basis_set_family_fallback or btype not in basis_set_family_fallback.keys():
                    raise ValidationError(("Basis Set Family {} does not contain basis sets for {}"
                                           ).format(family_name, ', '.join(missing_kinds)))

                fallback_basis_sets = (BasisSetFamily.query
                                       .filter_by(name=basis_set_family_fallback[btype])
                                       .one().basissets
                                       .filter(BasisSet.element.in_(missing_kinds))
                                       .all())

                if len(fallback_basis_sets) != len(missing_kinds):
                    missing_kinds = missing_kinds - set(b.element for b in fallback_basis_sets)

                    raise ValidationError(("Neither Basis Set Family {} nor {} contain basis sets for {}"
                                           ).format(family_name,
                                                    basis_set_family_fallback[btype],
                                                    ', '.join(missing_kinds)))

                # if we found the missing basis sets in the fallback family, add them to the list
                app.logger.info(
                    "found basis sets for %s for structure %s in fallback family", ', '.join(missing_kinds),
                    structure.name)
                all_basis_sets[btype] += fallback_basis_sets

        default_settings = None

        if test and code:
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

        return calculation

    @apiauth.login_required
    @use_kwargs(calculation_args)
    def post(self, **kwargs):
        calculation = self.new_calculation(**kwargs)
        db.session.add(calculation)
        db.session.commit()
        schema = CalculationSchema()
        return schema.jsonify(calculation)


class CalculationResource(Resource):
    def get(self, cid):
        calculation = (Calculation.query
                       .options(joinedload('tasks').joinedload('machine'))
                       .options(joinedload('basis_sets').load_only("id", "element"))
                       .options(joinedload('basis_set_associations'))
                       .options(joinedload('pseudos').load_only("id", "element"))
                       .options(joinedload('collection'))
                       .options(joinedload('code'))
                       .options(joinedload('structure').load_only("id", "name"))
                       .options(joinedload('test'))
                       .options(joinedload('testresults'))
                       .get_or_404(cid))
        schema = CalculationSchema()
        return schema.jsonify(calculation)


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


class Task2ListResource(Resource):
    # calculation is set when called as subcollection of a calculation
    @use_kwargs({
        'limit': fields.Int(missing=None),
        'machine': fields.Str(missing=None),
        'status': fields.DelimitedList(
            fields.Str(validate=must_exist_in_db(TaskStatus, 'name'),
                       missing=None)),
        })
    def get(self, limit, machine, status, calculation=None):
        schema = Task2ListSchema(many=True)
        query = (Task2.query
                 .join(TaskStatus)
                 .options(
                     contains_eager('status'),
                     ))

        if calculation is not None:
            query = query.filter(Task2.calculation == calculation)

        if status:
            query = query.filter(TaskStatus.name.in_(status))

        if machine is not None:
            query = query.filter(Task2.machine.has(shortname=machine))

        if limit is not None:
            query = query.limit(limit)

        return schema.jsonify(query.all())

    @apiauth.login_required
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
    @apiauth.login_required
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

    @apiauth.login_required
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

            if calc.code.name == 'CP2K':

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

            runner = task.machine.settings['runner']

            # ensure that we don't accidentally modify ORM objects
            commands = copy.deepcopy(command.commands)
            environment = copy.deepcopy(command.environment) if command.environment else {}
            machine_settings = copy.deepcopy(task.machine.settings)

            # if the command defines arguments for this runner, merge them over the machine settings
            # The machine has only one set of runner_args, namely for the runner it is supposed to run
            if runner in calc.settings.get('runner_args', {}).keys():
                machine_settings = dict(mergedicts(
                    machine_settings,
                    {'runner_args': calc.settings['runner_args'][runner]}
                    ))

            # merge command_args manually since they are not a dict (to preserve order)
            for name, args in calc.settings.get('command_args', {}).items():
                for cmd in commands:
                    if cmd['name'] == name:
                        cmd['args'] += args
                        break

            # if the user specifies a new modules list, the one from command will get overridden
            environment = dict(mergedicts(
                environment,
                calc.settings.get('command_environment', {})))

            # If the command defines arguments (currently only runner_args) for this machine,
            # merge them over the runner-specific settings.
            # This makes it possible to set an estimated runtime from the calculation object.
            if machine in calc.settings.get('machine_settings', {}).keys():
                machine_settings = dict(mergedicts(
                    machine_settings,
                    {k: v for k, v in calc.settings['machine_settings'][machine].items() if k != 'command_args'}
                    ))

                # merge command_args manually since they are not a dict (to preserve order)
                for name, args in calc.settings['machine_settings'][machine].get('command_args', {}).items():
                    for cmd in commands:
                        if cmd['name'] == name:
                            cmd['args'] += args
                            break

            # intentionally merge the task settings over the other settings
            # to give the user the possibility to specify settings at task
            # creation which take precedence over any other settings
            settings = {
                # define a task name usable on most OS and with chars directly usable in URLs
                'name': 'fatman.{}'.format(task.id),
                'machine': machine_settings,
                # we always export environment and commands for easier introspection, even
                # though certain runners already contain them in their batch script
                'environment': environment,
                'commands': commands,
                'output_artifacts': calc.settings['output_artifacts'],
                }

            task.settings = dict(mergedicts(
                settings,
                task.settings if task.settings else {}))

            # This is after the settings merging by intention and uses directly merged task values
            # A client could in principal generate this file instead based on the exported data,
            # but we decided to do it on the server for archival purposes.
            if runner == "slurm":
                inputs['runner'] = Artifact(name="run.sh",
                                            path=basepath+"{id}")

                bytebuf = BytesIO()
                generate_slurm_batch_script(
                    name=task.settings['name'],
                    commands=task.settings['commands'],
                    environment=task.settings['environment'],
                    sbatch_args=task.settings['machine'].get('runner_args', {}).get('sbatch'),
                    srun_args=task.settings['machine'].get('runner_args', {}).get('srun'),
                    output=bytebuf)
                bytebuf.seek(0)
                inputs['runner'].save(bytebuf)

            elif task.settings['machine']['runner'] == "direct":
                # a direct runner uses the commands provided and runs them directly (hence the name)
                pass
            elif task.settings['machine']['runner'] == "mpirun":
                # a mpi runner uses the commands provided and wraps them in mpirun
                pass
            else:
                app.logger.error("runner {} not (yet) supported", task.settings['machine']['runner'])
                abort(500)

            # now that we have all artifacts in place, add them to the task
            for artifact in inputs.values():
                db.session.add(Task2Artifact(artifact=artifact, task=task,
                                             linktype="input"))

        elif status in ['error', 'new', 'cancelled', 'running', 'deferred', 'done']:
            if data:
                if task.data:
                    task.data = mergedicts(task.data, data, True)
                else:
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

    @apiauth.login_required
    @use_kwargs(upload_args)
    @use_kwargs(file_args, locations=('files', ))
    def post(self, tid, name, data):
        task = (Task2.query
                .options(joinedload('calculation'))
                .get_or_404(tid))

        basepath = "fkup://results/{t.id}/".format(t=task)
        artifact = Artifact(name=name, path=basepath+"{id}")
        artifact.save(data)

        db.session.add(Task2Artifact(artifact=artifact, task=task,
                                     linktype="output"))
        db.session.commit()

        # start generating results and test results if succeeded
        generate_calculation_results.apply_async(
            (task.calculation.id,),
            link=generate_test_result.si(task.calculation.id))

        schema = ArtifactSchema()
        return schema.jsonify(artifact)


class StructureSchema(ma.ModelSchema):
    calculations = fields.Nested(CalculationListSchema, many=True,
                                 exclude=('structure', ))

    sets = fields.DelimitedList(fields.Str())
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('structureresource_v2', sid='<id>'),
        'collection': ma.AbsoluteURLFor('structurelistresource_v2'),
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

    @apiauth.login_required
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


class StructureSetSchema(ma.ModelSchema):
    _links = ma.Hyperlinks({
        'self': ma.AbsoluteURLFor('structuresetresource', name='<name>'),
        'collection': ma.AbsoluteURLFor('structuresetlistresource'),
        'calculations': ma.AbsoluteURLFor('structuresetcalculationslistresource', name='<name>'),
        })

    class Meta:
        model = StructureSet
        exclude = ('id', 'structures', )


class StructureSetListResource(Resource):
    def get(self):
        schema = StructureSetSchema(many=True)
        return schema.jsonify(StructureSet.query.all())


class StructureSetCalculationsListResource(Resource):
    def get(self, name):
        # validate the name
        StructureSet.query.filter_by(name=name).one()

        schema = CalculationListSchema(many=True)
        return schema.jsonify(Calculation.query
                              .join(Structure)
                              .join(StructureSet, Structure.sets)
                              .filter(StructureSet.name == name)
                              .all())

    @apiauth.login_required
    @use_kwargs({k: v for k, v in CalculationListResource.calculation_args.items() if k != 'structure'})
    def post(self, name, **kwargs):
        # validate the name
        StructureSet.query.filter_by(name=name).one()

        structures = (db.session.query(Structure.name)
                      .join(StructureSet, Structure.sets)
                      .filter(StructureSet.name == name)
                      .all())

        calculations = []
        errors = {}
        for structure in structures:
            try:
                calculations.append(CalculationListResource.new_calculation(structure=structure.name, **kwargs))
                db.session.add(calculations[-1])
            except ValidationError as exc:
                app.logger.exception("Creating calculation for structure %s failed", structure.name)
                errors[structure.name] = exc

        if errors:
            # get an original Flask exception and augment it with data
            try:
                flask.abort(422)
            except HTTPException as exc:
                exc.data = {
                    'errors': {
                        'basis_set_family': ['structure {}: {}'.format(s, e) for s, e in errors.items()],
                        },
                    }
                raise exc

        db.session.commit()
        schema = CalculationSchema(many=True)
        return schema.jsonify(calculations)


class StructureSetResource(Resource):
    def get(self, name):
        schema = StructureSetSchema()
        return schema.jsonify(StructureSet.query.filter_by(name=name).one())


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
api.add_resource(StructureSetListResource, '/structuresets')
api.add_resource(StructureSetResource, '/structuresets/<string:name>')
api.add_resource(StructureSetCalculationsListResource, '/structuresets/<string:name>/calculations')
api.add_resource(StructureListResource_v2, '/structures')
api.add_resource(StructureResource_v2, '/structures/<uuid:sid>')
api.add_resource(BasisSetListResource, '/basissets')
api.add_resource(BasisSetResource, '/basissets/<uuid:bid>')
