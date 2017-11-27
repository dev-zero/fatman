
import bz2
from urllib.parse import urlsplit
from os.path import basename
from io import TextIOWrapper, BytesIO, StringIO
import copy
import collections
import itertools

import flask
from flask import make_response, request, send_file, url_for
from flask_restful import Api, Resource
from webargs import fields, ValidationError
from webargs.flaskparser import (
    use_kwargs,
    parser,
    abort,
    )
from werkzeug.wrappers import Response
from werkzeug.exceptions import HTTPException
from sqlalchemy import and_, or_, cast, distinct, literal
from sqlalchemy.orm import contains_eager, joinedload, aliased
from sqlalchemy.dialects.postgresql import JSONB

from ase import io as ase_io, data as ase_data
import numpy as np

from . import app, db, resultfiles, apiauth, capp, calculation_finished
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
    TaskRuntimeSettings,
    Artifact,
    Machine,
    Command,
    TestResult2,
    TestResult2Collection,
    )
from .tools import json2atoms, atoms2json, mergedicts
from .tools.generators import generate_CP2K_inputs
from .tools.slurm import generate_slurm_batch_script
from .tools.webargs import nested_parser
from .tools.deltatest import calcDelta, ATOMIC_ELEMENTS

from .tasks import (
    generate_calculation_results,
    generate_all_calculation_results,
    generate_test_result,
    generate_all_test_results,
    )

from .schemas import (
    ArtifactSchema,
    BasisSetSchema,
    PseudopotentialSchema,
    CalculationCollectionSchema,
    CalculationListSchema,
    TestResultSchema,
    TestResultCollectionSchema,
    CalculationSchema,
    CalculationListActionSchema,
    Task2ListSchema,
    Task2Schema,
    StructureSchema,
    StructureListSchema,
    StructureSetSchema,
    TestResultListActionSchema,
    CodeListSchema,
    CodeCommandListSchema,
    CodeCommandSchema,
    CodeSchema,
    DeltatestComparisonSchema,
    TestListSchema,

    BoolValuedDict,
    )


def must_exist_in_db(model, field='id'):
    def check(model_id):
        if not model.query.filter_by(**{field: model_id}).first():
            raise ValidationError('The specified {} does not exist'.format(
                model.__name__))

    return check


class ArtifactListResource(Resource):
    def get(self):
        schema = ArtifactSchema(many=True)
        return schema.jsonify(Artifact.query.all())


class ArtifactResource(Resource):
    def get(self, aid):
        schema = ArtifactSchema()
        return schema.jsonify(Artifact.query.get_or_404(aid))


class ArtifactDownloadResource(Resource):
    def __init__(self, mode='download'):
        self._content_dispo = 'attachment' if mode == 'download' else 'inline'

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
                    response.headers["Content-Disposition"] = "{}; filename={}".format(
                        self._content_dispo, filename)
                    response.headers["Content-Type"] = "text/plain"
                    return response

            elif compressed is None:
                with open(resultfiles.path(path[1:]), 'rb') as infile:
                    response = make_response(infile.read())
                    response.headers["Content-Disposition"] = "{}; filename={}".format(self._content_dispo, filename)
                    response.headers["Content-Type"] = "text/plain"
                    return response
            else:
                app.logger.error("unknown compression scheme found: %s", compressed)
                abort(500)

        app.logger.warning("{} contains unknown path '{path}'".format(
            artifact, **artifact.__dict__))
        abort(500)


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


class Pseudopotential2ListResource(Resource):
    def get(self):
        schema = PseudopotentialSchema(exclude=('pseudo', ), many=True)
        return schema.jsonify(Pseudopotential.query.all())


class Pseudopotential2Resource(Resource):
    def get(self, pid):
        schema = PseudopotentialSchema()
        return schema.jsonify((Pseudopotential.query.get_or_404(pid)))


class CalculationCollectionListResource(Resource):
    def get(self):
        schema = CalculationCollectionSchema(many=True)
        return schema.jsonify(CalculationCollection.query.all())

    @apiauth.login_required
    @use_kwargs({'name': fields.Str(required=True), 'desc': fields.Str(required=True)})
    def post(self, name, desc):
        coll = CalculationCollection(name=name, desc=desc)

        db.session.add(coll)
        db.session.commit()

        schema = CalculationCollectionSchema()
        return schema.jsonify(coll)


class CalculationCollectionResource(Resource):
    def get(self, ccid):
        schema = CalculationCollectionSchema()
        return schema.jsonify((CalculationCollection.query.get_or_404(ccid)))


class CalculationListResource(Resource):
    calculation_list_args = {
        'collection': fields.Str(
            validate=must_exist_in_db(CalculationCollection, 'name')),
        'test': fields.String(
            validate=must_exist_in_db(Test, 'name')),
        'structure': fields.Str(),  # not validating here since we have a "contains" match for this
        'code': fields.String(
            validate=must_exist_in_db(Code, 'name')),
        'status': fields.String(
            validate=must_exist_in_db(TaskStatus, 'name')),
        'page': fields.Integer(required=False, missing=1, validate=lambda n: n > 0),
        'per_page': fields.Integer(required=False, missing=20, validate=lambda n: n > 0 and n <= 200),
        }

    @use_kwargs(calculation_list_args, locations=('querystring',))
    def get(self, collection, test, structure, code, status, page, per_page):
        schema = CalculationListSchema(many=True)

        # The following join is inspired by http://stackoverflow.com/a/2111420/1400465
        # to first join the last added tasks and then sort the list of Calculations by the mtime
        # of this latest task.
        # Alternatively we could exploit that the ORM is doing a de-dup and sort by mtime right away.
        t2 = aliased(Task2)
        calcs = (db.session.query(Calculation)
                 .options(joinedload('structure', innerjoin=True).load_only('id', 'name'))
                 .options(joinedload('code', innerjoin=True))
                 .options(joinedload('test', innerjoin=True))
                 .join(Task2)
                 .options(contains_eager('tasks').joinedload('machine', innerjoin=True))  # load selected Task together with Calculation
                 .outerjoin(t2, and_(Calculation.id == t2.calculation_id, Task2.ctime < t2.ctime))
                 .filter(t2.id == None))

        if collection:
            calcs = calcs.join(Calculation.collection).filter(CalculationCollection.name == collection)

        if test:
            calcs = calcs.join(Calculation.test).filter(Test.name == test)

        if structure:
            calcs = calcs.join(Calculation.structure).filter(Structure.name.contains(structure))

        if code:
            calcs = calcs.join(Calculation.code).filter(Code.name == code)

        if status:
            calcs = calcs.join(Task2.status).filter(TaskStatus.name == status)

        # can't use Flask-Marshmallow's paginate() here since it will cause a refetch for the complete set
        calc_total_count = calcs.count()
        pages = int(np.ceil(calc_total_count / float(per_page)))

        calcs = (calcs.order_by(Task2.mtime.desc())
                 .limit(per_page)
                 .offset((page - 1)*per_page))

        response = schema.jsonify(calcs)

        link_header = []

        if page < pages:
            link_header.append('<{}>; rel="next"'.format(url_for('calculationlistresource',
                                                                 page=page+1, _external=True)))
        link_header.append('<{}>; rel="last"'.format(url_for('calculationlistresource',
                                                             page=pages, _external=True)))

        if page > 1:
            link_header.append('<{}>; rel="prev"'.format(url_for('calculationlistresource',
                                                                 page=page-1, _external=True)))

        link_header.append('<{}>; rel="first"'.format(url_for('calculationlistresource',
                                                              page=1, _external=True)))

        response.headers['Link'] = ", ".join(link_header)
        response.headers['X-total-count'] = calc_total_count

        return response

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
        'settings': fields.Dict(missing={}),
        'restrictions': fields.Dict(missing=None),
        }

    @staticmethod
    def new_calculation(collection, test, structure, code,
                        pseudo_family, basis_set_family, basis_set_family_fallback,
                        settings, restrictions):

        # replace IDs by ORM objects where necessary:

        code = Code.query.filter_by(name=code).one()

        if test:
            test = Test.query.filter_by(name=test).one()

        collection = (CalculationCollection.query
                      .filter_by(name=collection)
                      .one())

        structure = (Structure.query
                     .filter(Structure.name == structure, Structure.replaced_by_id == None)
                     .one())
        atoms = json2atoms(structure.ase_structure)
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
                                   ).format(pseudo_family))

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
            # get the general per-code and -test settings first
            default_settings = (db.session
                                .query(CalculationDefaultSettings.settings)
                                .filter_by(code=code, test=test, structure=None, structure_set=None)
                                .scalar())

            if not default_settings:
                default_settings = {}
            else:
                app.logger.info("found base default settings for code '%s' and test '%s'", code.name, test.name)

            def walk_ssets(sset, depth=0):
                if depth > 100:
                    raise RuntimeError(
                        "possible circular reference detected in structure super sets for set {}".format(sset.name))

                yield sset.id

                if sset.superset:
                    yield from walk_ssets(sset.superset, depth+1)

            # for each structure set there may be a path of nested sets, follow it upwards
            sset_paths = []
            for sset in structure.sets:
                sset_paths.append([ss for ss in walk_ssets(sset)])

            # for each path, apply per-structure-set settings from top to bottom
            for sset_path in sset_paths:
                for sset_id in reversed(sset_path):
                    additional_settings = (db.session
                                           .query(CalculationDefaultSettings.settings)
                                           .filter_by(code=code, test=test, structure=None, structure_set_id=sset_id)
                                           .scalar())

                    if additional_settings:
                        app.logger.info(("found additional default settings for"
                                         " code '%s', test '%s' and structure set id '%i'"),
                                        code.name, test.name, sset_id)

                        default_settings = dict(mergedicts(default_settings, additional_settings))

            # finally check whether we have structure-specific settings
            additional_settings = (db.session
                                   .query(CalculationDefaultSettings.settings)
                                   .filter_by(code=code, test=test, structure=structure, structure_set=None)
                                   .scalar())

            if additional_settings:
                app.logger.info(("found additional default settings for"
                                 " code '%s', test '%s' and structure '%s'"),
                                code.name, test.name, structure)

                default_settings = dict(mergedicts(default_settings, additional_settings))

        if default_settings:
            # merge the settings specified by the user over the default_settings
            # to give the user the possibility to overwrite them, settings is at least an empty dict
            settings = dict(mergedicts(default_settings, settings))

        calculation = Calculation(collection=collection, test=test,
                                  structure=structure, code=code,
                                  settings=settings,
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


class CalculationListActionResource(Resource):
    @apiauth.login_required
    @use_kwargs(CalculationListActionSchema)
    def post(self, generateResults):
        if generateResults:
            async_result = generate_all_calculation_results.delay(generateResults['update'])

            return Response(status=202, headers={
                'Location': api.url_for(ActionResource, aid=async_result.id, _external=True)})

        abort(400)


class CalculationActionResource(Resource):
    @apiauth.login_required
    @use_kwargs(CalculationListActionSchema)
    def post(self, cid, generateResults):
        if generateResults:
            async_result = generate_calculation_results.delay(cid, generateResults['update'])

            return Response(status=202, headers={
                'Location': api.url_for(ActionResource, aid=async_result.id, _external=True)})

        abort(400)


class CalculationPreviewResource(Resource):
    def get(self, cid):
        calc = Calculation.query.get_or_404(cid)

        bsets = []
        for basis_set_assoc in calc.basis_set_associations:
            bsets.append((
                basis_set_assoc.btype,
                basis_set_assoc.basis_set.id,
                basis_set_assoc.basis_set.element,
                basis_set_assoc.basis_set.family.name,
                basis_set_assoc.basis_set.basis))

        if calc.code.name == "CP2K":
            inputs = generate_CP2K_inputs(
                calc.settings['input'],
                bsets,
                [(p.id, p.element, p.family.name, p.core_electrons, p.pseudo) for p in calc.pseudos],
                json2atoms(calc.structure.ase_structure),
                "AUTOGENERATED by FATMAN for preview",
                )
        else:
            abort(501)

        template = """
******************** BEGIN: {name} ********************
{content}
******************** END:   {name} ********************
"""
        preview = "".join(template.format(name=n, content=b.getvalue().decode('utf-8')) for n, b in inputs.items())
        response = make_response(preview)
        response.headers["Content-Disposition"] = "inline; filename=calculation_{}.preview.txt".format(calc.id)
        response.headers["Content-Type"] = "text/plain"
        return response


class CalculationResource(Resource):
    def get(self, cid):
        calculation = (Calculation.query
                       .options(joinedload('tasks').joinedload('machine'))
                       .options(joinedload('basis_sets').defer('basis'))
                       .options(joinedload('basis_set_associations'))
                       .options(joinedload('pseudos').defer('pseudo'))
                       .options(joinedload('collection'))
                       .options(joinedload('code'))
                       .options(joinedload('structure').load_only("id", "name"))
                       .options(joinedload('test'))
                       .options(joinedload('testresults'))
                       .get_or_404(cid))
        schema = CalculationSchema()
        return schema.jsonify(calculation)


class Task2ListResource(Resource):
    # calculation is set when called as subcollection of a calculation
    @use_kwargs({
        'limit': fields.Int(missing=None),
        'machine': fields.Str(missing=None),
        'status': fields.DelimitedList(
            fields.Str(validate=must_exist_in_db(TaskStatus, 'name'),
                       missing=None)),
        }, locations=('querystring',))
    @use_kwargs({
        'worker_machine': fields.Str(load_from='x-fatman-worker-machine', missing=None),
        }, locations=('headers',))
    def get(self, limit, machine, status, worker_machine, cid=None):
        schema = Task2ListSchema(many=True)
        query = (Task2.query
                 .join(TaskStatus)
                 .options(joinedload('machine').load_only('name'))
                 .options(contains_eager('status')))

        if cid is not None:
            query = query.join(Task2.calculation).filter(Calculation.id == cid)

        if status:
            query = query.filter(TaskStatus.name.in_(status))

        if machine is not None:
            query = query.filter(Task2.machine.has(shortname=machine))

        if worker_machine is not None:  # if we don't know the workers name, assume a browser and don't filter
            query = query.filter(or_(
                Task2.restrictions['machine'] == None,  # unrestricted tasks, or...
                literal(worker_machine).op("~")(Task2.restrictions['machine'].astext) # tasks with matching regex
                ))

        if limit is not None:
            query = query.limit(limit)

        return schema.jsonify(query.all())

    @apiauth.login_required
    @use_kwargs({
        'calculation': fields.UUID(required=False,
                                   validate=must_exist_in_db(Calculation)),
        'status': fields.String(required=False, missing=None,
                                validate=must_exist_in_db(TaskStatus, 'name')),
        'restrictions': fields.Dict(required=False, missing=None),
        })
    def post(self, calculation, status, restrictions, cid=None):
        if calculation:
            cid = calculation

        if not cid:  # an invalid calculation id will be caught by the DB
            app.logger.error("neither cid nor calculation defined for creating a task")
            flask.abort(500)  # raise an ISE since this should not happen

        if status not in ['new', 'deferred']:
            raise ValidationError("can create task only as new or deferred")

        calc = Calculation.query.get_or_404(cid)

        if calc.restrictions:
            if restrictions:  # if both the calc restrictions and additional restrictions are not empty, merge them
                restrictions = dict(mergedicts(calc.restrictions, restrictions))
            else:  # otherwise overwrite them
                restrictions = calc.restrictions

        task = Task2(cid, status, restrictions)
        db.session.add(task)
        db.session.commit()
        return Task2Schema().jsonify(task)


# dict of possible task state transitions:
TASK_STATES = {
    'new': ['pending', 'cancelled', 'deferred'],
    'pending': ['running', 'error'],
    'cancelled': [],
    'deferred': ['new', 'cancelled', 'done'],  # setting a task to done from deferred is for manual upload
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
        return schema.jsonify(Task2.query
                              .options(joinedload('infiles'))
                              .options(joinedload('outfiles'))
                              .options(joinedload('calculation').joinedload('code'))
                              .options(joinedload('calculation').joinedload('test'))
                              .options(joinedload('calculation')
                                       .joinedload('structure')
                                       .load_only('id', 'name'))
                              .get_or_404(tid))

    @apiauth.login_required
    @use_kwargs({
        'status': fields.Str(required=True,
                             validate=lambda k: k in TASK_STATES.keys()),
        'data': fields.Dict(),
        'machine': fields.Str(),
        }, validate=task_args_validate)
    def patch(self, tid, status, data, machine):
        task = (Task2.query
                # lock this task for editing and don't wait: if it is already locked,
                # chances are the client has to reconsider what to do any way since the state will have changed
                .with_for_update(of=Task2, nowait=True)
                .get_or_404(tid))

        # since the task is locked for writing and we can not change state
        # to the same state it already is, we get atomic behaviour here:
        if status not in TASK_STATES[task.status.name]:
            raise ValidationError("can not set task to {} from {}".format(
                status, task.status.name))

        calc = task.calculation

        # setting a task to pending will generate the input files
        # since this is the point where the machine is finally known
        if status == 'pending':
            task.machine = Machine.query.filter_by(shortname=machine).one()

            task_rt_settings = (db.session.query(TaskRuntimeSettings.settings)
                                .filter_by(machine=task.machine)
                                .filter_by(code=calc.code)
                                .filter_by(test=calc.test)
                                .scalar())

            if task_rt_settings is None:
                task_rt_settings = {}
            else:
                app.logger.info("found task runtime settings for machine '%s'", task.machine)

            # see that we have a command for this code and machine
            command = (Command.query
                       .filter_by(machine=task.machine, code=calc.code)
                       .one())

            artifacts = {}
            # we put all files in one task-subfolder
            # to avoid hitting a max-file-per-dir limit
            basepath = "fkup://results/{t.id}/".format(t=task)

            if calc.code.name == 'CP2K':
                bsets = []
                for basis_set_assoc in calc.basis_set_associations:
                    bsets.append((
                        basis_set_assoc.btype,
                        basis_set_assoc.basis_set.id,
                        basis_set_assoc.basis_set.element,
                        basis_set_assoc.basis_set.family.name,
                        basis_set_assoc.basis_set.basis))

                inputs = generate_CP2K_inputs(
                    task.calculation.settings['input'],
                    bsets,
                    [(p.id, p.element, p.family.name, p.core_electrons, p.pseudo) for p in calc.pseudos],
                    json2atoms(calc.structure.ase_structure),
                    "AUTOGENERATED by FATMAN for Task {t.id}".format(t=task),
                    task_rt_settings.get('input', {})
                    )

                for name, bytebuf in inputs.items():
                    artifacts[name] = Artifact(name=name, path=basepath+"{id}")
                    artifacts[name].save(bytebuf)

            else:
                app.logger.error("code {} not (yet) supported", calc.code.name)
                abort(500)

            # ensure that we don't accidentally modify ORM objects
            commands = copy.deepcopy(command.commands)
            environment = copy.deepcopy(command.environment) if command.environment else {}
            machine_settings = copy.deepcopy(task.machine.settings)
            runner = machine_settings['runner']

            # arguments for the runner or the machine can't come from the Calculation object since we don't
            # know the runner or the machine at the point of creation of the Calculation object

            # merge command_args manually since they are not a dict (to preserve order)
            for name, args in calc.settings.get('command_args', {}).items():
                for cmd in commands:
                    if cmd['name'] == name:
                        cmd['args'] += args
                        break

            # if the user specifies a new modules list, the one from command will get overridden.
            # But usually the user will only overwrite environment variables like OMP_NUM_THREADS
            # since everything else depends on the environment which is unknown at the point of
            # creation of the Calculation object
            environment = dict(mergedicts(
                environment,
                calc.settings.get('command_environment', {})))

            # Merge the Task Runtime Settings over the pre-machine selection settings
            # (note: code input merging already happened above)
            environment = dict(mergedicts(
                environment,
                task_rt_settings.get('command', {}).get('environment', {})))

            machine_settings = dict(mergedicts(
                machine_settings,
                task_rt_settings.get('machine', {})))

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
                artifacts['runner.slurm'] = Artifact(name="run.sh", path=basepath+"{id}")

                bytebuf = BytesIO()
                generate_slurm_batch_script(
                    name=task.settings['name'],
                    commands=task.settings['commands'],
                    environment=task.settings['environment'],
                    sbatch_args=task.settings['machine'].get('runner_args', {}).get('sbatch'),
                    srun_args=task.settings['machine'].get('runner_args', {}).get('srun'),
                    output=bytebuf)
                bytebuf.seek(0)
                artifacts['runner.slurm'].save(bytebuf)

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
            for artifact in artifacts.values():
                db.session.add(Task2Artifact(artifact=artifact, task=task,
                                             linktype="input"))

        elif status in ['error', 'new', 'cancelled', 'running', 'deferred', 'done']:
            if data:
                if task.data:
                    task.data = mergedicts(task.data, data)
                else:
                    task.data = data

        task.status = TaskStatus.query.filter_by(name=status).one()
        db.session.commit()

        if status in ['error', 'done']:
            calculation_finished.send(self, task=task, calculation=task.calculation)

        if status == 'done':
            # start generating results and then test results if succeeded
            generate_calculation_results.apply_async(
                (task.calculation.id,),
                link=generate_test_result.si(task.calculation.id))

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

        schema = ArtifactSchema()
        return schema.jsonify(artifact)


class StructureListResource_v2(Resource):
    filter_args = {
        'include_replaced': fields.Boolean(required=False, default=False),
        'limit': fields.Integer(required=False, missing=-1),
        }

    @use_kwargs(filter_args, locations=('query',))
    def get(self, include_replaced, limit):
        query = (Structure.query
                 .join(Structure.sets)
                 .options(contains_eager(Structure.sets))
                )

        if not include_replaced:
            query = query.filter(Structure.replaced_by_id == None)

        if limit > 0:
            query = query.limit(limit)

        schema = StructureListSchema(many=True)
        return schema.jsonify(query.all())

    structure_args = {
        'name': fields.Str(required=True),
        'sets': fields.DelimitedList(fields.Str(required=True, validate=must_exist_in_db(StructureSet, 'name'))),
        'pbc': fields.Boolean(required=False, missing=None),  # assumed to be true when autogenerating cell
        'charges': fields.DelimitedList(fields.Float(), required=False, missing=None),
        'cell': fields.DelimitedList(fields.Float(), required=False, missing=None),
        'magmoms': fields.DelimitedList(fields.Float(), required=False, missing=None),
        'gformat': fields.Str(required=True, load_from='format'),
        'cubic_cell': fields.Boolean(required=False, missing=False),  # whether to generate a cubic cell
        'center': fields.Boolean(required=False, missing=False),  # center coords in the cell (True if cell autogen)
        'replace_existing': fields.Boolean(required=False, missing=False),  # whether to replace an existing structure
        }
    file_args = {
        'geometry': fields.Field(required=True)
        }

    @apiauth.login_required
    @use_kwargs(structure_args)
    @use_kwargs(file_args, locations=('files', ))
    def post(self, name, sets,
             pbc, charges, cell, magmoms, gformat, geometry, cubic_cell, center,
             replace_existing):
        sets = (StructureSet.query
                .filter(StructureSet.name.in_(sets))
                .all())

        existing_structure = (Structure.query
                              .filter(Structure.name == name, Structure.replaced_by_id == None)
                              .one_or_none())

        if existing_structure and not replace_existing:
            try:
                flask.abort(422)
            except HTTPException as exc:
                exc.data = {
                    'errors': {
                        'name': "A structure with this name already exists and replace_existing is not true",
                        },
                    }
                raise exc

        struct = ase_io.read(TextIOWrapper(geometry), format=gformat)

        if pbc is not None:
            struct.set_pbc(pbc)

        if cell is not None:
            struct.set_cell(cell)

        if charges is not None:
            struct.set_initial_charges(charges)

        if magmoms is not None:
            struct.set_initial_magnetic_moments(magmoms)

        if not struct.get_cell().any():
            # if no cell is defined until here (either inside the format or by an explicit cell argument),
            # we generate a cell based on the moleculare boundary box,
            # where the molecule boundary box is defined as the minimal/maximal coordinates
            # over all atoms minus/plus their respective VdW radii plus a buffer of 2.5 Å,
            # resp. 5 for non-pbc on each side
            # So, for the periodic case the buffers "overlap" between images compared to the non-periodic case.
            # This results in each side of the box being > 5 Å for the periodic,
            # respectively > 10 Å for the non-periodic case

            if (pbc is None) or pbc:
                buf_size = 5.
            else:
                buf_size = 10.

            # The VdW radius if available, otherwise the covalent radius
            radii = [ase_data.vdw_radii[n]
                     if not np.isnan(ase_data.vdw_radii[n])
                     else ase_data.covalent_radii[n]
                     for n in struct.get_atomic_numbers()]

            radii_pos = list(zip(radii, struct.get_positions()))

            cell = [
                max(pos[i] + rad for rad, pos in radii_pos)
                - min(pos[i] - rad for rad, pos in radii_pos)
                + buf_size
                for i in range(3)
                ]

            # and round up to nearest full Angstrom
            cell = np.ceil(cell)

            # if a cubic box is requested, use the largest coordinate
            if cubic_cell:
                cell = [max(cell)] * 3

            struct.set_cell(cell, scale_atoms=False)

            # Finally, we center the atom in the cell (with the cell starting at 0/0/0)
            # to also accomodate for non-periodic calculations and cases where we do not want
            # the code to auto-center the coordinates in the box.
            struct.center()

        # also center if explictly requested by the user
        if center:
            struct.center()

        ase_structure = atoms2json(struct)

        structure = Structure(name=name, sets=sets,
                              ase_structure=ase_structure)
        db.session.add(structure)

        if existing_structure:
            existing_structure.replaced_by = structure

        db.session.commit()

        schema = StructureSchema()
        return schema.jsonify(structure)


class StructureResource_v2(Resource):
    def get(self, sid):
        schema = StructureSchema()
        return schema.jsonify(Structure.query.get_or_404(sid))

    @apiauth.login_required
    def delete(self, sid):
        # here we load the replaced relationship to have SQLA set the references to NULL
        structure = (Structure.query
                     .options(joinedload(Structure.replaced))
                     .get_or_404(sid))
        db.session.delete(structure)
        db.session.commit()
        return Response(status=204)  # return completely empty


class StructureDownloadResource(Resource):

    SUPPORTED_MIMETYPES = collections.OrderedDict([
            # mimetype: (file ending, Python ASE formatter)
            ("chemical/x-xyz", ("xyz", "xyz")),
            ("chemical/x-cif", ("cif", "cif")),
            ("chemical/x-pdb", ("pdb", "proteindatabank")),
            ])

    def get(self, sid):
        if not any(mt in request.accept_mimetypes for mt in self.SUPPORTED_MIMETYPES):
            abort(406)

        # this will return XYZ if the client sends */* due to the OrderedDict
        selected_mimetype = request.accept_mimetypes.best_match(self.SUPPORTED_MIMETYPES.keys())
        fileending, formatter = self.SUPPORTED_MIMETYPES[selected_mimetype]

        structure = Structure.query.get_or_404(sid)
        asestruct = json2atoms(structure.ase_structure)

        if 'key_value_pairs' in asestruct.info:
            # it seems Python ASE is unable to handle nested dicts in
            # the Atoms.info attribute when writing XYZ, even though it
            # creates it in the first place
            # see https://gitlab.com/ase/ase/issues/60
            asestruct.info = dict(mergedicts(
                {k: v for k, v in asestruct.info.items()
                 if k != 'key_value_pairs'},
                asestruct.info['key_value_pairs']))

        stringbuf = StringIO()

        ase_io.write(stringbuf, asestruct, format=formatter)

        return send_file(
            BytesIO(stringbuf.getvalue().encode('utf-8')),
            as_attachment=True,
            attachment_filename="{}.{}".format(structure.name, fileending),
            mimetype=selected_mimetype)


class StructureSetListResource(Resource):
    def get(self):
        schema = StructureSetSchema(many=True)
        return schema.jsonify(StructureSet.query.all())

    structure_set_args = {
        'name': fields.String(required=True),
        'description': fields.Str(required=False, missing=None),
        'superset': fields.Str(required=False, missing=None)
        }
    @apiauth.login_required
    @use_kwargs(structure_set_args)
    def post(self, name, description, superset):

        superset_id = None
        if superset:
            superset_id = db.session.query(StructureSet.id).filter_by(name=superset).scalar()
            if not superset_id:
                raise ValidationError("invalid superset '{}' specified".format(superset))

        structureset = StructureSet(name=name, description=description, superset_id=superset_id)
        db.session.add(structureset)
        db.session.commit()

        schema = StructureSetSchema()
        return schema.jsonify(structureset)


class StructureSetStructureListResource(Resource):
    def get(self, name):
        # validate the name
        StructureSet.query.filter_by(name=name).one()

        query = (Structure.query
                 .join(Structure.sets)
                 .options(contains_eager(Structure.sets))
                 .filter(StructureSet.name == name)
                 .filter(Structure.replaced_by_id == None))

        schema = StructureListSchema(many=True)
        return schema.jsonify(query.all())


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

    bulk_calculation_args = {
        k: v for k, v in CalculationListResource.calculation_args.items()
        if k != 'structure'
        }
    bulk_calculation_args['ignore_failed'] = fields.Boolean(required=False, default=False)

    @apiauth.login_required
    @use_kwargs(bulk_calculation_args)
    def post(self, name, ignore_failed, **kwargs):

        # validate the name
        sset = StructureSet.query.filter_by(name=name).one()

        # TODO: the following should be done using a Recursive CTE,
        #        but since our dataset is small enough, oh, well ;-)

        def get_set_ids(structure_set):
            yield structure_set.id

            for subset in structure_set.subsets:
                yield from get_set_ids(subset)

        sset_ids = set(i for i in get_set_ids(sset))

        # get the name of structures already calculated in this calculation collection
        calculated_structures = (db.session.query(Structure.name)
                                 .join(Structure.calculations)
                                 .join(Calculation.collection)
                                 .filter(CalculationCollection.name == kwargs['collection'])
                                 .all())

        structures = (db.session.query(distinct(Structure.name))
                      .join(StructureSet, Structure.sets)
                      .filter(Structure.sets.any(StructureSet.id.in_(sset_ids)),  # get structures in subsets
                              Structure.replaced_by_id == None,  # filter out structures replaced by newer revisions
                              ~Structure.name.in_(calculated_structures))  # ignore structures already calculated
                      .all())

        calculations = []
        errors = {}
        for structure in structures:
            try:
                calculations.append(CalculationListResource.new_calculation(structure=structure, **kwargs))
                db.session.add(calculations[-1])
            except ValidationError as exc:
                app.logger.exception("Creating calculation for structure %s failed", structure)
                errors[structure] = exc

        if errors and not ignore_failed:
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


class TestResultListResource(Resource):
    filter_args = {
        'test': fields.String(required=False),
        'structure': fields.String(required=False),
        # the nested parser unpacks data.foo=10 to data: { 'foo': 10 },
        # TODO: validate checks to be a dict of (str, bool)
        'data': fields.Nested({
            'checks': BoolValuedDict(),
            'coefficients': fields.Dict(),
            'status': fields.String(),
            'element': fields.String(),
            }),
        }

    @nested_parser.use_kwargs(filter_args, locations=('query',))
    def get(self, test, structure, data):
        schema_excluded_columns = ['data', ]

        # optimize the query by eager_loading calculations, structure and task
        query = (TestResult2.query
                 .options(
                     joinedload('calculations').
                         joinedload('structure')
                     )
                 .options(
                     joinedload('calculations').
                         joinedload('tasks')
                     )
                 .options(
                     joinedload('calculations').
                         joinedload('code')
                     )
                 )

        if test:
            query = query.join(TestResult2.test).filter(Test.name.contains(test))

        if structure:
            query = (query
                     .join(TestResult2.calculations)
                     .join(Calculation.structure)
                     .filter(Structure.name.contains(structure)))

        if data:
            # show the data when filtering by data
            schema_excluded_columns.remove('data')

            # here we are going to build the JSON query
            if 'checks' in data:
                # the checks are a series of key: bool
                for check, value in data['checks'].items():
                    query = query.filter(TestResult2.data[('checks', check, )] == cast(value, JSONB))

            if 'element' in data.keys():
                query = query.filter(TestResult2.data['element'] == cast(data['element'], JSONB))

            if 'coefficients' in data.keys():
                for coeff, value in data['coefficients'].items():
                    # we can either have coefficients.R0=10. for an equal match
                    # or something like coefficients.R0.lt=10. for a '< 10' match
                    if isinstance(value, collections.Mapping):
                        op, value = value.popitem()
                    else:
                        op = 'eq'

                    try:
                        value = float(value)
                    except ValueError:
                        raise ValidationError("invalid boolean value specified for 'coefficients.{}'".format(coeff))

                    # Note: PostgreSQL is faster when doing the op evaluation on the JSONB side rather than
                    #       explicitly converting the retrieved values to float and doing the comparison by itself.
                    if op == 'lt':
                        query = query.filter(TestResult2.data[('coefficients', coeff, )] < cast(value, JSONB))
                    elif op == 'le':
                        query = query.filter(TestResult2.data[('coefficients', coeff, )] <= cast(value, JSONB))
                    elif op == 'eq':
                        query = query.filter(TestResult2.data[('coefficients', coeff, )] == cast(value, JSONB))
                    elif op == 'gt':
                        query = query.filter(TestResult2.data[('coefficients', coeff, )] > cast(value, JSONB))
                    elif op == 'ge':
                        query = query.filter(TestResult2.data[('coefficients', coeff, )] >= cast(value, JSONB))
                    else:
                        raise ValidationError("invalid operator '{}' specified for 'coefficients.{}'".format(op, coeff))

        schema = TestResultSchema(many=True, exclude=schema_excluded_columns)
        return schema.jsonify(query.all())


class TestResultResource(Resource):
    def get(self, trid):
        schema = TestResultSchema()
        tresult = (TestResult2.query
                    .options(
                        joinedload('calculations').
                            joinedload('structure', 'tasks')
                        )
                    .get_or_404(trid))

        return schema.jsonify(tresult)


class TestResultListActionResource(Resource):
    @apiauth.login_required
    @use_kwargs(TestResultListActionSchema)
    def post(self, generate):
        if generate:
            async_result = generate_all_test_results.delay(generate['update'])

            return Response(status=202, headers={
                'Location': api.url_for(ActionResource, aid=async_result.id, _external=True)})

        abort(400)


class TestResultCollectionListResource(Resource):
    def get(self):
        query = TestResult2Collection.query

        schema = TestResultCollectionSchema(many=True, exclude=('testresults', ))
        return schema.jsonify(query.all())


    creation_args = {
        'name': fields.String(required=True),
        'desc': fields.String(required=True),
        'testresults': fields.List(fields.UUID(validate=must_exist_in_db(TestResult2)),
                                   required=True),
        }

    @apiauth.login_required
    @use_kwargs(creation_args)
    def post(self, name, desc, testresults):
        trcoll = TestResult2Collection(
            name=name,
            desc=desc,
            # validation already happened, each testresult (id) yields a db model entry:
            testresults=TestResult2.query.filter(TestResult2.id.in_(testresults)).all()
            )

        db.session.add(trcoll)
        db.session.commit()

        schema = TestResultCollectionSchema()
        return schema.jsonify(trcoll)


class TestResultCollectionResource(Resource):
    def get(self, trcid):
        trc = TestResult2Collection.query.get(trcid)

        schema = TestResultCollectionSchema()
        return schema.jsonify(trc)


    @apiauth.login_required
    def delete(self, trcid):
        trcoll = TestResult2Collection.query.get_or_404(trcid)
        db.session.delete(trcoll)
        db.session.commit()
        return Response(status=204)  # return completely empty


class ActionResource(Resource):
    def get(self, aid):
        async_result = capp.AsyncResult(str(aid))

        if not async_result.ready():
            return Response(status=202, headers={'Location': request.path})

        # TODO: figure out task name using Celery 4.x and add more infos
        return {'status': async_result.status}


class CodeListResource(Resource):
    def get(self):
        schema = CodeListSchema(many=True)
        return schema.jsonify(Code.query.all())


class CodeResource(Resource):
    def get(self, cid):
        schema = CodeSchema()
        return schema.jsonify((Code.query.get_or_404(cid)))


class CodeCommandListResource(Resource):
    def get(self, cid):
        schema = CodeCommandListSchema(many=True)
        return schema.jsonify(Command.query.filter_by(code_id=cid).all())


class CodeCommandResource(Resource):
    def get(self, cid, mid):
        schema = CodeCommandSchema()
        return schema.jsonify((Command.query.filter_by(code_id=cid, machine_id=mid).one()))

    codecommand_args = {
        'commands': fields.Nested({
            'name': fields.Str(required=True),
            'args': fields.List(fields.Str, required=True),
            'cmd': fields.Str(required=True),
            }, many=True, required=True),
        'environment': fields.Nested({
            'modules': fields.List(fields.Str, required=False),
            'variables': fields.Dict(required=False),
            }, required=True),
        }

    @apiauth.login_required
    @use_kwargs(codecommand_args)
    def post(self, cid, mid, commands, environment):
        command = (Command.query
                   .filter_by(code_id=cid, machine_id=mid)
                   .one())

        command.commands = commands
        command.environment = environment

        db.session.commit()


class ComparisonListResource(Resource):
    comparison_args = {
        'metric': fields.String(validate=lambda m: m in ['deltatest']),
        'testresult_collections': fields.DelimitedList(
            fields.UUID(validate=must_exist_in_db(TestResult2Collection)), required=True),
        }
    @use_kwargs(comparison_args)
    def post(self, metric, testresult_collections):

        trcollections = (TestResult2Collection.query
                .filter(TestResult2Collection.id.in_(testresult_collections))
                .options(joinedload("testresults"))
                .all())

        data = {
            "metric": metric,
            "testresult_collections": trcollections,
            }

        if metric == 'deltatest':
            # We are assuming here that the element per collection is unique, which may not be the case
            results = {trc.id: {tr.data['element']: tr for tr in trc.testresults} for trc in trcollections}

            # all elements for which at least one testresult is available in at least one collection
            data['elements'] = list(set(el for r in results.values() for el in r))

            # sort by atomic number
            data['elements'].sort(key=lambda el: ATOMIC_ELEMENTS[el]['num'])

            data['values'] = []

            # Generate a list of tuples of the form:
            #  (AtomicSymbol, CollectionA.id, CollectionB.id, delta)
            # with delta being None/null, if the element is missing in either CollectionA and/or CollectionB
            # This means that we may generate one delta multiple times of the same testresult is part of multiple
            # collections. Since the testresults are yielded as part of the TestResultCollections, the user
            # should be able to interprete the data himself.
            for element in data['elements']:
                for colla, collb in itertools.combinations(results.keys(), 2):
                    delta = None

                    # Ignore if an element is missing from one or both collections:
                    try:
                        trA = results[colla][element]
                        trB = results[collb][element]
                    except KeyError:
                        continue

                    try:
                        coeffs = ['V', 'B0', 'B1']
                        deltaValuesA = [trA.data['coefficients'][c] for c in coeffs]
                        deltaValuesB = [trB.data['coefficients'][c] for c in coeffs]
                        delta = calcDelta(deltaValuesA, deltaValuesB)

                        # the volume is already supposed to be per-atom
                        coeffDeltas = dict(zip(coeffs, np.abs(np.array(deltaValuesA) - np.array(deltaValuesB))))
                    except KeyError:
                        pass

                    data['values'].append(
                        {"element": element,
                         "collectionA": colla,
                         "collectionB": collb,
                         "testresultA": trA.id,
                         "testresultB": trB.id,
                         "delta": delta,
                         "coefficientDeltas": coeffDeltas,
                         })

        schema = DeltatestComparisonSchema()
        return schema.jsonify(data)


class TestListResource(Resource):
    def get(self):
        schema = TestListSchema(many=True)
        return schema.jsonify(Test.query.all())


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
api.add_resource(CalculationListActionResource, '/calculations/action')
api.add_resource(CalculationResource, '/calculations/<uuid:cid>')
api.add_resource(CalculationActionResource, '/calculations/<uuid:cid>/action')
api.add_resource(CalculationPreviewResource, '/calculations/<uuid:cid>/preview')
api.add_resource(Task2ListResource, '/calculations/<uuid:cid>/tasks', endpoint='calculationtask2listresource')
api.add_resource(Task2ListResource, '/tasks')
api.add_resource(Task2Resource, '/tasks/<uuid:tid>')
api.add_resource(Task2UploadResource, '/tasks/<uuid:tid>/uploads')
api.add_resource(ArtifactListResource, '/artifacts')
api.add_resource(ArtifactResource, '/artifacts/<uuid:aid>')
api.add_resource(ArtifactDownloadResource, '/artifacts/<uuid:aid>/download', resource_class_kwargs={'mode': 'download'})
api.add_resource(ArtifactDownloadResource, '/artifacts/<uuid:aid>/view', endpoint='artifactviewresource', resource_class_kwargs={'mode': 'view'})
api.add_resource(StructureSetListResource, '/structuresets')
api.add_resource(StructureSetResource, '/structuresets/<string:name>')
api.add_resource(StructureSetStructureListResource, '/structuresets/<string:name>/structures')
api.add_resource(StructureSetCalculationsListResource, '/structuresets/<string:name>/calculations')
api.add_resource(StructureListResource_v2, '/structures')
api.add_resource(StructureResource_v2, '/structures/<uuid:sid>')
api.add_resource(StructureDownloadResource, '/structures/<uuid:sid>/download')
api.add_resource(BasisSetListResource, '/basissets')
api.add_resource(BasisSetResource, '/basissets/<uuid:bid>')
api.add_resource(Pseudopotential2ListResource, '/pseudopotentials')
api.add_resource(Pseudopotential2Resource, '/pseudopotentials/<uuid:pid>')
api.add_resource(TestResultListResource, '/testresults')
api.add_resource(TestResultResource, '/testresults/<uuid:trid>')
api.add_resource(TestResultListActionResource, '/testresults/action')
api.add_resource(TestResultCollectionListResource, '/testresultcollections')
api.add_resource(TestResultCollectionResource, '/testresultcollections/<uuid:trcid>')
api.add_resource(ComparisonListResource, '/comparisons')
api.add_resource(ActionResource, '/actions/<uuid:aid>')
api.add_resource(CodeListResource, '/codes')
api.add_resource(CodeResource, '/codes/<uuid:cid>')
api.add_resource(CodeCommandListResource, '/codes/<uuid:cid>/commands')
api.add_resource(CodeCommandResource, '/codes/<uuid:cid>/commands/<uuid:mid>')
api.add_resource(TestListResource, '/tests')
