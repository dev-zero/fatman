#!/usr/bin/env python

import json
from os import path
import bz2
from io import BytesIO, StringIO
from uuid import UUID

from flask_restful import (Api, Resource,
                           abort, reqparse,
                           fields, marshal_with, marshal)
from flask import make_response, url_for, request
from werkzeug.datastructures import FileStorage
from werkzeug.wrappers import Response
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm import joinedload, defer, contains_eager
from sqlalchemy import func, case

from . import app, db, resultfiles, cache
from .models import (
    Structure,
    BasisSet,
    BasisSetFamily,
    Pseudopotential,
    PseudopotentialFamily,
    Method,
    Test,
    TaskStatus,
    Task,
    Result,
    TestResult,
    TestStructure,
)

from .utils import route_from
from .tools import calcDelta, eos, Json2Atoms, Atoms2Json
from .tasks import (postprocess_result_file,
                    postprocess_result_files,
                    postprocess_result,
                    )

import numpy as np


method_resource_fields = {
    'id': fields.String,
    'code': fields.Raw,
    'pseudopotential': fields.String(attribute='pseudopotential.name'),
    'basis_set': fields.String(attribute='basis_set.name'),
    'settings': fields.Raw,
    '_links': {'self': fields.Url('methodresource')},
    }

method_list_fields = {
    'id': fields.String,
    'code': fields.Raw,
    'pseudopotential': fields.String(attribute='pseudopotential.name'),
    'basis_set': fields.String(attribute='basis_set.name'),
    '_links': {'self': fields.Url('methodresource')},
    }

structure_resource_fields = {
    'id': fields.String,
    'name': fields.Raw,
    'ase_structure': fields.Raw,
    }

structure_list_fields = {
    'id': fields.String,
    'name': fields.Raw,
    }

task_resource_fields = {
    'id': fields.String,
    'ctime': fields.String,
    'mtime': fields.String,
    'machine': fields.Raw,
    'status': fields.String(attribute='status.name'),
    'method': fields.Nested(method_resource_fields),
    'structure': fields.Nested(structure_resource_fields),
    '_links': {'self': fields.Url('taskresource')},
    'priority': fields.Integer,
    }

task_list_fields = {
    'id': fields.String,
    'ctime': fields.String,
    'mtime': fields.String,
    'machine': fields.Raw,
    'status': fields.String(attribute='status.name'),
    'method': fields.Nested(method_list_fields),
    'structure': fields.Nested(structure_list_fields),
    '_links': {'self': fields.Url('taskresource')},
    }

pseudo_nested_fields = {
    'id': fields.String,
    'format': fields.Raw,
    '_links': {'self': fields.Url('pseudopotentialresource')},
    }

pseudo_list_fields = {
    'id': fields.String,
    'element': fields.Raw,
    'family': fields.String(attribute='family.name'),
    'format': fields.Raw,
    'converted_from': fields.Nested(pseudo_nested_fields, allow_null=True),
    '_links': {'self': fields.Url('pseudopotentialresource')},
    }

pseudo_resource_fields = {
    'id': fields.String,
    'element': fields.Raw,
    'family': fields.String(attribute='family.name'),
    'format': fields.Raw,
    'pseudo': fields.Raw,
    'converted_from': fields.Nested(pseudo_nested_fields, allow_null=True),
    '_links': {'self': fields.Url('pseudopotentialresource')},
    }

pseudofamily_list_fields = {
    'name': fields.Raw,
    }

result_resource_fields = {
    'id': fields.String,
    'energy': fields.Float,
    'task': fields.Nested(task_resource_fields),
    '_links': {'self': fields.Url('resultresource'),
               'file': fields.Url('resultfileresource')},
    'filename': fields.String,
    'data': fields.Raw,
    }

test_resource_fields = {
    'id': fields.Raw,
    'name': fields.String,
    }

teststructure_resource_fields = {
    'test': fields.Nested(test_resource_fields),
    'structure': fields.Nested(structure_resource_fields),
    }

testresult_resource_fields = {
    'id': fields.String,
    'test': fields.Nested(test_resource_fields),
    'method': fields.Nested(method_resource_fields),
    'ctime': fields.String,
    'data': fields.Raw(attribute='result_data'),
    }


class ParameterError(Exception):
    pass


class TaskResource(Resource):
    @marshal_with(task_resource_fields)
    def get(self, id):
        return (Task.query
                .options(joinedload("method"))
                .options(joinedload("structure"))
                .get_or_404(id))

    @marshal_with(task_resource_fields)
    def patch(self, id):
        parser = reqparse.RequestParser()
        parser.add_argument('status', type=str, required=False)
        parser.add_argument('machine', type=str, required=False)
        parser.add_argument('priority', type=int, required=False)
        args = parser.parse_args()

        # update the status and/or priority (mtime is updated by SQLA)
        task = Task.query.get_or_404(id)

        if args['status'] is not None:
            task.status = TaskStatus.query.filter_by(name=args['status']).one()
        if args['machine'] is not None:
            task.machine = args['machine']
        if args['priority'] is not None:
            task.priority = args['priority']

        db.session.commit()

        return task


class TaskList(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('status', type=str)
        parser.add_argument('limit', type=int)
        parser.add_argument('timeorder', type=bool, default=False)
        parser.add_argument('machine', type=str)
        parser.add_argument('structure', type=str)
        parser.add_argument('pseudofamily', type=str)
        args = parser.parse_args()

        q = (Task.query
             .options(joinedload("method"))
             .options(joinedload("structure"))
             .join(TaskStatus)
             .options(contains_eager("status")))

        if args['timeorder']:
            q = q.order_by(Task.mtime.desc())

        if args['status'] is not None:
            q = q.filter(TaskStatus.name == args['status'])

        if args['machine'] is not None:
            q = q.filter_by(machine=args['machine'])

        if args['structure'] is not None:
            # explicit join is required since SQLA separates
            # pre-loading of ORM objects and query-making
            q = q.join(Structure).filter(
                Structure.name.contains(args['structure']))

        if args['pseudofamily'] is not None:
            q = q.join(Method).filter(
                Method.pseudopotential.has(name=args['pseudofamily']))

        if args['limit'] is not None:
            q = q.limit(args['limit'])

        return [marshal(t, task_list_fields) for t in q]

    def post(self):
        ret = []
        parser = reqparse.RequestParser()
        parser.add_argument('structure', type=str, required=False)
        parser.add_argument('method', type=UUID, required=True)
        parser.add_argument('status', type=str, default="new")
        parser.add_argument('test', type=str, required=True)
        parser.add_argument('priority', type=int, default=0)
        args = parser.parse_args()

        m = Method.query.get_or_404(args['method'])
        s = TaskStatus.query.filter_by(name=args['status']).one()
        t = Test.query.filter_by(name=args['test']).one()

        if args['structure'] is not None:
            # not using .all() to get exception if nothing found
            structs = [Structure.query.filter_by(name=args['structure']).one()]
        else:
            structs = t.structures

        for struct in structs:
            task = (Task.query
                    .filter_by(structure=struct, method=m, test=t)
                    .first())

            if not task:
                task = Task(structure=struct, method=m, test=t, status=s,
                            machine='-', priority=args['priority'])
                db.session.add(task)
                db.session.commit()

            ret.append(str(task.id))

        return ret


class ResultResource(Resource):
    @marshal_with(result_resource_fields)
    def get(self, id):
        return Result.query.get_or_404(id)

    @marshal_with(result_resource_fields)
    def patch(self, id):
        parser = reqparse.RequestParser()
        parser.add_argument('data', type=str, required=False)
        parser.add_argument('energy', type=float, required=False)
        parser.add_argument('task', type=UUID, required=False)
        args = parser.parse_args()

        res = Result.query.get_or_404(id)

        if args['data'] is not None:
            res.data = json.loads(args['data'])
        if args['energy'] is not None:
            res.energy = args['energy']
        if args['task'] is not None:
            res.task = Task.get(Task.id == args['task'])

        db.session.commit()

        return res


class ResultFileResource(Resource):
    def post(self, id):
        result = Result.query.get_or_404(id)

        if result.filename is not None:
            abort(400, message="Data is already uploaded for this result")

        parser = reqparse.RequestParser()
        parser.add_argument('file', required=True,
                            type=FileStorage, location='files')
        args = parser.parse_args()

        filename = resultfiles.save(args['file'], folder=str(result.id))

        result.filename = filename

        db.session.commit()

        postprocess_result_file.delay(id)

        # use a raw werkzeug Response object to return 201
        # status code without a body
        return Response(status=201)

    def get(self, id):
        result = Result.query.get_or_404(id)

        if result.filename is None:
            abort(400, message="No data available")

        # return the file
        with bz2.BZ2File(resultfiles.path(result.filename)) as infile:
            response = make_response(infile.read())
            response.headers["Content-Disposition"] = \
                "attachment; filename={:}".format(
                        path.splitext(result.filename)[0])
            response.headers["Content-Type"] = "text/plain"
            return response


class ResultActionResource(Resource):
    def get(self, rid, action, tid):
        Result.query.get_or_404(rid)

        if action == 'doPostprocessing':
            async_result = postprocess_result_file.AsyncResult(str(tid))

            if not async_result.ready():
                return Response(status=202, headers={'Location': request.path})

            resp = {'status': async_result.status,
                    'parent': api.url_for(ResultResource, id=rid)}

            if async_result.successful():
                resp['result'] = {'updated': async_result.result}

            return resp

        raise ParameterError("invalid action specified: {}".format(action))


class ResultActionList(Resource):
    def post(self, rid):
        Result.query.get_or_404(rid)

        parser = reqparse.RequestParser()
        parser.add_argument('doPostprocessing', type=dict, required=False)
        args = parser.parse_args()

        if args['doPostprocessing'] is not None:
            if ('update' in args['doPostprocessing'].keys() and not
               type(args['doPostprocessing']['update']) == bool):
                raise ParameterError(
                        "the 'update' parameter must be a boolean")

            update = args['doPostprocessing'].get('update', False)
            async_result = postprocess_result_file.delay(rid, update)

            return Response(status=202, headers={'Location': api.url_for(
                        ResultActionResource,
                        rid=rid, action='doPostprocessing',
                        tid=async_result.id)})

        raise ParameterError("invalid action specified")


class ResultsActionResource(Resource):
    def get(self, action, tid):

        if action == 'doPostprocessing':
            group_result = postprocess_result_files.AsyncResult(str(tid))

            if not group_result.ready():
                return Response(status=202, headers={'Location': request.path})

            rids, results = group_result.result
            all_done = all(r.ready() for r in results)

            resp = [{'status': res.status,
                     '_links': {'parent': api.url_for(ResultResource, id=rid)},
                     'result': {'updated': res.result} if res.successful()
                     else {}}
                    for rid, res in zip(rids, results)]

            if all_done:
                return resp
            else:
                return resp, 202, {'Location': request.path}

        raise ParameterError("invalid action specified: {}".format(action))


class ResultsActionList(Resource):
    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('doPostprocessing', type=dict, required=False)
        args = parser.parse_args()

        if args['doPostprocessing'] is not None:
            if ('update' in args['doPostprocessing'].keys() and not
               type(args['doPostprocessing']['update']) == bool):
                raise ParameterError(
                        "the 'update' parameter must be a boolean")

            update = args['doPostprocessing'].get('update', False)

            async_result = postprocess_result_files.delay(update)

            return Response(status=202, headers={'Location': api.url_for(
                        ResultsActionResource, action='doPostprocessing',
                        tid=async_result.id)})

        raise ParameterError("invalid action specified")


class ResultList(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('test', type=str)
        parser.add_argument('method', type=UUID)
        parser.add_argument('structure', type=str)
        parser.add_argument('calculated_on', type=str)
        parser.add_argument('code', type=str)
        args = parser.parse_args()

        # ensure we load structure and method in one go
        q = (db.session.query(Result).join(Task)
             .options(joinedload("task").joinedload("structure"))
             .options(joinedload("task").joinedload("method"))
             )

        if args['test'] is not None:
            q = q.join(Test).filter_by(name=args['test'])

        if args['method'] is not None:
            q = q.join(Method).filter_by(id=args['method'])

        if args['structure']:
            q = q.join(Structure).filter_by(name=args['structure'])

        if args['calculated_on']:
            q = q.filter_by(calculated_on=args['method'])

        if args['code']:
            q = q.join(Method).filter_by(code=args['code'])

        return [marshal(x, result_resource_fields) for x in q]

    @marshal_with(result_resource_fields)
    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('energy', type=float, required=True)
        parser.add_argument('task', type=str)
        parser.add_argument('task_id', type=UUID)
        parser.add_argument('data', type=str)
        args = parser.parse_args()

        task_id = args['task_id']

        if args['task']:
            url, data = route_from(args['task'], 'GET')
            if url != 'taskresource':
                raise ParameterError("Invalid URL specified for task")
            task_id = data['id']

        if not task_id:
            raise ParameterError("Either task or task_id must be specified")

        if args['data'] is not None:
            extradata = json.loads(args['data'])
        else:
            extradata = None

        result = Result(task_id=task_id,
                        energy=args['energy'],
                        data=extradata)
        db.session.commit()

        postprocess_result.delay(result.id)

        return (result, 201,
                {'Location': api.url_for(ResultResource, id=result.id)})


class TestResource(Resource):
    def get(self, testname):
        teststructures = (db.session.query(TestStructure, Test, Structure)
                          .join(Structure).join(Test)
                          .filter_by(name=testname))
        return [marshal({'test': ts[2], 'structure': ts[3]},
                        teststructure_resource_fields)
                for ts in teststructures]


class TestList(Resource):
    def get(self):
        return [(t.id, str(t)) for t in Test.query]


class MethodList(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('test', type=str)
        args = parser.parse_args()

        q = Method.query

        if args['test']:
            q = (q.join(TestResult).join(Test)
                 .filter_by(name=args['test']))

        return [marshal(m, method_list_fields) for m in q]

    @marshal_with(method_resource_fields)
    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('code', type=str, required=True)
        # TODO: the column and the attribute are called pseudopotential,
        #       but are in fact family-references
        parser.add_argument('pseudopotential', type=str, required=True)
        parser.add_argument('basis_set', type=str, required=True)
        parser.add_argument('settings', type=str, required=True)
        args = parser.parse_args()

        bsfamily = BasisSetFamily.query.filter_by(name=args['basis_set']).one()
        pfamily = (PseudopotentialFamily.query
                   .filter_by(name=args['pseudopotential']).one())

        settings = json.loads(args['settings'])

        try:
            return (Method.query
                    .filter_by(code=args['code'])
                    .filter_by(basis_set=bsfamily)
                    .filter_by(pseudopotential=pfamily)
                    .filter_by(settings=settings)
                    .one())
        except NoResultFound:
            method = Method(code=args['code'],
                            basis_set=bsfamily,
                            pseudopotential=pfamily,
                            settings=settings)
            db.session.commit()

        return method


class MethodResource(Resource):
    """Return all the details for a method, or set them"""

    @marshal_with(method_resource_fields)
    def get(self, id):
        return Method.query.get_or_404(id)


class Basissets(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('family', type=str, required=True)
        parser.add_argument('element', type=str, required=True,
                            action="append")
        args = parser.parse_args()

        q = (BasisSet.query
             .filter_by(BasisSet.family.has(name=args['family']))
             .filter_by(BasisSet.element.in_(args['element'])))

        return {b.element: b.basis for b in q}


class PseudopotentialResource(Resource):
    @marshal_with(pseudo_resource_fields)
    def get(self, id):
        return (Pseudopotential.query
                .options((joinedload("converted_from")
                          .defer('pseudo')))
                .get_or_404(id))


class PseudopotentialFamilyList(Resource):
    """Get the list of pseudopotential families"""
    def get(self):
        q = (PseudopotentialFamily.query
             .order_by(PseudopotentialFamily.name))

        return [marshal(f, pseudofamily_list_fields) for f in q]


class PseudopotentialList(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('family', type=str)
        parser.add_argument('format', type=str)
        parser.add_argument('element', type=str, action="append")
        args = parser.parse_args()

        # never load the pseudo column for converted_from entries
        # (would take much more time and memory)
        q = (Pseudopotential.query
             .options((joinedload("converted_from")
                       .defer('pseudo'))))

        if args['family']:
            q = (q.join(PseudopotentialFamily)
                  .filter(PseudopotentialFamily.name == args['family']))

        if args['format']:
            q = q.filter_by(format=args['format'])

        if args['element']:
            q = q.filter(Pseudopotential.element.in_(args['element']))

        if any(args.values()):
            return [marshal(p, pseudo_resource_fields) for p in q]
        else:
            # defer loading of the pseudo column for the primary
            # pseudo as well (reduces loading time by 80%)
            q = q.options(defer('pseudo'))
            return [marshal(p, pseudo_list_fields) for p in q]

    @marshal_with(pseudo_resource_fields)
    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('family', type=str, required=True)
        parser.add_argument('element', type=str, required=True)
        parser.add_argument('pseudo', type=str, required=True)
        parser.add_argument('format', type=str, required=True)
        parser.add_argument('converted_from', type=UUID, required=False)
        parser.add_argument('overwrite', type=bool, default=False)

        args = parser.parse_args()

        family = (PseudopotentialFamily.query
                  .filter_by(name=args['family'])
                  .one())

        data = {'family': family,
                'element': args['element'],
                'format': args['format']}

        if args['converted_from']:
            data['converted_from'] = args['converted_from']

        try:
            pseudo = (Pseudopotential.query
                      .filter(**data)
                      .one())

            if args['overwrite']:
                pseudo.pseudo = args['pseudo']
                db.session.commit()
                return (pseudo, 200,
                        {'Location': api.url_for(PseudopotentialResource,
                                                 id=pseudo.id)})

        except NoResultFound:
            pseudo = Pseudopotential(pseudo=args['pseudo'], **args)
            db.session.commit()
            return (pseudo, 201,
                    {'Location': api.url_for(PseudopotentialResource,
                                             id=pseudo.id)})

        abort(400, message="Pseudo is already uploaded for this result")


class MachineStatus(Resource):
    """Return a dictionary with the number of running and total tasks per machine:
              MACHINE   RUN   TOTAL
       e.g.   machine1   4    124
              machine2   5    24
              machine3   0    242
     """
    def get(self):
        return (db.session.query(Task.machine,
                                 func.count(Task.id),
                                 func.count(
                                     case([(Task.status.has(name="running"),
                                            Task.id)]))
                                 )
                .group_by(Task.machine)
                .all())


class TestResultList(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('method', type=UUID, required=False)
        parser.add_argument('test', type=str, required=False)
        # method ID of the reference for the deltatest:
        parser.add_argument('deltaref', type=UUID, required=False)
        parser.add_argument('limit', type=int, required=False)
        args = parser.parse_args()

        query = (TestResult.query
                 .join(Test)
                 .join(Method)
                 .options(contains_eager(TestResult.method))
                 .options(contains_eager(TestResult.test))
                 .order_by(TestResult.ctime.desc()))

        if args['test']:
            query = query.filter(Test.name.contains(args["test"]))

        if args['method']:
            query = query.filter(Method.id == args["method"])

        if args['limit']:
            query = query.limit(args['limit'])

        results = [marshal(tr, testresult_resource_fields)
                   for tr in query]

        if args['deltaref']:
            # get the test result entries for the reference method:
            ref_res_query = (TestResult.query
                             .options(joinedload("test"))
                             .filter_by(method_id=args['deltaref']))

            # create a dictionary for faster lookup with the test name as key:
            refresults = {r.test.name: r for r in ref_res_query}

            for result in results:
                # note: the results here are now dicts, not ORM objects

                # ignore results with missing V,B0,B1 coefficients
                if not result['data']['_status'] == 'fitted':
                    continue

                try:
                    refdata = refresults[result['test']['name']].result_data

                    result['data']['deltavalue'] = calcDelta(
                        (refdata['V'], refdata['B0'], refdata['B1']),
                        (result['data']['V'],
                         result['data']['B0'],
                         result['data']['B1']))

                except KeyError:
                    # do not calculate delta if the ref method doesn't have it
                    continue

        return results


class Plot(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('method', type=UUID, action="append")
        parser.add_argument('test', type=str, required=True)
        args = parser.parse_args()

        test1 = Test.query.filter_by(name=args['test']).one()

        from matplotlib.backends.backend_agg import FigureCanvasAgg \
            as FigureCanvas
        import matplotlib.pyplot as plt
        import itertools
        colors = itertools.cycle([
            "#1D0662", "#8F1E00", "#00662C", "#8F7700",
            "#3E1BA7", "#F33D0D", "#0AAD51", "#F3CD0D",
            "#8E7BC5", "#FFAA93", "#75CA9A", "#FFED93"])

        fig = plt.figure(figsize=(12, 8), dpi=200, facecolor="#FFFFFF")
        fig.subplots_adjust(left=0.07, right=0.98, top=0.99, bottom=0.04)
        ax = plt.subplot(111)
        stride = int(min(35./len(args['method']), 7))
        label_xpos = 5
        warning_yshift = 0
        i = 0

        for mid in args['method']:
            i += 1

            method = Method.query.get_or_404(mid)

            # here the curve, based on the fitted data
            try:
                r = TestResult.query.filter_by(method=method, test=test1).one()
            except NoResultFound:
                ax.text(0.5, 0.5+warning_yshift,
                        "No data for method {:}".format(method),
                        transform=ax.transAxes,
                        color="#FF0000")
                warning_yshift += 0.05
                continue

            if not (r.result_data['_status'] == 'fitted' or
                    r.result_data['_status'] == 'reference'):
                ax.text(0.5, 0.5+warning_yshift,
                        "Fitting problem for method {:}".format(method),
                        transform=ax.transAxes,
                        color="#FF0000")
                warning_yshift += 0.05
                continue

            mycolor = next(colors)
            V0 = r.result_data['V']
            B0 = r.result_data['B0']
            B1 = r.result_data['B1']
            xfit, yfit = eos(V0, B0, B1)
            ax.plot(xfit, yfit, color=mycolor)

            # here the points
            q = (Result.query
                 .join(Task)
                 .join(Structure)
                 .filter(Task.method == method)
                 .filter(Task.test == test1))

            E = []
            V = []
            for x in q:
                atoms = Json2Atoms(x.task.structure.ase_structure)
                natom = len(atoms.get_masses())
                V.append(atoms.get_volume()/natom)
                E.append(x.energy/natom)
            x2 = np.array(V)
            y2 = np.array(E)

            if 'E0' in r.result_data.keys():
                E0 = r.result_data['E0']
                y2 -= E0

            ax.plot(x2, y2, 's', color=mycolor)

            ax.annotate("Method {:}".format(method),
                        xy=(xfit[label_xpos], yfit[label_xpos]),
                        xytext=(-5, -30 if i != 3 else 20),
                        textcoords="offset points",
                        fontsize=10,
                        arrowprops={'facecolor': 'black',
                                    'shrink': 0.05,
                                    'headwidth': 3,
                                    'width': 1},
                        horizontalalignment='right',
                        color=mycolor)
            label_xpos += stride

        canvas = FigureCanvas(fig)
        png_output = BytesIO()
        canvas.print_png(png_output)
        response = make_response(png_output.getvalue())
        response.headers['Content-Type'] = 'image/png'
        return response


class StructureResource(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('id', type=UUID, required=False)
        parser.add_argument('test', type=str, required=False)
        parser.add_argument('repeat', type=int, required=False, default=1)
        parser.add_argument('viewer', type=bool, required=False, default=False)
        parser.add_argument('size', type=int, required=False, default=300)
        parser.add_argument('format', type=str, required=False,
                            default="xyz", choices=("xyz", "json"))
        args = parser.parse_args()

        if args['id'] and args['test']:
            raise ParameterError(
                    "exactly one of id or test parameter is required")

        if args['id'] is not None:
            s = Structure.query.get_or_404(args['id'])
        elif args['test']:
            s = (Structure.query
                 .filter(Structure.tests.any(name=args['test']))
                 .first())
        else:
            raise ParameterError("exactly one of id or test parameter is required")

        atoms = Json2Atoms(s.ase_structure)

        if args['viewer']:
            from ase.io import cif
            mycif = StringIO()
            cif.write_cif(mycif, atoms)
            viewer_html = """
                <html>
                    <head>
                        <link rel="stylesheet" href="https://hub.chemdoodle.com/cwc/latest/ChemDoodleWeb.css" type="text/css">
                        <script type="text/javascript" src="https://hub.chemdoodle.com/cwc/latest/ChemDoodleWeb.js"></script>
                    </head>
                    <body>
                        <script>
                            var t = new ChemDoodle.TransformCanvas3D('transformBallAndStick', {0:d}, {0:d});
                            t.specs.projectionPerspective_3D = false;
                            t.specs.atoms_font_size_2D = 12;
                            t.specs.atoms_useVDWDiameters_3D = true;
                            t.specs.atoms_vdwMultiplier_3D = 8;
                            t.specs.set3DRepresentation('Ball and Stick');
                            t.specs.backgroundColor = 'white';
                            t.specs.atoms_displayLabels_3D = true;
                            cif=`{1:}`;
                            var molecule = ChemDoodle.readCIF(cif, {2:d},{2:d},{2:d});
                            t.loadContent([molecule.molecule],[molecule.unitCell]);
                        </script><br/>
                        <a href="https://web.chemdoodle.com/"><small>Powered by ChemDoodle Web Components (GPL)</small></a>
                    </body></html>""".format(args['size'], mycif.getvalue(), args['repeat'])
            return make_response(viewer_html)

        if args['format'] == 'json':
            atoms_rep = atoms.repeat(args['repeat'])
            return Atoms2Json(atoms_rep)

        atoms_rep = atoms.repeat(args['repeat'])
        xyz_output = ""
        xyz_output += "{:d}\n".format(len(atoms_rep.numbers))
        xyz_output += "CELL: {:}\n".format(str(atoms.get_cell()).replace('\n', ';'))
        xyz_output += "\n".join(["{:}   {:10.6f} {:10.6f} {:10.6f}" \
                .format(x[0], x[1][0], x[1][1], x[1][2]) for x in
                                 zip(atoms_rep.get_chemical_symbols(),
                                     atoms_rep.get_positions())
                                ])
        response = make_response(xyz_output)
        response.headers["Content-Type"] = "text/plain"

        return response


class Comparison(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('method1', type=UUID, required=True)
        parser.add_argument('method2', type=UUID, required=False)
        parser.add_argument('test', type=str, action="append")
        args = parser.parse_args()

        # COMPARE all tests for 1 reference method
        # if neither method2 nor test are not specified
        mode = "1method"

        method1 = Method.query.get_or_404(args["method1"])

        if args["method2"] is not None:
            # COMPARE 2 METHODS
            mode = "2methods"
            method2 = Method.query.get_or_404(args["method2"])
        elif args["method2"] is None and args["test"] is not None:
            # COMPARE all tests for 1 reference method
            mode = "methodbytest"
            if len(args["test"]) > 1:
                raise ParameterError(("Specify either 2 methods "
                                      "and optionally a list of tests, "
                                      "or 1 method and 1 test"))

        ret = {"test": {}, "methods": [], "method": {}, "summary": {}}

        if mode == "methodbytest":
            ret["methods"] = [str(method1.id)]

            testname = args["test"][0]

            for method2 in Method.query.all():
                result_data = self._getResultData(method1.id, method2.id,
                                                  testname)
                if result_data:
                    ret["method"][str(method2.id)] = [str(method2)] + result_data

        elif mode == "2methods":
            return self._get_2methods(method1.id, method2.id, args['test'])
        elif mode == "1method":
            return self._get_1method(method1.id)

        return ret

    @staticmethod
    @cache.memoize(timeout=3600)
    def _get_2methods(method1_id, method2_id, testlist=None):
        ret = {"test": {}, "methods": [], "method": {}, "summary": {}}

        if not testlist:
            testlist = [t.name for t in Test.query]

        all_delta = []
        for testname in testlist:
            result_data = Comparison._getResultData(method1_id, method2_id,
                                                    testname)
            if result_data:
                ret["test"][testname] = result_data
                all_delta.append(result_data[-1])

        all_delta = np.array(all_delta)
        ret["summary"] = {"avg": np.average(all_delta),
                          "stdev": np.std(all_delta),
                          "N": len(all_delta)}
        ret["methods"] = [str(method1_id), str(method2_id)]

        return ret

    @staticmethod
    @cache.memoize(timeout=3600)
    def _get_1method(method1_id):
        """We need the ID here and not the object to get a proper cache key"""

        ret = {"test": {}, "methods": [], "method": {}, "summary": {}}

        testlist = [t.name for t in Test.query.all()]

        # loop over method2
        q2 = Method.query.order_by(Method.id).all()
        for method2 in q2:
            ret["methods"].append(str(method2.id))
            all_delta = []
            for testname in testlist:
                result_data = Comparison._getResultData(method1_id, method2.id,
                                                        testname)
                if result_data:
                    all_delta.append(result_data[-1])

            if len(all_delta) > 0:
                all_delta = np.array(all_delta)
                ret['method'][method2.id] = [str(method2),
                                             np.average(all_delta),
                                             np.std(all_delta),
                                             len(all_delta)]
            else:
                ret['method'][method2.id] = [str(method2), -1, -1, 0]

        return ret

    @staticmethod
    def _getResultData(method1_id, method2_id, testname):
        """Fetch test results from the database and calculate some metric.

        This returns a list of values.
        The meaning of those values depends on the respective test:

        For deltatest-tests: [V_m1, B0_m1, B1_m1, V_m2, B0_m2, B1_m2, Δ].
        For GMTKN: [Δ]

        In the following cases False is returned:
        * the specified test is missing for at least one of the methods
        * the specified test is unknown
        * the calculation the relevant metric fails
        """

        try:
            r1 = (TestResult.query
                  .filter_by(method_id=method1_id)
                  .join(Test)
                  .filter(Test.name == testname)
                  .one())
            r2 = (TestResult.query
                  .filter_by(method_id=method2_id)
                  .join(Test)
                  .filter(Test.name == testname)
                  .one())
        except Exception as e:
            app.logger.warning(('Invalid method ids specified (%d, %d) or '
                                'test %s does not exist for both methods: '
                                '%s'), method1_id, method2_id, testname, e)
            return False

        try:
            if "deltatest" in testname:
                data_f = [r1.result_data["V"],
                          r1.result_data["B0"],
                          r1.result_data["B1"]]
                data_w = [r2.result_data["V"],
                          r2.result_data["B0"],
                          r2.result_data["B1"]]

                delta = calcDelta(data_f, data_w)

                return data_f + data_w + [delta]

            elif "GMTKN" in testname:
                n = 0
                mad = 0.

                for e1, e2 in zip(r1.result_data["energies"],
                                  r2.result_data["energies"]):
                    mad += abs(e1-e2)
                    n += 1

                delta = mad/n

                return [delta]

        except Exception as e:
            app.logger.error(('Calculating the metric between methods '
                              '%d and %d failed for %s: %s'),
                             method1_id, method2_id, testname, e)
            return False

        app.logger.warning('Invalid testname specified to get results: %s',
                           testname)
        return False


class StatsResource(Resource):
    """Return a list of summary stats, including links to go to the details:
       e.g. for tasks
         [ { "status": "running", "count": 7, "_links": { "self": "/tasks?status=running" }},
           { "status": "error", "count": 2, "_links": { "self": "/tasks?status=error" }},
           { "status": "total", "count": 9, "_links": { "self": "/tasks" }} ]
    """
    def get(self, what):
        if what == "tasks":

            q = (db.session.query(TaskStatus.name,
                                  func.count(Task.id).label("count"))
                 .outerjoin(TaskStatus.tasks)
                 .group_by(TaskStatus.name)
                 .order_by(TaskStatus.name))

            stats = [{"status": e.name,
                      "count": e.count,
                      "_links": {"self": url_for("tasklist", status=e.name)},
                      } for e in q]

            return stats

        abort(404)

# Catch common exceptions in the REST dispatcher
errors = {
        'TaskStatusDoesNotExist': {
            'message': "Invalid status specified.",
            'status': 404,
            },
        'TaskDoesNotExist': {
            'message': "Task with the specified ID does not exist.",
            'status': 404,
            },
        'ResultDoesNotExist': {
            'message': "Result with the specified ID does not exist.",
            'status': 404,
            },
        'ParameterError': {
            'message': "Specify either 2 methods and optionally a list of tests, or 1 method and 1 test.",
            'status': 404,
            },
        }


api = Api(app, errors=errors)

api.add_resource(TaskResource, '/tasks/<uuid:id>')
api.add_resource(TaskList, '/tasks')
api.add_resource(ResultResource, '/results/<uuid:id>')
api.add_resource(ResultsActionResource,
                 '/results/action/<string:action>:<uuid:tid>')
api.add_resource(ResultsActionList, '/results/action')
api.add_resource(ResultFileResource, '/results/<uuid:id>/file')
api.add_resource(ResultActionList,
                 '/results/<uuid:rid>/action')
api.add_resource(ResultActionResource,
                 '/results/<uuid:rid>/action/<string:action>:<uuid:tid>')
api.add_resource(ResultList, '/results')
api.add_resource(Basissets, '/basis')
api.add_resource(PseudopotentialResource, '/pseudos/<uuid:id>')
api.add_resource(PseudopotentialList, '/pseudos')
api.add_resource(PseudopotentialFamilyList, '/pseudofamilies')
api.add_resource(MachineStatus, '/machinestatus')
api.add_resource(TestResultList, '/testresults')
api.add_resource(Comparison, '/compare')
api.add_resource(MethodResource, '/methods/<uuid:id>')
api.add_resource(MethodList, '/methods')
api.add_resource(TestList, '/tests')
api.add_resource(TestResource, '/tests/<string:testname>')
api.add_resource(Plot, '/plot')
api.add_resource(StructureResource, '/structure')
api.add_resource(StatsResource, '/stats/<string:what>')

#  vim: set ts=4 sw=4 tw=0 :
