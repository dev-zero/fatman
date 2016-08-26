#!/usr/bin/env python

from datetime import datetime as dt
import json
from os import path
import bz2
from io import BytesIO, StringIO

from flask_restful import Api, Resource, abort, reqparse, fields, marshal_with, marshal
from flask import make_response, url_for, request
from playhouse.shortcuts import model_to_dict
from werkzeug.datastructures import FileStorage
from werkzeug.wrappers import Response

from . import app, resultfiles, cache
from .models import *
from .utils import route_from
from .tools import calcDelta, eos, Json2Atoms, Atoms2Json
from .tasks import postprocess_result_file, postprocess_result_files

import numpy as np

method_resource_fields = {
    'id': fields.Raw,
    'code': fields.Raw,
    'pseudopotential': fields.String(attribute='pseudopotential.name'),
    'basis_set': fields.String(attribute='basis_set.name'),
    'settings': fields.Raw,
    '_links': {'self': fields.Url('methodresource')},
    }

method_list_fields = {
    'id': fields.Raw,
    'code': fields.Raw,
    'pseudopotential': fields.String(attribute='pseudopotential.name'),
    'basis_set': fields.String(attribute='basis_set.name'),
    '_links': {'self': fields.Url('methodresource')},
    }

structure_resource_fields = {
    'id': fields.Raw,
    'name': fields.Raw,
    'ase_structure' : fields.Raw, 
    }

structure_list_fields = {
    'id': fields.Raw,
    'name': fields.Raw,
    }

task_resource_fields = {
    'id': fields.Raw,
    'ctime': fields.String,
    'mtime': fields.String,
    'machine': fields.Raw,
    'status': fields.String(attribute='status.name'),
    'method': fields.Nested(method_resource_fields),
    'structure': fields.Nested(structure_resource_fields),
    '_links': {'self': fields.Url('taskresource')},
    'priority' : fields.Integer,
    }

task_list_fields = {
    'id': fields.Raw,
    'ctime': fields.String,
    'mtime': fields.String,
    'machine': fields.Raw,
    'status': fields.String(attribute='status.name'),
    'method': fields.Nested(method_list_fields),
    'structure': fields.Nested(structure_list_fields),
    '_links': {'self': fields.Url('taskresource')},
    }

test_resource_fields = {
    'id': fields.Raw,
    'test': fields.Raw,
    'structure': fields.Nested(structure_resource_fields),
    }

pseudo_nested_fields = {
    'id': fields.Raw,
    'format': fields.Raw,
    '_links': {'self': fields.Url('pseudopotentialresource')},
    }

pseudo_list_fields = {
    'id': fields.Raw,
    'element': fields.Raw,
    'family': fields.String(attribute='family.name'),
    'format': fields.Raw,
    'converted_from': fields.Nested(pseudo_nested_fields, default={}),
    '_links': {'self': fields.Url('pseudopotentialresource')},
    }

pseudo_resource_fields = {
    'id': fields.Raw,
    'element': fields.Raw,
    'family': fields.String(attribute='family.name'),
    'format': fields.Raw,
    'pseudo': fields.Raw,
    'converted_from': fields.Nested(pseudo_nested_fields, default={}),
    '_links': {'self': fields.Url('pseudopotentialresource')},
    }

pseudofamily_list_fields = {
    'name': fields.Raw,
    }

result_resource_fields = {
    'id': fields.Raw,
    'energy': fields.Float,
    'task': fields.Nested(task_resource_fields),
    '_links': {'self': fields.Url('resultresource'),
               'file': fields.Url('resultfileresource')},
    'filename': fields.String,
    'data': fields.Raw,
    }


class ParameterError(Exception):
    pass


class TaskResource(Resource):
    @marshal_with(task_resource_fields)
    def get(self, id):
        task = Task.select(Task, TaskStatus, Method, Structure, Test, PseudopotentialFamily, BasissetFamily) \
            .where(Task.id == id) \
            .join(TaskStatus).switch(Task) \
            .join(Method) \
               .join(PseudopotentialFamily).switch(Method) \
               .join(BasissetFamily).switch(Method) \
               .switch(Task) \
            .join(Structure).switch(Task) \
            .join(Test).switch(Task) \
            .get()

        return model_to_dict(task)

    @marshal_with(task_resource_fields)
    def patch(self, id):
        parser = reqparse.RequestParser()
        parser.add_argument('status', type=str, required=True)
        parser.add_argument('machine', type=str, required=False)
        parser.add_argument('priority', type=int, required=False)
        args = parser.parse_args()

        # update the status and/or priority and reset the modification time
        task = Task.get(Task.id == id)
        task.status = TaskStatus.get(TaskStatus.name == args['status']).id
        task.mtime = dt.now()
        if 'machine' in args.keys() and args['machine'] is not None:
            task.machine = args['machine'] 
        if 'priority' in args.keys() and args['priority'] is not None:
            task.priority = args['priority'] 
        task.save()

        return model_to_dict(task)

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

        q = Task.select(Task, TaskStatus, Method, Structure, Test, PseudopotentialFamily, BasissetFamily) \
            .join(TaskStatus).switch(Task) \
            .join(Method) \
               .join(PseudopotentialFamily).switch(Method) \
               .join(BasissetFamily).switch(Method) \
               .switch(Task) \
            .join(Structure).switch(Task) \
            .join(Test).switch(Task) \
            .order_by(Task.priority.desc())

        if args['timeorder']:
            q = q.order_by(Task.mtime.desc())

        if args['status'] is not None:
            status = TaskStatus.get(TaskStatus.name == args['status'])
            q = q.where(Task.status == status)

        if args['machine'] is not None:
            q = q.where(Task.machine == args['machine'])

        if args['structure'] is not None:
            q = q.where(Structure.name.contains(args['structure']))

        if args['pseudofamily'] is not None:
            q = q.where(PseudopotentialFamily.name == args['pseudofamily'])

        if args['limit'] is not None:
            q = q.limit(args['limit'])

        return [marshal(model_to_dict(t), task_list_fields) for t in q]

    def post(self):
        ret = []
        parser = reqparse.RequestParser()
        parser.add_argument('structure', type=str, required=False)
        parser.add_argument('method', type=int, required=True)
        parser.add_argument('status', type=str, default="new")
        parser.add_argument('test', type=str, required=True)
        parser.add_argument('priority', type=int, default=0)
        args = parser.parse_args()

        m = Method.get(Method.id == args['method'])
        s = TaskStatus.get(TaskStatus.name == args['status'])
        t = Test.get(Test.name == args['test'])

        if 'structure' in args.keys() and args['structure'] is not None:
            struct = Structure.get(Structure.name == args['structure'])
            idlist = [struct.id]
        else:
            q = TestStructure.select().where(TestStructure.test == t)
            q.execute()
            idlist = [x.structure.id for x in q]

        for id in idlist:
            struct = Structure.get(Structure.id == id)
            ta, created = Task.get_or_create(structure = struct,
                                             method = m, 
                                             test = t,
                                             defaults = dict(ctime = dt.now(),
                                                             mtime = dt.now(),
                                                             status = s,
                                                             machine = '-',
                                                             priority = args['priority']))
            ret.append(ta.id)

        return ret


class ResultResource(Resource):
    @marshal_with(result_resource_fields)
    def get(self, id):
        return model_to_dict(Result.get(Result.id == id))

    @marshal_with(result_resource_fields)
    def patch(self, id):
        parser = reqparse.RequestParser()
        parser.add_argument('data', type=str, required=False)
        parser.add_argument('energy', type=float, required=False)
        parser.add_argument('task', type=int, required=False)
        args = parser.parse_args()

        res = Result.get(Result.id==id)

        if args['data'] is not None:
            res.data = json.loads(args['data'])
        if args['energy'] is not None:
            res.energy = args['energy']
        if args['task'] is not None:
            res.task = Task.get(Task.id == args['task'])

        res.save()

        return model_to_dict(res)


class ResultFileResource(Resource):
    def post(self, id):
        result = Result.get(Result.id == id)

        if result.filename is not None:
            abort(400, message="Data is already uploaded for this result")

        parser = reqparse.RequestParser()
        parser.add_argument('file', required=True,
                            type=FileStorage, location='files')
        args = parser.parse_args()

        filename = resultfiles.save(args['file'], folder=str(result.id))

        result.filename = filename
        result.save()

        postprocess_result_file.delay(id)

        # use a raw werkzeug Response object to return 201 status code without a body
        return Response(status=201)

    def get(self, id):
        result = Result.get(Result.id == id)

        if result.filename is None:
            abort(400, message="No data available")

        # return the file
        with bz2.BZ2File(resultfiles.path(result.filename)) as infile:
            response = make_response(infile.read())
            response.headers["Content-Disposition"] = "attachment; filename={:}".format(path.splitext(result.filename)[0])
            response.headers["Content-Type"] = "text/plain"
            return response


class ResultActionResource(Resource):
    def get(self, rid, action, tid):
        Result.get(Result.id == rid)

        if action == 'doPostprocessing':
            async_result = postprocess_result_file.AsyncResult(tid)

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
        Result.get(Result.id == rid)

        parser = reqparse.RequestParser()
        parser.add_argument('doPostprocessing', type=dict, required=False)
        args = parser.parse_args()

        if args['doPostprocessing'] is not None:
            if ('update' in args['doPostprocessing'].keys() and not
               type(args['doPostprocessing']['update']) == bool):
                raise ParameterError("the 'update' parameter must be a boolean")

            update = args['doPostprocessing'].get('update', False)
            async_result = postprocess_result_file.delay(rid, update)

            return Response(status=202, headers={'Location': api.url_for(
                        ResultActionResource,
                        rid=rid, action='doPostprocessing',
                        tid=async_result.id)})

        raise ParameterError("invalid action specified: {}".format(action))


class ResultsActionResource(Resource):
    def get(self, action, tid):

        if action == 'doPostprocessing':
            group_result = postprocess_result_files.AsyncResult(tid)

            if not group_result.ready():
                return Response(status=202, headers={'Location': request.path})

            rids, results = group_result.result
            all_done = all(r.ready() for r in results)

            resp = [{'status': res.status,
                     '_links': {'parent': api.url_for(ResultResource, id=rid)},
                     'result': {'updated': res.result} if res.successful() else {}}
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
                raise ParameterError("the 'update' parameter must be a boolean")

            update = args['doPostprocessing'].get('update', False)

            async_result = postprocess_result_files.delay(update)

            return Response(status=202, headers={'Location': api.url_for(
                        ResultsActionResource, action='doPostprocessing',
                        tid=async_result.id)})

        raise ParameterError("invalid action specified: {}".format(action))


class ResultList(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('test', type=str)
        parser.add_argument('method', type=int)
        parser.add_argument('structure', type=str)
        parser.add_argument('calculated_on', type=str)
        parser.add_argument('code', type=str)
        args = parser.parse_args()

        q = Result.select(Result, Task, TaskStatus, Method, Structure, Test, PseudopotentialFamily, BasissetFamily) \
            .join(Task) \
               .join(TaskStatus).switch(Task) \
               .join(Method) \
                  .join(PseudopotentialFamily).switch(Method) \
                  .join(BasissetFamily).switch(Method) \
                  .switch(Task) \
               .join(Structure).switch(Task) \
               .join(Test).switch(Task) \
               .switch(Result) \
            .order_by(Result.id.asc())

        if args['test'] is not None:
            t = Test.get(Test.name == args['test'])
            q = q.where(Task.test == t)

        if args['method'] is not None:
            m = Method.get(Method.id == args['method'])
            q = q.where(Task.method == m)
        
        if args['structure']:
            q = q.where(Structure.name == args['structure'])

        if args['calculated_on']:
            q = q.where(Task.machine == args['calculated_on'])

        if args['code']:
            q = q.where(Method.code == args['code'])

        return [marshal(model_to_dict(x), result_resource_fields) for x in q]


    @marshal_with(result_resource_fields)
    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('energy', type=float, required=True)
        parser.add_argument('task', type=str)
        parser.add_argument('task_id', type=int)
        parser.add_argument('data', type=str)
        args = parser.parse_args()

        if args['task'] is not None:
            url, data = route_from(args['task'], 'GET')
            if url != 'taskresource':
                abort(400, message="Invalid URL specified for task")
            task_id = data['id']
        else:
            task_id = args['task_id']

        if task_id is None:
            abort(400, message="Either task or task_id must be specified")

        if args['data'] is not None:
            extradata = json.loads(args['data'])
        else:
            extradata = None

        result = Result(task=Task.get(Task.id==task_id), energy=args['energy'], data=extradata)
        result.save()

        return model_to_dict(result), 201, {'Location': api.url_for(ResultResource, id=result.id)}

class TestResource(Resource):
    @marshal_with(test_resource_fields)
    def get(self, testname):
        st = TestStructure.select(TestStructure, Test, Structure) \
                .join(Test).switch(TestStructure) \
                .join(Structure).switch(Structure) \
                .where(Test.name == testname)

        return [marshal(model_to_dict(s), test_resource_fields) for s in st]

class TestList(Resource):
    def get(self):
        return [(t.id, str(t)) for t in Test.select()]

class MethodList(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('test', type=str)
        args = parser.parse_args()

        q = Method.select(Method, PseudopotentialFamily, BasissetFamily) \
            .join(PseudopotentialFamily).switch(Method) \
            .join(BasissetFamily).switch(Method)

        if args['test']:
            q = q.join(TestResult) \
                    .join(Test).switch(Method) \
                .where(Test.name == args['test']) \
                .order_by(TestResult.method)
        else:
            q = q.order_by(Method.id)

        return [marshal(model_to_dict(m), method_list_fields) for m in q]

    @marshal_with(method_resource_fields)
    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('code', type=str, required=True)
        # TODO: the column and the attribute are called pseudopotential, but are in fact family-references
        parser.add_argument('pseudopotential', type=str, required=True)
        parser.add_argument('basis_set', type=int, required=True)
        parser.add_argument('settings', type=str, required=True)
        args = parser.parse_args()

        b = BasissetFamily.get(BasissetFamily.id == args['basis_set'])
        p = PseudopotentialFamily.get(PseudopotentialFamily.name == args['pseudopotential'])

        # the client relies on this method doing deduplication by using get_or_create
        m, _ = Method.get_or_create(code=args['code'],
                                    basis_set=b, pseudopotential=p,
                                    settings=json.loads(args['settings']))

        return model_to_dict(m)

class MethodResource(Resource):
    """Return all the details for a method, or set them"""
    @marshal_with(method_resource_fields)
    def get(self, id):
        q = Method.select(Method, PseudopotentialFamily, BasissetFamily) \
            .where(Method.id == id) \
            .join(PseudopotentialFamily).switch(Method) \
            .join(BasissetFamily).switch(Method) \
            .get()

        return model_to_dict(q)

class Basissets(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('family', type=str, required=True)
        parser.add_argument('element', type=str, required=True, action="append")
        args = parser.parse_args()

        q = BasisSet \
                .select(BasisSet, BasissetFamily) \
                .join(BasissetFamily) \
                .where((BasissetFamily.name == args['family']) & (BasisSet.element << args['element']))

        return {b.element: b.basis for b in q}

def pseudo_to_dict(pseudo):
    '''Convert an empty dict to None to trigger the default for fields.Nested'''

    pdict = model_to_dict(pseudo)

    if not pdict['converted_from']:
        pdict['converted_from'] = None

    return pdict

class PseudopotentialResource(Resource):
    @marshal_with(pseudo_resource_fields)
    def get(self, id):
        ConvertedPseudo = Pseudopotential.alias()
        q = (Pseudopotential
             .select(Pseudopotential, PseudopotentialFamily.name, ConvertedPseudo)
             .join(PseudopotentialFamily).switch(Pseudopotential)
             .join(ConvertedPseudo, JOIN.LEFT_OUTER,
                   on=(Pseudopotential.converted_from == ConvertedPseudo.id)).switch(Pseudopotential)
             .where(Pseudopotential.id == id)
             .get())

        return pseudo_to_dict(q)

class PseudopotentialFamilyList(Resource):
    """Get the list of pseudopotential families"""
    def get(self):
        q = (PseudopotentialFamily
             .select(PseudopotentialFamily)
             .order_by(PseudopotentialFamily.name))

        return [marshal(model_to_dict(f), pseudofamily_list_fields) for f in q]

class PseudopotentialList(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('family', type=str)
        parser.add_argument('format', type=str)
        parser.add_argument('element', type=str, action="append")
        args = parser.parse_args()

        ConvertedPseudo = Pseudopotential.alias()
        q = (Pseudopotential
             .select(Pseudopotential, PseudopotentialFamily.name, ConvertedPseudo)
             .join(PseudopotentialFamily).switch(Pseudopotential)
             .join(ConvertedPseudo, JOIN.LEFT_OUTER,
                   on=(Pseudopotential.converted_from == ConvertedPseudo.id)).switch(Pseudopotential)
             .order_by(Pseudopotential.id.asc()))

        if args['family']:
            q = q.where(PseudopotentialFamily.name == args['family'])

        if args['format']:
            q = q.where(Pseudopotential.format == args['format'])

        if args['element']:
            q = q.where(Pseudopotential.element << args['element'])

        if any(args.values()):
            return [marshal(pseudo_to_dict(p), pseudo_resource_fields) for p in q]
        else:
            return [marshal(pseudo_to_dict(p), pseudo_list_fields) for p in q]

    @marshal_with(pseudo_resource_fields)
    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('family', type=str, required=True)
        parser.add_argument('element', type=str, required=True)
        parser.add_argument('pseudo', type=str, required=True)
        parser.add_argument('format', type=str, required=True)
        parser.add_argument('converted_from', type=int, required=False)
        parser.add_argument('overwrite', type=bool, default=False)

        args = parser.parse_args()

        f, _ = PseudopotentialFamily.get_or_create(name=args['family'])

        data = {'family': f,
                'element': args['element'],
                'format': args['format']}

        if args['converted_from']:
            print('here')
            data['converted_from'] = args['converted_from']

        p, created = Pseudopotential.get_or_create(defaults=dict(pseudo=args['pseudo']), **data)

        if created:
            return pseudo_to_dict(p), 201, {'Location': api.url_for(PseudopotentialResource, id=p.id)}

        # the user can specify overwriting of the pseudo
        if args['overwrite']:
            q = (Pseudopotential
                    .update(pseudo=args['pseudo'])
                    .where(Pseudopotential.id == p.id)
                    .returning(Pseudopotential)
                    .execute())
            return pseudo_to_dict(next(q)), 200, {'Location': api.url_for(PseudopotentialResource, id=p.id)}
        else:
            abort(400, message="Pseudo is already uploaded for this result")

class MachineStatus(Resource):
    """Return a dictionary with the number of running and total tasks per machine:
              MACHINE   RUN   TOTAL
       e.g.   machine1   4    124
              machine2   5    24
              machine3   0    242
     """
    def get(self):
        ret = {}

        # see https://github.com/coleifer/peewee/issues/1010
        q = Task.select(Task.machine,
                    fn.COUNT(Task.id).alias('total'),
                    fn.COUNT(SQL('CASE WHEN t2.name = %s THEN t1.id END', 'running')).alias('running')) \
                .join(TaskStatus).alias('ts').switch(Task) \
                .group_by(Task.machine)

        return [(m.machine, m.running, m.total) for m in q]

class TestResultResource(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('method', type=int, required=False)
        parser.add_argument('test', type=str, required=True)
        args = parser.parse_args()

        q = TestResult.select(Test, TestResult) \
                .join(Test) \
                .where(Test.name == args["test"])

        if args['method']:
            q = q.where(TestResult.method == args["method"])
            return [{tr.test.name: tr.result_data} for tr in q]
        else:
            return [{tr.method_id: tr.result_data} for tr in q]

class Plot(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('method', type=int, action="append")
        parser.add_argument('test', type=str, required=True)
        args = parser.parse_args()

        test1 = Test.get(Test.name==args["test"])


        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
       #from matplotlib.figure import Figure
       #from matplotlib.dates import DateFormatter
        import matplotlib.pyplot as plt
        import itertools
        colors = itertools.cycle(["#1D0662", "#8F1E00", "#00662C", "#8F7700", "#3E1BA7", "#F33D0D", "#0AAD51", "#F3CD0D", "#8E7BC5", "#FFAA93", "#75CA9A", "#FFED93"])


        fig = plt.figure(figsize=(12,8),dpi=200, facecolor="#FFFFFF")
        fig.subplots_adjust(left=0.07,right=0.98,top=0.99,bottom=0.04)
        ax = plt.subplot(111)
        stride = int(min(35./len(args['method']),7))
        label_xpos = 5
        warning_yshift = 0
        i=0

        for m in args['method']:
            i+=1
            meth = Method.get(Method.id==m)

            #here the curve, based on the fitted data
            try:
                r = TestResult.get((TestResult.method==meth) & (TestResult.test==test1))
            except TestResult.DoesNotExist:
                ax.text(0.5,0.5+warning_yshift, "No data for method {:}".format(m), transform=ax.transAxes, color="#FF0000")
                warning_yshift+=0.05
                continue

            if not (r.result_data['_status']=='fitted' or r.result_data['_status']=='reference'):
                ax.text(0.5,0.5+warning_yshift, "Fitting problem for method {:}".format(m), transform=ax.transAxes, color="#FF0000")
                warning_yshift+=0.05
                continue

            mycolor = next(colors)
            V0 = r.result_data['V']
            B0 = r.result_data['B0']
            B1 = r.result_data['B1']
            xfit,yfit = eos(V0,B0,B1)
            ax.plot(xfit,yfit, color=mycolor)


            #here the points
            q = Result.select(Result, Task, TaskStatus, Method, Structure, Test) \
                .join(Task) \
                   .join(TaskStatus).switch(Task) \
                   .join(Method).switch(Task) \
                   .join(Structure).switch(Task) \
                   .join(Test).switch(Task) \
                   .switch(Result) \
                .order_by(Result.id.asc()) \
                .where(Task.method == meth) \
                .where(Task.test == test1)

            E = []
            V = []
            for x in q:
                natom = len(Json2Atoms(x.task.structure.ase_structure).get_masses())
                V.append(Json2Atoms(x.task.structure.ase_structure).get_volume()/natom)
                E.append(x.energy/natom)
            x2 = np.array(V)
            y2 = np.array(E)

            if 'E0' in r.result_data.keys():
                E0 = r.result_data['E0']
                y2 -= E0

            ax.plot(x2,y2, 's', color=mycolor)
             

            ax.annotate("Method {:}".format(m), xy= (xfit[label_xpos], yfit[label_xpos]), xytext = (-5,-30 if i!=3 else 20), textcoords = "offset points", fontsize=10, arrowprops=dict(facecolor='black', shrink=0.05, headwidth=3, width=1), horizontalalignment='right', color=mycolor)
            label_xpos += stride

        canvas=FigureCanvas(fig)
        png_output = BytesIO()
        canvas.print_png(png_output)
        response=make_response(png_output.getvalue())
        response.headers['Content-Type'] = 'image/png'
        return response


class StructureResource(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('id', type=int, required=False)
        parser.add_argument('test', type=str, required=False)
        parser.add_argument('repeat', type=int, required=False, default=1)
        parser.add_argument('viewer', type=bool, required=False, default=False)
        parser.add_argument('size', type=int, required=False, default=300)
        parser.add_argument('format', type=str, required=False,
                            default="xyz", choices=("xyz", "json"))
        args = parser.parse_args()
    
        if args['id'] and args['test']:
            raise ParameterError("exactly one of id or test parameter is required")

        if args['id'] is not None:
            s = Structure.get(Structure.id == args['id'])
        elif args['test']:
            s = Structure.select() \
                    .join(TestStructure).join(Test) \
                    .where(Test.name == args['test']) \
                    .get()
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
                                 zip(atoms_rep.get_chemical_symbols(), atoms_rep.get_positions())
                                ])
        response = make_response(xyz_output)
        response.headers["Content-Type"] = "text/plain"

        return response

class Comparison(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('method1', type=int, required=True)
        parser.add_argument('method2', type=int, required=False)
        parser.add_argument('test', type=str, action="append")
        args = parser.parse_args()

        all_delta = []

        # COMPARE all tests for 1 reference method
        # if neither method2 nor test are not specified
        mode = "1method"

        method1 = Method.get(Method.id == args["method1"])

        if args["method2"] is not None:
            # COMPARE 2 METHODS
            mode = "2methods"
            method2 = Method.get(Method.id == args["method2"])
        elif args["method2"] is None and args["test"] is not None:
            # COMPARE all tests for 1 reference method
            mode = "methodbytest"
            if len(args["test"]) > 1:
                raise ParameterError(("Specify either 2 methods "
                                      "and optionally a list of tests, "
                                      "or 1 method and 1 test"))

        ret = {"test": {}, "methods": [], "method": {}, "summary": {}}

        if mode == "methodbytest":
            ret["methods"] = [method1.id]

            testname = args["test"][0]

            for method2 in Method.select():
                result_data = self._getResultData(method1.id, method2.id,
                                                  testname)
                if result_data:
                    ret["method"][method2.id] = [str(method2)] + result_data

        elif mode == "2methods":
            if args["test"] is not None:
                testlist = args["test"]
            else:
                testlist = [t.name for t in Test.select()]

            for testname in testlist:
                result_data = self._getResultData(method1.id, method2.id,
                                                  testname)
                if result_data:
                    ret["test"][testname] = result_data
                    all_delta.append(result_data[-1])

            all_delta = np.array(all_delta)
            ret["summary"] = {"avg": np.average(all_delta),
                              "stdev": np.std(all_delta),
                              "N": len(all_delta)}
            ret["methods"] = [method1.id, method2.id]

        elif mode == "1method":
            return self._get_1method(method1.id)

        return ret

    @staticmethod
    @cache.memoize(timeout=3600)
    def _get_1method(method1_id):
        """We need the ID here and not the object to get a proper cache key"""

        ret = {"test": {}, "methods": [], "method": {}, "summary": {}}

        testlist = [t.name for t in Test.select()]

        # loop over method2
        q2 = Method.select().order_by(Method.id)
        for method2 in q2:
            ret["methods"].append(method2.id)
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
            r1 = (TestResult.select(TestResult)
                  .join(Test)
                  .where((TestResult.method == method1_id) &
                         (Test.name == testname))
                  .get())
            r2 = (TestResult.select(TestResult)
                  .join(Test)
                  .where((TestResult.method == method2_id) &
                         (Test.name == testname))
                  .get())
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
            q = (TaskStatus
                 .select(TaskStatus.name, fn.COUNT(Task.id))
                 .join(Task, JOIN.LEFT_OUTER)
                 .group_by(TaskStatus.name)
                 .order_by(TaskStatus.name))

            stats = [{"status": e.name,
                      "count": e.count,
                      "_links": {"self": url_for("tasklist", status=e.name)}
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

api.add_resource(TaskResource, '/tasks/<int:id>')
api.add_resource(TaskList, '/tasks')
api.add_resource(ResultResource, '/results/<int:id>')
api.add_resource(ResultsActionResource,
                 '/results/action/<string:action>:<string:tid>')
api.add_resource(ResultsActionList, '/results/action')
api.add_resource(ResultFileResource, '/results/<int:id>/file')
api.add_resource(ResultActionList,
                 '/results/<int:rid>/action')
api.add_resource(ResultActionResource,
                 '/results/<int:rid>/action/<string:action>:<string:tid>')
api.add_resource(ResultList, '/results')
api.add_resource(Basissets, '/basis')
api.add_resource(PseudopotentialResource, '/pseudos/<int:id>')
api.add_resource(PseudopotentialList, '/pseudos')
api.add_resource(PseudopotentialFamilyList, '/pseudofamilies')
api.add_resource(MachineStatus, '/machinestatus')
api.add_resource(TestResultResource, '/testresult')
api.add_resource(Comparison, '/compare')
api.add_resource(MethodResource, '/methods/<int:id>')
api.add_resource(MethodList, '/methods')
api.add_resource(TestList, '/tests')
api.add_resource(TestResource, '/tests/<string:testname>')
api.add_resource(Plot, '/plot')
api.add_resource(StructureResource, '/structure')
api.add_resource(StatsResource, '/stats/<string:what>')

#  vim: set ts=4 sw=4 tw=0 :
