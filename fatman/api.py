#!/usr/bin/env python

from flask_restful import Api, Resource, abort, reqparse, fields, marshal_with, marshal
from flask import make_response
from playhouse.shortcuts import model_to_dict
from werkzeug.datastructures import FileStorage
from werkzeug.wrappers import Response

from datetime import datetime

from fatman import app, resultfiles
from fatman.models import *
from fatman.utils import route_from
from fatman.tools import calcDelta, eos, Json2Atoms

import numpy as np

import json

method_resource_fields = {
    'id': fields.Raw,
    'code': fields.Raw,
    'pseudopotential': fields.String(attribute='pseudopotential.name'),
    'basis_set': fields.String(attribute='basis_set.name'),
    'settings': fields.Raw,
    '_links': { 'self': fields.Url('methodresource') },
    }

method_list_fields = {
    'id': fields.Raw,
    'code': fields.Raw,
    'pseudopotential': fields.String(attribute='pseudopotential.name'),
    'basis_set': fields.String(attribute='basis_set.name'),
    '_links': { 'self': fields.Url('methodresource') },
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
    'id': fields.Integer,
    'ctime': fields.String,
    'mtime': fields.String,
    'machine': fields.Raw,
    'status': fields.String(attribute='status.name'),
    'method': fields.Nested(method_resource_fields),
    'structure': fields.Nested(structure_resource_fields),
    '_links': { 'self': fields.Url('taskresource') },
    'priority' : fields.Integer,
    }

task_list_fields = {
    'id': fields.Integer,
    'ctime': fields.String,
    'mtime': fields.String,
    'machine': fields.Raw,
    'status': fields.String(attribute='status.name'),
    'method': fields.Nested(method_list_fields),
    'structure': fields.Nested(structure_list_fields),
    '_links': { 'self': fields.Url('taskresource') },
    }


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
        task.mtime = datetime.now()
        if 'machine' in args.keys() and args['machine'] is not None:
            task.machine = args['machine'] 
        if 'priority' in args.keys() and args['priority'] is not None:
            task.priority = args['priority'] 
        task.save()

        return model_to_dict(task)

    def post(self, id):
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
                                             defaults = dict(ctime = datetime.now(),
                                                             mtime = datetime.now(),
                                                             status = s,
                                                             machine = '-',
                                                             priority = args['priority']))
            ret.append(ta.id)

        return ret

class TaskList(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('status', type=str)
        parser.add_argument('limit', type=int)
        parser.add_argument('timeorder', type=bool, default=False)
        parser.add_argument('machine', type=str)
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

        if args['limit'] is not None:
            q = q.limit(args['limit'])

        return [marshal(model_to_dict(t), task_list_fields) for t in q]


result_resource_fields = {
   'id': fields.Raw,
   'energy': fields.Float,
   'task': fields.Nested(task_resource_fields),
   '_links': { 'self': fields.Url('resultresource') },
   }

class ResultResource(Resource):
    @marshal_with(result_resource_fields)
    def get(self, id):
        return model_to_dict(Result.get(Result.id == id))

class ResultFileResource(Resource):
    def post(self, id):
        result = Result.get(Result.id == id)

        if result.filename is not None:
            abort(400, message="Data is already uploaded for this result")

        parser = reqparse.RequestParser()
        parser.add_argument('file', type=FileStorage, location='files', required=True)
        args = parser.parse_args()
        filename = resultfiles.save(args['file'], folder=str(result.id))

        result.filename = filename
        result.save()

        # use a raw werkzeug Response object to return 201 status code without a body
        return Response(status=201)

class ResultList(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('test', type=str)
        parser.add_argument('method', type=int)
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
        
        return [marshal(model_to_dict(x), result_resource_fields) for x in q]


    @marshal_with(result_resource_fields)
    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('energy', type=float, required=True)
        parser.add_argument('task', type=str)
        parser.add_argument('task_id', type=int)
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

        result = Result(task=Task.get(Task.id==task_id), energy=args['energy'])
        result.save()

        return model_to_dict(result), 201, {'Location': api.url_for(ResultResource, id=result.id)}

class TestList(Resource):
    def get(self):
        return [(t.id, str(t)) for t in Test.select()]

class MethodList(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('test', type=str)
        args = parser.parse_args()

        if args['test'] is not None:
            test = Test.get(Test.name==args['test'])
            tr = TestResult.select().where(TestResult.test == test)
            
            return [marshal(model_to_dict(x.method), method_list_fields) for x in tr]

        else:
            q = Method.select(Method, PseudopotentialFamily, BasissetFamily) \
                .join(PseudopotentialFamily).switch(Method) \
                .join(BasissetFamily).switch(Method) \
                .order_by(Method.id)

            return [marshal(model_to_dict(m), method_list_fields) for m in q]

    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('code', type=str, required=True)
        parser.add_argument('pseudopotential', type=int, required=True)
        parser.add_argument('basis_set', type=int, required=True)
        parser.add_argument('settings', type=str, required=True)
        args = parser.parse_args()

        b = BasissetFamily.get(BasissetFamily.id==args['basis_set'])
        p = PseudopotentialFamily.get(PseudopotentialFamily.id==args['pseudopotential'])

        m,created = Method.get_or_create(code = args['code'],
                                         basis_set = b,
                                         pseudopotential = p,
                                         settings = json.loads(args['settings']))

        return {'id': m.id}

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
        parser.add_argument('family', type=str)
        parser.add_argument('element', type=str, action="append")
        args = parser.parse_args()

        ret = {}
        for element in args['element']:
            basis = BasisSet.get(BasisSet.family==BasissetFamily.get(BasissetFamily.name==args['family']),BasisSet.element==element)
            ret[element] = basis.basis

        return ret

class Pseudopotentials(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('family', type=str)
        parser.add_argument('element', type=str, action="append")
        args = parser.parse_args()

        ret = {}
        if not args['element'] is None:
            for element in args['element']:
                f = PseudopotentialFamily.get(PseudopotentialFamily.name==args['family'])
                pseudo = Pseudopotential.get((Pseudopotential.family==f) & (Pseudopotential.element==element))
                ret[element] = pseudo.pseudo
        else:
            f = PseudopotentialFamily.get(PseudopotentialFamily.name==args['family'])
            q =  Pseudopotential.select().where(Pseudopotential.family==f)
            for ps in q:
                ret[ps.element] = ps.pseudo

        return ret

    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('family', type=str)
        parser.add_argument('element', type=str)
        parser.add_argument('pseudo', type=str)
        parser.add_argument('overwrite', type=bool, default=False)
        
        args = parser.parse_args()

        f,created = PseudopotentialFamily.get_or_create(name=args['family'])
        p,created = Pseudopotential.get_or_create(family=f, element=args['element'], defaults=dict(pseudo=args['pseudo']))

        #the user can specify overwriting of the pseudo
        if args['overwrite'] and not created:
            q = Pseudopotential.update(pseudo=args['pseudo']).where(Pseudopotential.id == p.id)
            q.execute()

        if not created and not args['overwrite']:
            abort(400, message="Pseudo is already uploaded for this result")
        
        return {'id': f.id}

class MachineStatus(Resource):
    """Return a dictionary with the number of running and total tasks per machine:
              MACHINE   RUN   TOTAL
       e.g.   machine1   4    124
              machine2   5    24
              machine3   0    242
     """
    def get(self):
        ret = {}

        q = Task.select(Task.machine,fn.Count(Task.machine).alias('count')).where(Task.status==TaskStatus.get(TaskStatus.name=="running")).group_by(Task.machine)
        for x in q:
            ret[x.machine] = [x.count]

        q = Task.select(Task.machine,fn.Count(Task.machine).alias('count')).group_by(Task.machine)
        for x in q:
            if x.machine in ret.keys():
                ret[x.machine] = [ret[x.machine][0],x.count]
            else:
                ret[x.machine] = [0,x.count]
        
        return [ [k, v[0], v[1]] for k,v in ret.items() ]

class CalcStatus(Resource):
    """Return a dictionary of task statuses and the number of tasks with that status
       e.g.   running   7
              error     2
              total    124
    """
    def get(self):
        q = TaskStatus.select().annotate(Task)
        ret = dict([(x.name,x.count) for x in q])
        
        return ret

class TestResultResource(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('method', type=int, required=False)
        parser.add_argument('test', type=str, required=True)
        args = parser.parse_args()

        test1 = Test.get(Test.name==args["test"])

        if args['method'] is not None:
            method1 = Method.get(Method.id==args["method"])

            r = TestResult.get((TestResult.method==method1) & (TestResult.test==test1))
            ret=[{r.test.name: r.result_data}]
        else:
            q = TestResult.select().where((TestResult.test==test1))
            ret=[]
            for r in q:
                ret.append({r.method_id: r.result_data})

        return ret

class Plot(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('method', type=int, action="append")
        parser.add_argument('test', type=str, required=True)
        args = parser.parse_args()

        test1 = Test.get(Test.name==args["test"])

        from io import BytesIO

        from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
       #from matplotlib.figure import Figure
       #from matplotlib.dates import DateFormatter
        import matplotlib.pyplot as plt
        import itertools
        colors = itertools.cycle(["#1D0662", "#8F1E00", "#00662C", "#8F7700", "#3E1BA7", "#F33D0D", "#0AAD51", "#F3CD0D", "#8E7BC5", "#FFAA93", "#75CA9A", "#FFED93"])


        fig = plt.figure(figsize=(12,8),dpi=200, facecolor="#FFFFFF")
        ax = plt.subplot(111)
        stride = min(35./len(args['method']),7)
        label_xpos = 5
        warning_yshift = 0

        for m in args['method']:
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
             

            ax.annotate("Method {:}".format(m), xy= (xfit[label_xpos], yfit[label_xpos]), xytext = (-5,-30), textcoords = "offset points", fontsize=10, arrowprops=dict(facecolor='black', shrink=0.05, headwidth=3, width=1), horizontalalignment='right', color=mycolor)
            label_xpos += stride

        canvas=FigureCanvas(fig)
        png_output = BytesIO()
        canvas.print_png(png_output)
        response=make_response(png_output.getvalue())
        response.headers['Content-Type'] = 'image/png'
        return response

    


class Comparison(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('method1', type=int, required=True)
        parser.add_argument('method2', type=int, required=False)
        parser.add_argument('test', type=str, action="append")
        args = parser.parse_args()

        by_test = False
        all_delta = []

        #need at least one method, which will be the 'reference' method
        method1 = Method.get(Method.id==args["method1"])

        #if a second method is specified, we compare all tests for the pair of methods
        if args["method2"] is not None:
            method2 = Method.get(Method.id==args["method2"])
        else:
            by_test = True

        ret={"test": {}, "methods": [], "method":{}}

        if by_test:
            ret["methods"] = [method1.id]
            if args["test"] is not None and len(args["test"])==1:
                testname = args["test"][0]
            else:
                return errors["ParameterError"]

            test = Test.get(Test.name==testname)

            q = Method.select()

            for method2 in q:
                result_data = self._getResultData(method1, method2, test)
                if not result_data == False:
                    ret["method"][method2.id] = [str(method2)] + result_data
                    all_delta.append(result_data[-1])

        else:
            if args["test"] is not None:
                testlist = args["test"]
            else:
                testlist = [t.name for t in Test.select()]

            for testname in testlist:
                test = Test.get(Test.name==testname)
                result_data = self._getResultData(method1, method2, test)
                if not result_data == False:
                    ret["test"][testname] = result_data
                    all_delta.append(result_data[-1])

        all_delta = np.array(all_delta)
        ret["summary"] = {"avg": np.average(all_delta), "stdev": np.std(all_delta), "N": len(all_delta)}
        return ret
    
    def _getResultData(self, method1, method2, test):
        dontadd = False
        try:
            r1 = TestResult.get((TestResult.method==method1) & (TestResult.test==test))
            r2 = TestResult.get((TestResult.method==method2) & (TestResult.test==test))
        except:
            return False
            #ret["test"][test.name] = "N/A"

        data_f = []
        data_w = []
        if "deltatest" in test.name:
            try:
                data_f = [r1.result_data["V"], r1.result_data["B0"], r1.result_data["B1"]]
                data_w = [r2.result_data["V"], r2.result_data["B0"], r2.result_data["B1"]]
            
                delta = calcDelta(data_f, data_w)
            except:
                dontadd = True

        elif "GMTKN" in test.name:
            try:
                n=0
                mad=0.
                for e1, e2 in zip(r1.result_data["energies"], r2.result_data["energies"]):
                    mad+= abs(e1-e2)
                    n+=1

                delta=mad/n
            except:
                dontadd = True

        if not dontadd:
            return data_f + data_w +  [delta]
        else:
            return False


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
api.add_resource(ResultFileResource, '/results/<int:id>/file')
api.add_resource(ResultList, '/results')
api.add_resource(Basissets, '/basis')
api.add_resource(Pseudopotentials, '/pseudo')
api.add_resource(CalcStatus, '/calcstatus')
api.add_resource(MachineStatus, '/machinestatus')
api.add_resource(TestResultResource, '/testresult')
api.add_resource(Comparison, '/compare')
api.add_resource(MethodResource, '/methods/<int:id>')
api.add_resource(MethodList, '/methods')
api.add_resource(TestList, '/tests')
api.add_resource(Plot, '/plot')

#  vim: set ts=4 sw=4 tw=0 :
