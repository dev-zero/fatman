#!/usr/bin/env python

from flask_restful import Api, Resource, abort, reqparse, fields, marshal_with, marshal
from playhouse.shortcuts import model_to_dict
from werkzeug.datastructures import FileStorage
from werkzeug.wrappers import Response

from datetime import datetime

from fatman import app, resultfiles
from fatman.models import *
from fatman.utils import route_from
from fatman.tools import calcDelta

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
    'structure': fields.Nested(structure_resource_fields),
    '_links': { 'self': fields.Url('taskresource') },
    }


class TaskResource(Resource):
    @marshal_with(task_resource_fields)
    def get(self, id):
        return model_to_dict(Task.get(Task.id == id))

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

        if args['status'] is not None:
            status = TaskStatus.get(TaskStatus.name == args['status'])
            q = q.where(Task.status == status)

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

        q = Result.select() \
            .join(Task) \
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

class MethodList(Resource):
    def get(self):
        return [(m.id, str(m)) for m in Method.select()]

class TestList(Resource):
    def get(self):
        return [(t.id, str(t)) for t in Test.select()]

class Methods(Resource):
    """Return all the details for a method, or set them"""
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('id', type=int)
        args = parser.parse_args()

        m = Method.get(Method.id==args['id'])
        ret = {'id'             : m.id,
               'code'           : m.code,
               'pseudopotential': m.pseudopotential.id,
               'basis_set'      : m.basis_set.id,
               'settings'       : m.settings }
               

        return ret

    def post(self):
        parser = reqparse.RequestParser()
        parser.add_argument('code', type=str)
        parser.add_argument('pseudopotential', type=int)
        parser.add_argument('basis_set', type=int)
        parser.add_argument('settings', type=str)
        args = parser.parse_args()

        b = BasissetFamily.get(BasissetFamily.id==args['basis_set'])
        p = PseudopotentialFamily.get(PseudopotentialFamily.id==args['pseudopotential'])

        print(args)
        print(json.loads(args['settings']))
        m,created = Method.get_or_create(code = args['code'],
                                         basis_set = b,
                                         pseudopotential = p,
                                         settings = json.loads(args['settings']))

        return {'id': m.id}
        

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
        print("HERE")
        print(args)
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
        
        return ret

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
        parser.add_argument('method', type=int, required=True)
        parser.add_argument('test', type=str)
        args = parser.parse_args()

        method1 = Method.get(Method.id==args["method"])

        test1 = Test.get(Test.name==args["test"])
        r = TestResult.get((TestResult.method==method1) & (TestResult.test==test1))
        ret={r.test.name: r.result_data}

        return ret

class Comparison(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('method1', type=int, required=True)
        parser.add_argument('method2', type=int, required=True)
        parser.add_argument('test', type=str, action="append")
        args = parser.parse_args()

        method1 = Method.get(Method.id==args["method1"])
        method2 = Method.get(Method.id==args["method2"])

        ret={"test": {}, "methods": [method1.id, method2.id]}
        all_delta = []

        if args["test"] is not None:
            testlist = args["test"]
        else:
            testlist = [t.name for t in Test.select()]

        for testname in testlist:
            dontadd = False
            test = Test.get(Test.name==testname)
            print(test.name)

            try:
                r1 = TestResult.get((TestResult.method==method1) & (TestResult.test==test))
                r2 = TestResult.get((TestResult.method==method2) & (TestResult.test==test))
            except:
                continue
                #ret["test"][test.name] = "N/A"
            print(r1, r2)

            if "deltatest" in testname:
                try:
                    data_f = [r1.result_data["V"], r1.result_data["B0"], r1.result_data["B1"]]
                    data_w = [r2.result_data["V"], r2.result_data["B0"], r2.result_data["B1"]]
                
                    delta = calcDelta(data_f, data_w)
                except:
                    dontadd = True

            elif "GMTKN" in testname:
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
                ret["test"][test.name] = data_f + data_w +  [delta]
                all_delta.append(delta)

            
        all_delta = np.array(all_delta)
        ret["summary"] = {"avg": np.average(all_delta), "stdev": np.std(all_delta), "N": len(all_delta)}
        return ret

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
api.add_resource(MethodList, '/methods')
api.add_resource(Methods, '/method')
api.add_resource(TestList, '/tests')

