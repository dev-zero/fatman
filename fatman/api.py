
from flask_restful import Api, Resource, abort, reqparse, fields, marshal_with, marshal
from playhouse.shortcuts import model_to_dict
from werkzeug.datastructures import FileStorage
from werkzeug.wrappers import Response

from datetime import datetime

from fatman import app, resultfiles
from fatman.models import *
from fatman.utils import route_from

method_resource_fields = {
    'id': fields.Raw,
    'code': fields.Raw,
    'pseudopotential': fields.String(attribute='pseudopotential.family.name'),
    'basis_set': fields.String(attribute='basis_set.family.name'),
    'settings': fields.Raw,
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
        args = parser.parse_args()

        # update the status and reset the modification time
        task = Task.get(Task.id == id)
        task.status = TaskStatus.get(TaskStatus.name == args['status']).id
        task.mtime = datetime.now()
        if 'machine' in args.keys():
            task.machine = args['machine'] 
        task.save()

        return model_to_dict(task)


class TaskList(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('status', type=str)
        parser.add_argument('limit', type=int)
        args = parser.parse_args()

        q = Task.select() \
            .join(TaskStatus).switch(Task) \
            .join(Method).switch(Task) \
            .join(Structure).switch(Task) \
            .order_by(Task.id.desc())

        if args['status'] is not None:
            status = TaskStatus.get(TaskStatus.name == args['status'])
            q = q.where(Task.status == status)

        if args['limit'] is not None:
            q = q.limit(args['limit'])

        return [marshal(model_to_dict(t), task_resource_fields) for t in q]


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
        return [marshal(model_to_dict(r), result_resource_fields) for r in Result.select().join(Task)]

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


class Basissets(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('family', type=str)
        parser.add_argument('element', type=str, action="append")
        args = parser.parse_args()

        ret = {}
        for element in args['element']:
            basis = BasisSet.get(BasisSet.family.name==args['family'],BasisSet.element==element)
            ret[element] = basis.basis

        return ret

class Pseudopotentials(Resource):
    def get(self):
        parser = reqparse.RequestParser()
        parser.add_argument('family', type=str)
        parser.add_argument('element', type=str, action="append")
        args = parser.parse_args()

        ret = {}
        for element in args['element']:
            pseudo = Pseudopotential.get(Pseudopotential.family==args['family'] and Pseudopotential.element==element)
            ret[element] = pseudo.pseudo

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

api = Api(app, prefix=app.config['APPLICATION_ROOT'], errors=errors)

api.add_resource(TaskResource, '/tasks/<int:id>')
api.add_resource(TaskList, '/tasks')
api.add_resource(ResultResource, '/results/<int:id>')
api.add_resource(ResultFileResource, '/results/<int:id>/file')
api.add_resource(ResultList, '/results')
api.add_resource(Basissets, '/basis')
api.add_resource(Pseudopotentials, '/pseudo')
