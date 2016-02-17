
from flask_restful import Api, Resource, abort, reqparse, fields, marshal_with, marshal
from playhouse.shortcuts import model_to_dict
from datetime import datetime

from fatman import app
from fatman.models import *

task_resource_fields = {
        'id': fields.Raw,
        'ctime': fields.String,
        'mtime': fields.String,
        'machine': fields.Raw,
        'status': fields.String,
        'method': fields.String,
        'structure': fields.String,
        }

class TaskResource(Resource):
    @marshal_with(task_resource_fields)
    def get(self, task_id):
        return Task.get(Task.id == task_id)

    @marshal_with(task_resource_fields)
    def patch(self, task_id):
        parser = reqparse.RequestParser()
        parser.add_argument('status', type=str, required=True)
        args = parser.parse_args()

        task = Task.get(Task.id == task_id)

        try:
            status_id = TaskStatus.get(TaskStatus.name == args['status']).id
        except DoesNotExist:
            abort(404, message="Invalid status {}".format(args['status']))

        # update the status and reset the modification time
        task.status = status_id
        task.mtime = datetime.now()
        task.save()

        return task


class TaskList(Resource):
    def get(self):
        return [
            marshal(t, task_resource_fields)
            for t in Task.select()
            .join(TaskStatus).switch(Task)
            .join(Method).switch(Task)
            .join(Structure).switch(Task)
            .order_by(Task.id.desc())
            ]


result_resource_fields = {
        'id': fields.Raw,
        'energy': fields.Float,
        'task': fields.String,
        }

class ResultResource(Resource):
    @marshal_with(result_resource_fields)
    def get(self, result_id):
        return Result.get(Result.id == result_id)

class ResultList(Resource):
    def get(self):
        return [marshal(r, result_resource_fields) for r in Result.select().join(Task)]

# Catch common exceptions in the REST dispatcher
errors = {
        'DoesNotExist': {
            'message': "A resource with that ID does not exist.",
            'status': 404,
            },
        }

api = Api(app, prefix=app.config['APPLICATION_ROOT'], errors=errors)
api.add_resource(TaskResource, '/tasks/<int:task_id>')
api.add_resource(TaskList, '/tasks')
api.add_resource(ResultResource, '/results/<int:result_id>')
api.add_resource(ResultList, '/results')
