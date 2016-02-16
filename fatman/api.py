
from flask_restful import Api, Resource, abort, reqparse, fields, marshal_with
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
        try:
            return Task.get(Task.id == task_id)
        except DoesNotExist:
            abort(404, message="Task {} does not exist".format(task_id))

    @marshal_with(task_resource_fields)
    def put(self, task_id):
        parser = reqparse.RequestParser()
        parser.add_argument('status', type=str, required=True)
        args = parser.parse_args()

        try:
            task = Task.get(Task.id == task_id)
        except DoesNotExist:
            abort(404, message="Task {} does not exist".format(task_id))

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
        return {
            t.id: {
                'ctime': str(t.ctime),
                'mtime': str(t.mtime),
                'machine': t.machine,
                'status': str(t.status),
                'method': str(t.method),
                'structure': str(t.structure),
                }
            for t in Task.select()
            .join(TaskStatus).switch(Task)
            .join(Method).switch(Task)
            .join(Structure).switch(Task)
            .order_by(Task.id.desc())
            }


class ResultResource(Resource):
    def get(self, result_id):
        try:
            return model_to_dict(Result.get(Result.id == result_id))
        except DoesNotExist:
            abort(404, message="Result {} does not exist".format(result_id))

class ResultList(Resource):
    def get(self):
        return [model_to_dict(r) for r in Result.select(Result.id)]

api = Api(app)
api.add_resource(TaskResource, '/tasks/<int:task_id>')
api.add_resource(TaskList, '/tasks')
api.add_resource(ResultResource, '/results/<int:result_id>')
api.add_resource(ResultList, '/results')
