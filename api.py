
from flask_restful import Api, Resource, abort
from playhouse.shortcuts import model_to_dict

from app import app
from models import Task, Result

class TaskResource(Resource):
    def get(self, task_id):
        try:
            return model_to_dict(Task.get(Task.id == task_id))
        except DoesNotExist:
            abort(404, message="Task {} does not exist".format(task_id))

class TaskList(Resource):
    def get(self):
        return [model_to_dict(t) for t in Task.select(Task.id, Task.status)]

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
