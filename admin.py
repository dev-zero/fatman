
from flask_admin import Admin
from flask_admin.contrib.peewee import ModelView

from app import app
from models import *

admin = Admin(app, name='Pseudo-Tests DB', template_mode='bootstrap3')

admin.add_view(ModelView(Structure))
#admin.add_view(ModelView(Method))
admin.add_view(ModelView(Test))
admin.add_view(ModelView(TestStructure))
admin.add_view(ModelView(TaskStatus))
admin.add_view(ModelView(Task))
admin.add_view(ModelView(Result))
