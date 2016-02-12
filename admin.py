
from flask_admin import Admin
from flask_admin.contrib.peewee import ModelView

from app import app
from models import *

admin = Admin(app, name='Pseudo-Tests DB', template_mode='bootstrap3')

class BaseView(ModelView):
    # I want the id to be displayed by default
    column_display_pk = True

admin.add_view(BaseView(Structure))
#admin.add_view(BaseView(Method))
admin.add_view(BaseView(Test))
admin.add_view(BaseView(TestStructure))
admin.add_view(BaseView(Task))
admin.add_view(BaseView(Result))
