
from flask_admin import Admin
from flask_admin.contrib.peewee import ModelView

from fatman import app
from fatman.models import *

admin = Admin(app, name='FATMAN', template_mode='bootstrap3')

class BaseView(ModelView):
    # I want the id to be displayed by default
    column_display_pk = True
    # permit CSV export everywhere
    can_export = True

class MethodView(BaseView):
    # wtf-peewee does not know about BSJONField (yet)
    form_excluded_columns = ['settings']

admin.add_view(BaseView(Structure))
admin.add_view(MethodView(Method))
admin.add_view(BaseView(Test))
admin.add_view(BaseView(TestStructure))
admin.add_view(BaseView(Task))
admin.add_view(BaseView(Result))
