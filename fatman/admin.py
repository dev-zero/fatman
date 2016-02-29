
from flask_admin import Admin, AdminIndexView
from flask_admin.contrib.peewee import ModelView

# Wire-up JSONB fields from peewee-playhouse to WTForms
from flask_admin.contrib.peewee.form import CustomModelConverter
from playhouse.postgres_ext import BinaryJSONField
from wtforms import fields
CustomModelConverter.defaults[BinaryJSONField] = fields.TextAreaField

from fatman import app
from fatman.models import *

admin = Admin(app,
              name='FATMAN', template_mode='bootstrap3',
              index_view=AdminIndexView(url=app.config["APPLICATION_ROOT"] + '/admin')
             )

class BaseView(ModelView):
    # I want the id to be displayed by default
    column_display_pk = True
    # permit CSV export everywhere
    can_export = True

admin.add_view(BaseView(User))
admin.add_view(BaseView(Role))
admin.add_view(BaseView(UserRole))
admin.add_view(BaseView(Structure))
admin.add_view(BaseView(Method))
admin.add_view(BaseView(BasisSet))
admin.add_view(BaseView(BasissetFamily))
admin.add_view(BaseView(Pseudopotential))
admin.add_view(BaseView(PseudopotentialFamily))
admin.add_view(BaseView(Test))
admin.add_view(BaseView(TestStructure))
admin.add_view(BaseView(Task))
admin.add_view(BaseView(Result))
admin.add_view(BaseView(TestResult))
