
from flask_admin import Admin, AdminIndexView
from flask_admin.contrib.peewee import ModelView
from flask.ext.security import current_user

# Wire-up JSONB fields from peewee-playhouse to WTForms
from flask_admin.contrib.peewee.form import CustomModelConverter
from playhouse.postgres_ext import BinaryJSONField
from wtforms import fields
CustomModelConverter.defaults[BinaryJSONField] = fields.TextAreaField

from fatman import app
from fatman.models import *

admin = Admin(app,
              name='FATMAN', template_mode='bootstrap3',
              index_view=AdminIndexView(url='/admin')
             )

class BaseDataView(ModelView):
    # I want the id to be displayed by default
    column_display_pk = True
    # permit CSV export everywhere
    can_export = True

class BaseManagementView(ModelView):
    """
    Base View for non-data models
    """
    def is_accessible(self):
        """
        Make them only accessible to users with the admin role
        """
        return current_user.has_role('admin')

admin.add_view(BaseManagementView(User))
admin.add_view(BaseManagementView(Role))
admin.add_view(BaseManagementView(UserRole))

admin.add_view(BaseDataView(Structure))
admin.add_view(BaseDataView(Method))
admin.add_view(BaseDataView(BasisSet))
admin.add_view(BaseDataView(BasissetFamily))
admin.add_view(BaseDataView(Pseudopotential))
admin.add_view(BaseDataView(PseudopotentialFamily))
admin.add_view(BaseDataView(Test))
admin.add_view(BaseDataView(TestStructure))
admin.add_view(BaseDataView(Task))
admin.add_view(BaseDataView(Result))
admin.add_view(BaseDataView(TestResult))
