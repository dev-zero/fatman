
from flask import redirect, abort, url_for, request
from flask_admin import Admin
from flask_admin import helpers as admin_helpers
from flask_admin.contrib.peewee import ModelView
from flask_security import current_user

# Wire-up JSONB fields from peewee-playhouse to WTForms
from flask_admin.contrib.peewee.form import CustomModelConverter
from playhouse.postgres_ext import BinaryJSONField
from wtforms import fields

from . import security, app
from .models import *

# Make flask_admin render the BinaryJSONFields as TextAreas
CustomModelConverter.defaults[BinaryJSONField] = fields.TextAreaField


class BaseDataView(ModelView):
    """
    Base View for non-data models
    """
    # I want the id to be displayed by default
    column_display_pk = True
    # permit CSV export everywhere
    can_export = True
    page_size = 150
    can_view_details = True

    def is_accessible(self):
        """
        Make them only accessible to active authenticated users
        """
        return current_user.is_active and current_user.is_authenticated

    def _handle_view(self, name, **kwargs):
        """
        Override builtin _handle_view in order to redirect users
        when a view is not accessible.
        """
        if not self.is_accessible():
            if current_user.is_authenticated:
                # permission denied
                abort(403)
            else:
                # login
                return redirect(url_for('security.login', next=request.url))


class BaseManagementView(ModelView):
    """
    Base View for non-data models
    """
    def is_accessible(self):
        """
        Make them only accessible to users with the admin role
        """
        return (current_user.is_active and
                current_user.is_authenticated and
                current_user.has_role('admin'))

    def _handle_view(self, name, **kwargs):
        """
        Override builtin _handle_view in order to redirect users
        when a view is not accessible.
        """
        if not self.is_accessible():
            if current_user.is_authenticated:
                # permission denied
                abort(403)
            else:
                # login
                return redirect(url_for('security.login', next=request.url))


class StructureSetView(BaseDataView):
    column_list = ('id', 'name', 'description')


class StructureView(BaseDataView):
    column_list = ('id', 'name', 'ase_structure')


class StructureSetStructureView(BaseDataView):
    column_display_pk = False
    column_list = ('set', 'structure')


class MethodView(BaseDataView):
    column_list = ('id', 'code', 'pseudopotential', 'basis_set', 'settings')


class BasisSetView(BaseDataView):
    column_list = ('id', 'element', 'family', 'basis')
    column_formatters = {'basis': lambda v, c, m, p: m.basis[:80]+"..."}


class BasissetFamilyView(BaseDataView):
    column_display_pk = False
    column_list = ('name',)


class PseudopotentialView(BaseDataView):
    column_list = ('id', 'element', 'family', 'format',
                   'pseudo', 'converted_from', )
    column_formatters = {'pseudo': lambda v, c, m, p: m.pseudo[:80]+"..."}


class PseudopotentialFamilyView(BaseDataView):
    column_display_pk = False
    column_list = ('name',)


class TestView(BaseDataView):
    column_list = ('id', 'name', 'description')


class TestStructureView(BaseDataView):
    column_display_pk = False
    column_list = ('test', 'structure')


class TaskView(BaseDataView):
    column_list = ('id', 'status', 'test', 'structure', 'method',
                   'machine', 'priority', 'ctime', 'mtime',)


class ResultView(BaseDataView):
    column_list = ('id', 'task', 'energy', 'filename', 'data',)


class TestResultView(BaseDataView):
    column_list = ('id', 'test', 'method', 'ctime', 'result_data',)


admin = Admin(app,
              name='FATMAN', template_mode='bootstrap3',
              base_template='my_master.html',
              )

admin.add_view(BaseManagementView(User))
admin.add_view(BaseManagementView(Role))
admin.add_view(BaseManagementView(UserRole))

admin.add_view(StructureView(Structure))
admin.add_view(StructureSetView(StructureSet))
admin.add_view(StructureSetStructureView(StructureSetStructure))
admin.add_view(MethodView(Method))
admin.add_view(BasisSetView(BasisSet))
admin.add_view(BasissetFamilyView(BasissetFamily))
admin.add_view(PseudopotentialView(Pseudopotential))
admin.add_view(PseudopotentialFamilyView(PseudopotentialFamily))
admin.add_view(TestView(Test))
admin.add_view(TestStructureView(TestStructure))
admin.add_view(TaskView(Task))
admin.add_view(ResultView(Result))
admin.add_view(ResultView(ResultWithoutTestResult))
admin.add_view(TestResultView(TestResult))


# define a context processor for merging flask-admin's template
# context into the flask-security views.
@security.context_processor
def security_context_processor():
    return dict(
        admin_base_template=admin.base_template,
        admin_view=admin.index_view,
        h=admin_helpers,
    )

#  vim: set ts=4 sw=4 tw=0 :
