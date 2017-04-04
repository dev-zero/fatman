
from flask import redirect, abort, url_for, request
from flask_admin import Admin
from flask_admin import helpers as admin_helpers
from flask_admin.contrib.sqla import ModelView
from flask_security import current_user

from . import security, app, db
from .models import (
    User,
    Role,
    Structure,
    StructureSet,
    BasisSet,
    BasisSetFamily,
    Pseudopotential,
    PseudopotentialFamily,
    Method,
    Test,
    Task,
    Result,
    TestResult,
    Calculation,
    CalculationCollection,
    CalculationDefaultSettings,
    Task2,
    TestResult2,
    TestResult2Collection,
    Code,
    Machine,
    Artifact,
    Command,
    TaskRuntimeSettings
)


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


class UserView(BaseManagementView):
    column_list = ('email', 'active', 'confirmed_at')


class StructureSetView(BaseDataView):
    column_list = ('id', 'name', 'description', 'superset', )
    form_columns = ('name', 'description', 'superset', 'structures', 'subsets', )
    column_details_list = ('id', ) + form_columns


class StructureView(BaseDataView):
    column_list = ('id', 'name')
    column_searchable_list = ('name',)
    column_filters = ('name',)


class StructureSetStructureView(BaseDataView):
    column_display_pk = False
    column_list = ('set', 'structure')


class MethodView(BaseDataView):
    column_list = (
        'id',
        'code',
        'pseudopotential',
        'basis_set',
        'settings',
        )
    column_filters = ('code',)


class BasisSetView(BaseDataView):
    column_list = ('id', 'element', 'family', 'basis')
    column_formatters = {'basis': lambda v, c, m, p: m.basis[:80]+"..."}
    column_filters = ('element', 'family',)


class BasisSetFamilyView(BaseDataView):
    column_display_pk = False
    column_list = ('name',)


class PseudopotentialView(BaseDataView):
    column_list = ('id', 'element', 'family', 'format',
                   'pseudo', 'core_electrons', 'converted_from', )
    column_formatters = {'pseudo': lambda v, c, m, p: m.pseudo[:80]+"..."}
    column_filters = ('element', 'family', 'format',)


class PseudopotentialFamilyView(BaseDataView):
    column_display_pk = False
    column_list = ('name',)


class TestView(BaseDataView):
    column_display_pk = False
    column_list = ('name', 'description')


class TestStructureView(BaseDataView):
    column_display_pk = False
    column_list = ('test', 'structure')


class TaskView(BaseDataView):
    column_list = ('id', 'status', 'test', 'structure', 'method',
                   'machine', 'priority', 'ctime', 'mtime',)
    column_filters = ('status',)


class ResultView(BaseDataView):
    column_list = ('id', 'task', 'filename', 'data',)


class TestResultView(BaseDataView):
    column_list = ('id', 'test', 'method', 'ctime', 'result_data',)


class CalculationView(BaseDataView):
    column_list = ('id', 'collection', 'test', 'structure',
                   'settings', 'restrictions', 'results')
    column_filters = ('test', 'collection', 'structure.name', )


class CalculationCollectionView(BaseDataView):
    column_list = ('id', 'name', 'desc', )
    form_columns = ('name', 'desc', )


class CalculationDefaultSettingsView(BaseDataView):
    column_list = ('code', 'test', 'structure', 'structure_set', )


class Task2View(BaseDataView):
    column_list = ('id', 'calculation', 'status', 'restrictions',
                   'machine', 'priority', 'ctime', 'mtime',)
    column_filters = ('status', 'machine', )


class TestResult2View(BaseDataView):
    column_list = ('id', 'test', 'calculations', 'data', )


class TestResult2CollectionView(BaseDataView):
    column_list = ('id', 'name', 'desc', )


class CodeView(BaseDataView):
    column_list = ('id', 'name', 'pseudo_format', )
    form_columns = ('id', 'name', 'pseudo_format', )


class MachineView(BaseDataView):
    column_list = ('id', 'shortname', 'name', )
    form_columns = ('shortname', 'name', 'settings', 'commands', )
    inline_models = [
        (TaskRuntimeSettings, dict(form_columns=['id', 'code', 'test', 'settings'])),
        ]


class CommandView(BaseDataView):
    column_exclude_list = ('environment', 'commands', )


class ArtifactView(BaseDataView):
    column_list = ('id', 'name', 'path')


admin = Admin(app,
              name='FATMAN', template_mode='bootstrap3',
              base_template='my_master.html')

admin.add_view(UserView(User, db.session))
admin.add_view(BaseManagementView(Role, db.session))

admin.add_view(StructureView(Structure, db.session))
admin.add_view(StructureSetView(StructureSet, db.session))
admin.add_view(MethodView(Method, db.session))
admin.add_view(BasisSetView(BasisSet, db.session))
admin.add_view(BasisSetFamilyView(BasisSetFamily, db.session))
admin.add_view(PseudopotentialView(Pseudopotential, db.session))
admin.add_view(PseudopotentialFamilyView(PseudopotentialFamily, db.session))
admin.add_view(TestView(Test, db.session))
admin.add_view(TaskView(Task, db.session))
admin.add_view(ResultView(Result, db.session))
admin.add_view(TestResultView(TestResult, db.session))

admin.add_view(CalculationView(Calculation, db.session))
admin.add_view(CalculationCollectionView(CalculationCollection, db.session))
admin.add_view(CalculationDefaultSettingsView(CalculationDefaultSettings, db.session))
admin.add_view(Task2View(Task2, db.session))
admin.add_view(TestResult2View(TestResult2, db.session))
admin.add_view(TestResult2CollectionView(TestResult2Collection, db.session))
admin.add_view(CodeView(Code, db.session))
admin.add_view(CommandView(Command, db.session))
admin.add_view(MachineView(Machine, db.session))
admin.add_view(ArtifactView(Artifact, db.session))


# define a context processor for merging flask-admin's template
# context into the flask-security views.
@security.context_processor
def security_context_processor():
    return dict(
        admin_base_template=admin.base_template,
        admin_view=admin.index_view,
        h=admin_helpers,
        get_url=url_for
    )

#  vim: set ts=4 sw=4 tw=0 :
