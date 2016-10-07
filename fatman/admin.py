
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
    column_list = ('id', 'name', 'description')


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
