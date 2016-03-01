
from flask import redirect, abort, url_for, request
from flask.ext.admin import Admin
from flask.ext.admin import helpers as admin_helpers
from flask.ext.admin.contrib.peewee import ModelView
from flask.ext.security import current_user

# Wire-up JSONB fields from peewee-playhouse to WTForms
from flask.ext.admin.contrib.peewee.form import CustomModelConverter
from playhouse.postgres_ext import BinaryJSONField
from wtforms import fields
CustomModelConverter.defaults[BinaryJSONField] = fields.TextAreaField

from fatman import app, security
from fatman.models import *

admin = Admin(app,
              name='FATMAN', template_mode='bootstrap3',
              base_template='my_master.html',
             )

class BaseDataView(ModelView):
    """
    Base View for non-data models
    """
    # I want the id to be displayed by default
    column_display_pk = True
    # permit CSV export everywhere
    can_export = True

    def is_accessible(self):
        """
        Make them only accessible to active authenticated users
        """
        return current_user.is_active and current_user.is_authenticated

    def _handle_view(self, name, **kwargs):
        """
        Override builtin _handle_view in order to redirect users when a view is not accessible.
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
        return current_user.is_active and current_user.is_authenticated \
                and current_user.has_role('admin')

    def _handle_view(self, name, **kwargs):
        """
        Override builtin _handle_view in order to redirect users when a view is not accessible.
        """
        if not self.is_accessible():
            if current_user.is_authenticated:
                # permission denied
                abort(403)
            else:
                # login
                return redirect(url_for('security.login', next=request.url))

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

# define a context processor for merging flask-admin's template context into the
# flask-security views.
@security.context_processor
def security_context_processor():
    return dict(
        admin_base_template=admin.base_template,
        admin_view=admin.index_view,
        h=admin_helpers,
    )
