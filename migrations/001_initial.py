"""Peewee migrations: ::

    > Model = migrator.orm['name']

    > migrator.sql(sql)
    > migrator.python(func, *args, **kwargs)
    > migrator.create_model(Model)
    > migrator.remove_model(Model, cascade=True)
    > migrator.add_fields(Model, **fields)
    > migrator.change_fields(Model, **fields)
    > migrator.remove_fields(Model, *field_names, cascade=True)
    > migrator.rename_field(Model, old_field_name, new_field_name)
    > migrator.rename_table(Model, new_table_name)
    > migrator.add_index(Model, *col_names, unique=False)
    > migrator.drop_index(Model, index_name)
    > migrator.add_not_null(Model, *field_names)
    > migrator.drop_not_null(Model, *field_names)
    > migrator.add_default(Model, field_name, default)

"""

import datetime as dt
import peewee as pw
import playhouse.postgres_ext as pw_pext

def migrate(migrator, database, fake=False, **kwargs):
    """Write your migrations here."""

    @migrator.create_model
    class BasissetFamily(pw.Model):
        name = pw.CharField(max_length=255, unique=True)

    @migrator.create_model
    class BasisSet(pw.Model):
        family = pw.ForeignKeyField(db_column='family_id', rel_model=BasissetFamily, to_field='id')
        element = pw.CharField(max_length=255)
        basis = pw.TextField()

    @migrator.create_model
    class PseudopotentialFamily(pw.Model):
        name = pw.CharField(max_length=255, unique=True)

    @migrator.create_model
    class Pseudopotential(pw.Model):
        family = pw.ForeignKeyField(db_column='family_id', rel_model=PseudopotentialFamily, to_field='id')
        element = pw.CharField(max_length=255)
        pseudo = pw.TextField()

    @migrator.create_model
    class Method(pw.Model):
        code = pw.CharField(max_length=255)
        pseudopotential = pw.ForeignKeyField(db_column='pseudopotential_id', rel_model=PseudopotentialFamily, to_field='id')
        basis_set = pw.ForeignKeyField(db_column='basis_set_id', rel_model=BasissetFamily, to_field='id')
        settings = pw_pext.BinaryJSONField(index=True, null=True)

    @migrator.create_model
    class Role(pw.Model):
        name = pw.CharField(max_length=255, unique=True)
        description = pw.TextField(null=True)

    @migrator.create_model
    class Structure(pw.Model):
        name = pw.CharField(max_length=255, unique=True)
        ase_structure = pw.TextField()

    @migrator.create_model
    class TaskStatus(pw.Model):
        name = pw.CharField(max_length=255, unique=True)

    @migrator.create_model
    class Test(pw.Model):
        name = pw.CharField(max_length=255, unique=True)
        description = pw.TextField(null=True)

    @migrator.create_model
    class Task(pw.Model):
        structure = pw.ForeignKeyField(db_column='structure_id', rel_model=Structure, to_field='id')
        method = pw.ForeignKeyField(db_column='method_id', rel_model=Method, to_field='id')
        status = pw.ForeignKeyField(db_column='status_id', rel_model=TaskStatus, to_field='id')
        test = pw.ForeignKeyField(db_column='test_id', rel_model=Test, to_field='id')
        ctime = pw.DateTimeField()
        mtime = pw.DateTimeField()
        machine = pw.CharField(max_length=255)
        priority = pw.IntegerField(default=0)

    @migrator.create_model
    class ResultWithoutTestResult(pw.Model):
        energy = pw.DoubleField()
        task = pw.ForeignKeyField(db_column='task_id', rel_model=Task, to_field='id')
        filename = pw.CharField(max_length=255, null=True)

    @migrator.create_model
    class Result(pw.Model):
        energy = pw.DoubleField()
        task = pw.ForeignKeyField(db_column='task_id', rel_model=Task, to_field='id')
        filename = pw.CharField(max_length=255, null=True)
        data = pw_pext.BinaryJSONField(index=True, null=True)

    @migrator.create_model
    class TestResult(pw.Model):
        ctime = pw.DateTimeField()
        test = pw.ForeignKeyField(db_column='test_id', rel_model=Test, to_field='id')
        method = pw.ForeignKeyField(db_column='method_id', rel_model=Method, to_field='id')
        result_data = pw_pext.BinaryJSONField(index=True, null=True)

    @migrator.create_model
    class TestStructure(pw.Model):
        structure = pw.ForeignKeyField(db_column='structure_id', rel_model=Structure, to_field='id')
        test = pw.ForeignKeyField(db_column='test_id', rel_model=Test, to_field='id')

    @migrator.create_model
    class User(pw.Model):
        email = pw.CharField(max_length=255, unique=True)
        password = pw.TextField()
        active = pw.BooleanField(default=False)
        confirmed_at = pw.DateTimeField(null=True)

    @migrator.create_model
    class UserRole(pw.Model):
        user = pw.ForeignKeyField(db_column='user_id', rel_model=User, to_field='id')
        role = pw.ForeignKeyField(db_column='role_id', rel_model=Role, to_field='id')



def rollback(migrator, database, fake=False, **kwargs):
    """Write your rollback migrations here."""

    migrator.remove_model('userrole')

    migrator.remove_model('user')

    migrator.remove_model('teststructure')

    migrator.remove_model('testresult')

    migrator.remove_model('result')

    migrator.remove_model('resultwithouttestresult')

    migrator.remove_model('task')

    migrator.remove_model('test')

    migrator.remove_model('taskstatus')

    migrator.remove_model('structure')

    migrator.remove_model('role')

    migrator.remove_model('method')

    migrator.remove_model('pseudopotential')

    migrator.remove_model('pseudopotentialfamily')

    migrator.remove_model('basisset')

    migrator.remove_model('basissetfamily')
