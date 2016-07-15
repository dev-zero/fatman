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


def migrate(migrator, database, fake=False, **kwargs):
    """Write your migrations here."""

    @migrator.create_model
    class StructureSet(pw.Model):
        name = pw.CharField(max_length=255, unique=True)
        description = pw.TextField(null=True)

    @migrator.create_model
    class StructureSetStructure(pw.Model):
        structure = pw.ForeignKeyField(db_column='structure_id', rel_model=migrator.orm['structure'], to_field='id')
        set = pw.ForeignKeyField(db_column='set_id', rel_model=StructureSet, to_field='id')

    # see https://github.com/klen/peewee_migrate/issues/19
    migrator.add_index('structuresetstructure', 'structure', 'set', unique=True)

    migrator.sql("""
        INSERT INTO structureset (name)
            SELECT DISTINCT upper(substring(name from '([^_]+)_.+')) AS set FROM structure;""")
    migrator.sql("""
        INSERT INTO structuresetstructure (structure_id, set_id)
            SELECT structure.id AS structure_id, structureset.id AS set_id
            FROM structure, structureset
            WHERE upper(structure.name) LIKE structureset.name || '_%%';""")

def rollback(migrator, database, fake=False, **kwargs):
    """Write your rollback migrations here."""

    migrator.remove_model('structuresetstructure')
    migrator.remove_model('structureset')
