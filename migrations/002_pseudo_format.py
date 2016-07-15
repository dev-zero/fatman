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

    migrator.add_fields(
        'pseudopotential',

        format = pw.CharField(max_length=255, null=True),
        converted_from = pw.ForeignKeyField(db_column='converted_from_id',
                                            rel_model=migrator.orm['pseudopotential'],
                                            to_field='id', null=True)
        )

    # initialize the format attribute
    migrator.sql('''
        UPDATE pseudopotential
            SET format = CASE WHEN (pseudopotentialfamily.name LIKE %s) THEN %s ELSE %s END
            FROM pseudopotentialfamily
            WHERE pseudopotential.family_id = pseudopotentialfamily.id;''',
                 '%-UPF', 'UPF', 'CP2K')

    migrator.add_not_null('pseudopotential', 'format')

    # set the converted_from_ids
    migrator.sql('''
        UPDATE pseudopotential
            SET converted_from_id=subquery.from_id
            FROM (SELECT t.id AS to_id, f.id AS from_id
                FROM pseudopotential AS f JOIN pseudopotentialfamily AS ff ON f.family_id=ff.id,
                     pseudopotential AS t JOIN pseudopotentialfamily AS tf ON t.family_id=tf.id
                WHERE (ff.name NOT LIKE %s)
                    AND (f.element = t.element)
                    AND (tf.name = ff.name || %s))
                AS subquery
            WHERE subquery.to_id = pseudopotential.id;''', '%-UPF', '-UPF')

    # make the combination (family,element,format) unique
    migrator.add_index('pseudopotential', 'family', 'element', 'format', unique=True)

    # have only one Family of pseudos, but more than one pseudo per family and element

    # 1. update the pseudo entries for that
    migrator.sql('''
        UPDATE pseudopotential
            SET family_id=s.fid
            FROM (SELECT t.id as id, t.family_id as fid
                FROM pseudopotential AS f JOIN pseudopotential AS t ON f.converted_from_id=t.id)
                AS s
            WHERE (converted_from_id IS NOT NULL)
                AND (converted_from_id = s.id);''')

    # 2. update the family in method as well
    migrator.sql('''
        UPDATE method
            SET pseudopotential_id=s.to_fid
            FROM (SELECT f.id AS from_fid, t.id AS to_fid
                FROM pseudopotentialfamily AS f,
                     pseudopotentialfamily AS t
                WHERE (f.name LIKE '%%-UPF')
                    AND (f.name = t.name || '-UPF'))
                AS s
            WHERE pseudopotential_id = s.from_fid;''')

    # 3. drop all UPF pseudo families
    migrator.sql("DELETE FROM pseudopotentialfamily WHERE name LIKE '%%-UPF';")

def rollback(migrator, database, fake=False, **kwargs):
    """Write your rollback migrations here."""

    migrator.remove_fields('pseudopotential', 'format', 'converted_from')
