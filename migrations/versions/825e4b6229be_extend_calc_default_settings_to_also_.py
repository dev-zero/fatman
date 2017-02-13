"""extend calc default settings to also use structur ids as key

Revision ID: 825e4b6229be
Revises: 77dda2cb3acf
Create Date: 2017-02-10

"""

# revision identifiers, used by Alembic.
revision = '825e4b6229be'
down_revision = '77dda2cb3acf'

import uuid
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql
from sqlalchemy.sql import column, func, text


def upgrade():
    # the sub/super-set functionality in the structure set:
    op.add_column('structure_set', sa.Column('superset_id', sa.Integer(), nullable=True))
    op.create_foreign_key(None, 'structure_set', 'structure_set', ['superset_id'], ['id'])

    # convert the PK for the settings from a composite to a single-column using an ID:
    op.add_column(
        'calculation_default_settings',
        sa.Column('id', postgresql.UUID(as_uuid=True), server_default=sa.text('gen_random_uuid()'), nullable=False))

    # drop the current primary key and create a new one
    op.drop_constraint('calculation_default_settings_pkey', 'calculation_default_settings')
    op.create_primary_key(None, 'calculation_default_settings', ['id',])


    # additional optional columns for specifiying settings:
    op.add_column(
        'calculation_default_settings',
        sa.Column('structure_id', postgresql.UUID(as_uuid=True), nullable=True))
    op.create_foreign_key(None, 'calculation_default_settings', 'structure', ['structure_id'], ['id'])

    op.add_column(
        'calculation_default_settings',
        sa.Column('structure_set_id', sa.Integer(), nullable=True))
    op.create_foreign_key(None, 'calculation_default_settings', 'structure_set', ['structure_set_id'], ['id'])

    # create the "not-UUID-0", resp. "not-0" constraints for structure and structure set FKs:
    op.create_check_constraint(None, 'calculation_default_settings', column('structure_id') != str(uuid.UUID(int=0)))
    op.create_check_constraint(None, 'calculation_default_settings', column('structure_set_id') != 0)

    # create the special unique constraint:
    op.create_index(
        'unique_calculation_default_settings_idx',
        'calculation_default_settings',
        [
            'code_id',
            'test_id',
            text("coalesce(structure_id, '00000000-0000-0000-0000-000000000000')"),
            text("coalesce(structure_set_id, 0)"),
            ]
        )

    # and the constraint that not both structure and structure set can be set at the same time:
    op.create_check_constraint(
        None,
        'calculation_default_settings',
        sa.or_(column('structure_set_id') == None, column('structure_id') == None)
        )


def downgrade():
    op.drop_constraint('calculation_default_settings_check', 'calculation_default_settings')

    op.drop_index('unique_calculation_default_settings_idx', 'calculation_default_settings')
    op.drop_constraint('calculation_default_settings_structure_id_check', 'calculation_default_settings')
    op.drop_constraint('calculation_default_settings_structure_set_id_check', 'calculation_default_settings')

    op.drop_constraint('calculation_default_settings_structure_set_id_fkey',
                       'calculation_default_settings', type_='foreignkey')
    op.drop_column('calculation_default_settings', 'structure_set_id')

    op.drop_constraint('calculation_default_settings_structure_id_fkey',
                       'calculation_default_settings', type_='foreignkey')
    op.drop_column('calculation_default_settings', 'structure_id')

    # drop the single column id and recreate the combined primary key
    op.drop_column('calculation_default_settings', 'id')
    op.create_primary_key(None, 'calculation_default_settings', ['code_id', 'test_id',])

    # drop the sub/super-set columns from the structure
    op.drop_constraint('structure_set_superset_id_fkey',
                       'structure_set', type_='foreignkey')
    op.drop_column('structure_set', 'superset_id')
