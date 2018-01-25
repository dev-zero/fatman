"""make default_settings priority based with more conditions

Revision ID: 554eb9aa8c51
Revises: 7e123be69564
Create Date: 2018-01-25 14:02:10.891044

"""

# revision identifiers, used by Alembic.
revision = '554eb9aa8c51'
down_revision = '7e123be69564'

import uuid
from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql
from sqlalchemy.sql import text, column

def upgrade():
    op.add_column('calculation_default_settings', sa.Column('basis_set_family_id', sa.Integer(), nullable=True))
    op.add_column('calculation_default_settings', sa.Column('basis_set_id', postgresql.UUID(as_uuid=True), nullable=True))
    op.add_column('calculation_default_settings', sa.Column('priority', sa.Integer(), server_default='0', nullable=False))
    op.add_column('calculation_default_settings', sa.Column('pseudopotential_family_id', sa.Integer(), nullable=True))
    op.add_column('calculation_default_settings', sa.Column('pseudopotential_id', postgresql.UUID(as_uuid=True), nullable=True))
    op.create_foreign_key(None, 'calculation_default_settings', 'pseudopotential', ['pseudopotential_id'], ['id'])
    op.create_foreign_key(None, 'calculation_default_settings', 'basis_set_family', ['basis_set_family_id'], ['id'])
    op.create_foreign_key(None, 'calculation_default_settings', 'basis_set', ['basis_set_id'], ['id'])
    op.create_foreign_key(None, 'calculation_default_settings', 'pseudopotential_family', ['pseudopotential_family_id'], ['id'])

    # create the "not-UUID-0", resp. "not-0" constraints for basissets and pseudos FKs:
    op.create_check_constraint(None, 'calculation_default_settings', column('basis_set_id') != str(uuid.UUID(int=0)))
    op.create_check_constraint(None, 'calculation_default_settings', column('basis_set_family_id') != 0)
    op.create_check_constraint(None, 'calculation_default_settings', column('pseudopotential_id') != str(uuid.UUID(int=0)))
    op.create_check_constraint(None, 'calculation_default_settings', column('pseudopotential_family_id') != 0)

    # drop and re-create the special unique constraint:
    op.drop_index('unique_calculation_default_settings_idx', 'calculation_default_settings')
    op.create_index(
        'unique_calculation_default_settings_idx',
        'calculation_default_settings',
        [
            'code_id',
            'test_id',
            text("coalesce(structure_id, '00000000-0000-0000-0000-000000000000')"),
            text("coalesce(structure_set_id, 0)"),
            text("coalesce(basis_set_id, '00000000-0000-0000-0000-000000000000')"),
            text("coalesce(basis_set_family_id, 0)"),
            text("coalesce(pseudopotential_id, '00000000-0000-0000-0000-000000000000')"),
            text("coalesce(pseudopotential_family_id, 0)"),
            ]
        )


def downgrade():
    # drop and re-create the special unique constraint:
    op.drop_index('unique_calculation_default_settings_idx', 'calculation_default_settings')

    op.drop_constraint('calculation_default_settings_basis_set_id_check', 'calculation_default_settings')
    op.drop_constraint('calculation_default_settings_basis_set_family_id_check', 'calculation_default_settings')
    op.drop_constraint('calculation_default_settings_pseudopotential_id_check', 'calculation_default_settings')
    op.drop_constraint('calculation_default_settings_pseudopotential_family_id_check', 'calculation_default_settings')

    op.drop_constraint('calculation_default_settings_basis_set_id_fkey', 'calculation_default_settings', type_='foreignkey')
    op.drop_constraint('calculation_default_settings_basis_set_family_id_fkey', 'calculation_default_settings', type_='foreignkey')
    op.drop_constraint('calculation_default_settings_pseudopotential_id_fkey', 'calculation_default_settings', type_='foreignkey')
    op.drop_constraint('calculation_default_settings_pseudopotential_family_id_fkey', 'calculation_default_settings', type_='foreignkey')
    op.drop_column('calculation_default_settings', 'pseudopotential_id')
    op.drop_column('calculation_default_settings', 'pseudopotential_family_id')
    op.drop_column('calculation_default_settings', 'priority')
    op.drop_column('calculation_default_settings', 'basis_set_id')
    op.drop_column('calculation_default_settings', 'basis_set_family_id')

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
