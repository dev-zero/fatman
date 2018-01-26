"""augmented basis sets support and desc for BasisSets/Pseudo families

Revision ID: cbbe8d742100
Revises: 554eb9aa8c51
Create Date: 2018-01-26 14:49:42.249302

"""

# revision identifiers, used by Alembic.
revision = 'cbbe8d742100'
down_revision = '554eb9aa8c51'

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

def upgrade():
    op.add_column('basis_set', sa.Column('augmented_basis_set_id', postgresql.UUID(as_uuid=True), nullable=True))
    op.create_foreign_key(None, 'basis_set', 'basis_set', ['augmented_basis_set_id'], ['id'])
    op.add_column('basis_set_family', sa.Column('augmented_family_id', sa.Integer(), nullable=True))
    op.add_column('basis_set_family', sa.Column('description', sa.Text(), nullable=True))
    op.create_foreign_key(None, 'basis_set_family', 'basis_set_family', ['augmented_family_id'], ['id'])
    op.drop_column('pseudopotential', 'description')
    op.add_column('pseudopotential_family', sa.Column('description', sa.Text(), nullable=True))


def downgrade():
    op.drop_column('pseudopotential_family', 'description')
    op.add_column('pseudopotential', sa.Column('description', sa.TEXT(), autoincrement=False, nullable=True))
    op.drop_constraint('basis_set_family_augmented_family_id_fkey', 'basis_set_family', type_='foreignkey')
    op.drop_column('basis_set_family', 'description')
    op.drop_column('basis_set_family', 'augmented_family_id')
    op.drop_constraint('basis_set_augmented_basis_set_id_fkey', 'basis_set', type_='foreignkey')
    op.drop_column('basis_set', 'augmented_basis_set_id')
