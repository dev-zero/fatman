"""models: introduce metadata fields for Calculations and TestResults

Revision ID: c4089fd3de22
Revises: fe45a2e57d78
Create Date: 2017-06-20 15:17:16.657793

"""

# revision identifiers, used by Alembic.
revision = 'c4089fd3de22'
down_revision = 'fe45a2e57d78'

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

def upgrade():
    op.add_column('calculation', sa.Column('metadata', postgresql.JSONB(astext_type=sa.Text()), nullable=False))
    op.add_column('test_result2', sa.Column('metadata', postgresql.JSONB(astext_type=sa.Text()), nullable=False))


def downgrade():
    op.drop_column('test_result2', 'metadata')
    op.drop_column('calculation', 'metadata')
