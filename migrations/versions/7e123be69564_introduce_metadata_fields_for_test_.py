"""introduce metadata fields for test result2 and calculation

Revision ID: 7e123be69564
Revises: fe45a2e57d78
Create Date: 2017-06-20 16:51:26.137576

"""

# revision identifiers, used by Alembic.
revision = '7e123be69564'
down_revision = 'fe45a2e57d78'

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

def upgrade():
    op.add_column('calculation', sa.Column('metadata', postgresql.JSONB(astext_type=sa.Text()), server_default=sa.text("'{}'::json"), nullable=False))
    op.add_column('test_result2', sa.Column('metadata', postgresql.JSONB(astext_type=sa.Text()), server_default=sa.text("'{}'::json"), nullable=False))


def downgrade():
    op.drop_column('test_result2', 'metadata')
    op.drop_column('calculation', 'metadata')
