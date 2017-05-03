"""db-level CASCADE on testresult removal

Revision ID: 66322a6bc199
Revises: 825e4b6229be
Create Date: 2017-05-03 11:15:55.575692

"""

# revision identifiers, used by Alembic.
revision = '66322a6bc199'
down_revision = '825e4b6229be'

from alembic import op
import sqlalchemy as sa

def upgrade():
    op.drop_constraint('test_result2_calculation_test_id_fkey', 'test_result2_calculation', type_='foreignkey')
    op.create_foreign_key(None, 'test_result2_calculation', 'test_result2', ['test_id'], ['id'], ondelete='CASCADE')


def downgrade():
    op.drop_constraint('test_result2_calculation_test_id_fkey', 'test_result2_calculation', type_='foreignkey')
    op.create_foreign_key('test_result2_calculation_test_id_fkey', 'test_result2_calculation', 'test_result2', ['test_id'], ['id'])
