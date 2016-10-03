"""move result.energy to result.data.total_energy

Revision ID: 61a24179b161
Revises: 03bd30ec0f3f
Create Date: 2016-10-03 23:03:58.293807

"""

# revision identifiers, used by Alembic.
revision = '61a24179b161'
down_revision = '03bd30ec0f3f'

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

def upgrade():
    op.execute("UPDATE result SET data = jsonb_set(data, '{total_energy}', to_jsonb(energy), true)")
    op.drop_column('result', 'energy')
