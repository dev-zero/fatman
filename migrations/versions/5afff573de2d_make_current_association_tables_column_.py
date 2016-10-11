"""make current association tables column-combi a PK

Revision ID: 5afff573de2d
Revises: d94b8d7da665
Create Date: 2016-10-11 12:48:57.011473

"""

# revision identifiers, used by Alembic.
revision = '5afff573de2d'
down_revision = 'd94b8d7da665'

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql

def upgrade():
    op.create_primary_key(
        "user_role_pkey", "user_role",
        ["user_id", "role_id"]
        )

    op.create_primary_key(
        "structure_set_structure_pkey", "structure_set_structure",
        ["structure_id", "set_id"]
        )

    op.create_primary_key(
        "test_structure_pkey", "test_structure",
        ["test_id", "structure_id"]
        )


def downgrade():
    op.drop_constraint("user_role_pkey", "user_role")
    op.drop_constriant("structure_set_structure_pkey",
                       "structure_set_structure")
    op.drop_constriant("test_structure_pkey", "test_structure")
