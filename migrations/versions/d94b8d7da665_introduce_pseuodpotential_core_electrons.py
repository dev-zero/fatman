"""introduce pseuodpotential.core_electrons

Revision ID: d94b8d7da665
Revises: 61a24179b161
Create Date: 2016-10-04 13:22:20.108499

"""

# revision identifiers, used by Alembic.
revision = 'd94b8d7da665'
down_revision = '61a24179b161'

from alembic import op
import sqlalchemy as sa

from sqlalchemy.sql import table, column
from sqlalchemy.schema import ForeignKey
from sqlalchemy.dialects.postgresql import UUID
from sqlalchemy.orm import sessionmaker

Session = sessionmaker()

def upgrade():
    op.alter_column('basis_set', 'family_id',
               existing_type=sa.INTEGER(),
               nullable=True)
    op.add_column('pseudopotential',
        sa.Column('core_electrons', sa.Integer(), nullable=True))
    op.add_column('pseudopotential',
        sa.Column('description', sa.Text(), nullable=True))
    op.alter_column('pseudopotential', 'family_id',
               existing_type=sa.INTEGER(),
               nullable=True)
    op.create_unique_constraint(None, 'pseudopotential', ['element', 'family_id', 'core_electrons', 'format'])

    pseudopotential_table = table(
        'pseudopotential',
        column('id', UUID(as_uuid=True)),
        column('element', sa.String(255)),
        column('pseudo', sa.Text),
        column('format', sa.String(255)),
        column('converted_from_id', UUID(as_uuid=True),
               ForeignKey('pseudopotential.id'))
    )

    sess = Session(bind=op.get_bind())

    cp2k_pseudos = (sess.query(pseudopotential_table)
                    .filter(pseudopotential_table.c.format == 'CP2K'))

    for pseudo in cp2k_pseudos:
        op.execute("UPDATE pseudopotential SET core_electrons = {q} WHERE id = '{id}'".format(
                   q=sum(map(int, pseudo.pseudo.splitlines()[0].split())), id=pseudo.id))

    op.execute(("UPDATE pseudopotential SET core_electrons = subquery.core_electrons"
                " FROM (SELECT t.id AS tid, s.core_electrons AS core_electrons"
                "       FROM pseudopotential AS t JOIN pseudopotential AS s"
                "       ON s.id = t.converted_from_id) AS subquery"
                " WHERE id = subquery.tid"))

    op.alter_column('pseudopotential', 'core_electrons',
               existing_type=sa.INTEGER(),
               nullable=False)

def downgrade():
    op.drop_constraint(None, 'pseudopotential', type_='unique')
    op.alter_column('pseudopotential', 'family_id',
               existing_type=sa.INTEGER(),
               nullable=False)
    op.drop_column('pseudopotential', 'description')
    op.drop_column('pseudopotential', 'core_electrons')
    op.alter_column('basis_set', 'family_id',
               existing_type=sa.INTEGER(),
               nullable=False)
