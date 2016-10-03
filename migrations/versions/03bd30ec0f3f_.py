"""empty message

Revision ID: 03bd30ec0f3f
Revises: None
Create Date: 2016-10-03 13:26:59.338489

"""

# revision identifiers, used by Alembic.
revision = '03bd30ec0f3f'
down_revision = None

from alembic import op
import sqlalchemy as sa
from sqlalchemy.dialects import postgresql
from sqlalchemy.sql import table, column

def upgrade():
    op.execute('CREATE EXTENSION IF NOT EXISTS pgcrypto')

    op.create_table('basis_set_family',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('name', sa.String(length=255), nullable=False),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('name')
    )
    op.create_table('pseudopotential_family',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('name', sa.String(length=255), nullable=False),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('name')
    )
    op.create_table('role',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('name', sa.String(length=255), nullable=False),
        sa.Column('description', sa.Text(), nullable=True),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('name')
    )
    op.create_table('structure',
        sa.Column('id', postgresql.UUID(as_uuid=True), server_default=sa.text('gen_random_uuid()'), nullable=False),
        sa.Column('name', sa.String(), nullable=False),
        sa.Column('ase_structure', sa.Text(), nullable=False),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('name')
    )
    op.create_table('structure_set',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('name', sa.String(length=255), nullable=False),
        sa.Column('description', sa.Text(), nullable=True),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('name')
    )
    op.create_table('task_status',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('name', sa.String(length=255), nullable=False),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('name')
    )
    op.create_table('test',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('name', sa.String(length=255), nullable=False),
        sa.Column('description', sa.Text(), nullable=True),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('name')
    )
    op.create_table('user',
        sa.Column('id', sa.Integer(), nullable=False),
        sa.Column('email', sa.String(length=255), nullable=False),
        sa.Column('password', sa.Text(), nullable=False),
        sa.Column('active', sa.Boolean(), nullable=False),
        sa.Column('confirmed_at', sa.DateTime(), nullable=True),
        sa.PrimaryKeyConstraint('id'),
        sa.UniqueConstraint('email')
    )
    op.create_table('basis_set',
        sa.Column('id', postgresql.UUID(as_uuid=True), server_default=sa.text('gen_random_uuid()'), nullable=False),
        sa.Column('element', sa.String(length=255), nullable=False),
        sa.Column('family_id', sa.Integer(), nullable=False),
        sa.Column('basis', sa.Text(), nullable=False),
        sa.ForeignKeyConstraint(['family_id'], ['basis_set_family.id'], ),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_table('method',
        sa.Column('id', postgresql.UUID(as_uuid=True), server_default=sa.text('gen_random_uuid()'), nullable=False),
        sa.Column('code', sa.String(length=255), nullable=False),
        sa.Column('pseudopotential_id', sa.Integer(), nullable=False),
        sa.Column('basis_set_id', sa.Integer(), nullable=False),
        sa.Column('settings', postgresql.JSONB(), nullable=True),
        sa.ForeignKeyConstraint(['basis_set_id'], ['basis_set_family.id'], ),
        sa.ForeignKeyConstraint(['pseudopotential_id'], ['pseudopotential_family.id'], ),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_table('pseudopotential',
        sa.Column('id', postgresql.UUID(as_uuid=True), server_default=sa.text('gen_random_uuid()'), nullable=False),
        sa.Column('element', sa.String(length=255), nullable=False),
        sa.Column('pseudo', sa.Text(), nullable=False),
        sa.Column('family_id', sa.Integer(), nullable=False),
        sa.Column('format', sa.String(length=255), nullable=False),
        sa.Column('converted_from_id', postgresql.UUID(as_uuid=True), nullable=True),
        sa.ForeignKeyConstraint(['converted_from_id'], ['pseudopotential.id'], ),
        sa.ForeignKeyConstraint(['family_id'], ['pseudopotential_family.id'], ),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_table('structure_set_structure',
        sa.Column('structure_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('set_id', sa.Integer(), nullable=False),
        sa.ForeignKeyConstraint(['set_id'], ['structure_set.id'], ),
        sa.ForeignKeyConstraint(['structure_id'], ['structure.id'], )
    )
    op.create_table('test_structure',
        sa.Column('test_id', sa.Integer(), nullable=False),
        sa.Column('structure_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.ForeignKeyConstraint(['structure_id'], ['structure.id'], ),
        sa.ForeignKeyConstraint(['test_id'], ['test.id'], )
    )
    op.create_table('user_role',
        sa.Column('user_id', sa.Integer(), nullable=False),
        sa.Column('role_id', sa.Integer(), nullable=False),
        sa.ForeignKeyConstraint(['role_id'], ['role.id'], ),
        sa.ForeignKeyConstraint(['user_id'], ['user.id'], )
    )
    op.create_table('task',
        sa.Column('id', postgresql.UUID(as_uuid=True), server_default=sa.text('gen_random_uuid()'), nullable=False),
        sa.Column('structure_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('method_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('status_id', sa.Integer(), nullable=False),
        sa.Column('test_id', sa.Integer(), nullable=False),
        sa.Column('ctime', sa.DateTime(), nullable=False),
        sa.Column('mtime', sa.DateTime(), nullable=False),
        sa.Column('machine', sa.String(length=255), nullable=False),
        sa.Column('priority', sa.Integer(), nullable=True),
        sa.ForeignKeyConstraint(['method_id'], ['method.id'], ),
        sa.ForeignKeyConstraint(['status_id'], ['task_status.id'], ),
        sa.ForeignKeyConstraint(['structure_id'], ['structure.id'], ),
        sa.ForeignKeyConstraint(['test_id'], ['test.id'], ),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_table('test_result',
        sa.Column('id', postgresql.UUID(as_uuid=True), server_default=sa.text('gen_random_uuid()'), nullable=False),
        sa.Column('ctime', sa.DateTime(), nullable=False),
        sa.Column('test_id', sa.Integer(), nullable=False),
        sa.Column('method_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('result_data', postgresql.JSONB(), nullable=True),
        sa.ForeignKeyConstraint(['method_id'], ['method.id'], ),
        sa.ForeignKeyConstraint(['test_id'], ['test.id'], ),
        sa.PrimaryKeyConstraint('id')
    )
    op.create_table('result',
        sa.Column('id', postgresql.UUID(as_uuid=True), server_default=sa.text('gen_random_uuid()'), nullable=False),
        sa.Column('energy', sa.Float(), nullable=False),
        sa.Column('task_id', postgresql.UUID(as_uuid=True), nullable=False),
        sa.Column('filename', sa.String(length=255), nullable=True),
        sa.Column('data', postgresql.JSONB(), nullable=True),
        sa.ForeignKeyConstraint(['task_id'], ['task.id'], ),
        sa.PrimaryKeyConstraint('id')
    )

    task_status_table = table('task_status',
        column('id', sa.Integer),
        column('name', sa.String)
        )

    op.bulk_insert(task_status_table,
            [
                {'id': 1, 'name': "new"},
                {'id': 2, 'name': "pending"},
                {'id': 3, 'name': "running"},
                {'id': 4, 'name': "done"},
                {'id': 5, 'name': "error"},
                {'id': 7, 'name': "running-remote"},
                {'id': 6, 'name': "deferred"},
                {'id': 8, 'name': "cancelled"},
            ]
        )

    basis_set_family_table = table(
        'basis_set_family',
        column('id', sa.Integer),
        column('name', sa.String)
    )

    op.bulk_insert(
        basis_set_family_table,
        [
            {'id':   1, 'name': "SZV-GTH"},
            {'id':   2, 'name': "DZV-GTH"},
            {'id':   3, 'name': "DZVP-GTH"},
            {'id':   4, 'name': "TZVP-GTH"},
            {'id':   5, 'name': "TZV2P-GTH"},
            {'id':   6, 'name': "QZV2P-GTH"},
            {'id':   7, 'name': "QZV3P-GTH"},
            {'id':   8, 'name': "aug-DZVP-GTH"},
            {'id':   9, 'name': "aug-TZVP-GTH"},
            {'id':  10, 'name': "aug-TZV2P-GTH"},
            {'id':  11, 'name': "aug-QZV2P-GTH"},
            {'id':  12, 'name': "aug-QZV3P-GTH"},
            {'id':  13, 'name': "6-31G*"},
            {'id':  14, 'name': "6-311ppG3f2d"},
            {'id':  15, 'name': "6-31ppG3f2d"},
            {'id':  16, 'name': "TZVP-pob"},
            {'id':  17, 'name': "DZVP-MOLOPT-SR-GTH"},
            {'id':  18, 'name': "DZVP-MOLOPT-GTH"},
            {'id': 257, 'name': "TZV2PX-MOLOPT-GTH"},
            {'id': 258, 'name': "pc-1"},
            {'id': 259, 'name': "pc-2"},
            {'id': 260, 'name': "pc-3"},
            {'id': 261, 'name': "pc-4"},
            {'id': 202, 'name': "TZVP-MOLOPT-GTH"},
            {'id':  55, 'name': "planewave"},
            {'id': 164, 'name': "DZ-ANO"},
            {'id': 165, 'name': "DZVP-ALL"},
        ]
    )

    pseudopotential_family_table = table(
        'pseudopotential_family',
        column('id', sa.Integer),
        column('name', sa.String)
    )

    op.bulk_insert(
        pseudopotential_family_table,
        [
            {'id':  1, 'name': "GTH-PBE"},
            {'id':  2, 'name': "GTH-NLCC-PBE"},
            {'id':  3, 'name': "GTH-NLCC2015-PBE"},
            {'id':  4, 'name': "ALL"},
            {'id': 13, 'name': "all-electron"},
        ]
    )


def downgrade():
    op.drop_table('result')
    op.drop_table('test_result')
    op.drop_table('task')
    op.drop_table('user_role')
    op.drop_table('test_structure')
    op.drop_table('structure_set_structure')
    op.drop_table('pseudopotential')
    op.drop_table('method')
    op.drop_table('basis_set')
    op.drop_table('user')
    op.drop_table('test')
    op.drop_table('task_status')
    op.drop_table('structure_set')
    op.drop_table('structure')
    op.drop_table('role')
    op.drop_table('pseudopotential_family')
    op.drop_table('basis_set_family')
