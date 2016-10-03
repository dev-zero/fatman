
# Design rationale wrt foreign key relationship specifications:
# Whenever referencing by integer keys (internal object references)
# we specify the ORM loading directly in the model using lazy='joined'
# and innerjoin=True
# All other references shall be loaded when during the query and as-needed
from datetime import datetime as dt

from flask_security import UserMixin, RoleMixin

from sqlalchemy import text
from sqlalchemy import Column, ForeignKey
from sqlalchemy import Integer, Float, String, Boolean, DateTime, Text
from sqlalchemy.dialects.postgresql import UUID, JSONB

# Flask-SQLAlchemy wraps only a part of all the attributes from SQLAlchemy.ORM
# and is missing the PostgreSQL-specific types. To make the model definitions
# as # similar as possible to pure SQLAlchemy, we directly use SQLAlchemy types
# where possible and reuse the wrapped types from Flask-SQLAlchemy
# where necessary.
from . import db
Base = db.Model
Table = db.Table
relationship = db.relationship
relation = db.relation
backref = db.backref
dynamic_loader = db.dynamic_loader

# some shorthands


def UUIDPKColumn():
    """Constructs a Primary Key UUID Column"""
    return Column(UUID(as_uuid=True), server_default=text("gen_random_uuid()"),
                  primary_key=True)


# model definitions


class Role(Base, RoleMixin):
    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True, nullable=False)
    description = Column(Text)

    def __repr__(self):
        return "<Role(name='{}')>".format(self.name)

    def __str__(self):
        return self.name


UserRole = Table('user_role', Base.metadata,
                 Column('user_id', Integer, ForeignKey('user.id'),
                        nullable=False),
                 Column('role_id', Integer, ForeignKey('role.id'),
                        nullable=False))


class User(Base, UserMixin):
    id = Column(Integer, primary_key=True)
    email = Column(String(255), unique=True, nullable=False)
    password = Column(Text, nullable=False)
    active = Column(Boolean, default=False, nullable=False)
    confirmed_at = Column(DateTime)
    roles = relationship('Role', secondary="user_role",
                         backref=backref('users', lazy='dynamic'))

    def __repr__(self):
        return "<User(email='{}', active={})>".format(
                self.email, self.active)

    def __str__(self):
        return self.email


class Structure(Base):
    id = UUIDPKColumn()
    name = Column(String, unique=True, nullable=False)
    # ase_structure = Column(JSONB, nullable=False)
    ase_structure = Column(Text, nullable=False)

    def __repr__(self):
        return "<Structure(name='{}')>".format(self.name)

    def __str__(self):
        return self.name


StructureSetStructure = Table('structure_set_structure', Base.metadata,
                              Column('structure_id', UUID(as_uuid=True),
                                     ForeignKey('structure.id'),
                                     nullable=False),
                              Column('set_id', Integer,
                                     ForeignKey('structure_set.id'),
                                     nullable=False))


class StructureSet(Base):
    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True, nullable=False)
    description = Column(Text)
    structures = relationship('Structure', secondary="structure_set_structure",
                              backref=backref('sets', lazy='dynamic'))

    def __repr__(self):
        return "<StructureSet(name='{}')>".format(self.name)

    def __str__(self):
        return self.name


class BasisSet(Base):
    id = UUIDPKColumn()
    element = Column(String(255), nullable=False)
    family_id = Column(Integer, ForeignKey("basis_set_family.id"),
                       nullable=False)
    family = relationship("BasisSetFamily",
                          foreign_keys="BasisSet.family_id",
                          lazy='joined')
    basis = Column(Text, nullable=False)

    def __repr__(self):
        return "<BasisSet(id='{}', element='{}', family_id={})>".format(
                self.id, self.element, self.family_id)

    def __str__(self):
        return "{} ({})".format(self.element, self.id)


class BasisSetFamily(Base):
    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True, nullable=False)
    basissets = relationship("BasisSet", backref="families")

    __mapper_args__ = {
        "order_by": name,
    }

    def __repr__(self):
        return "<BasisSetFamily(name='{}')>".format(self.name)

    def __str__(self):
        return self.name


class Pseudopotential(Base):
    id = UUIDPKColumn()
    element = Column(String(255), nullable=False)
    pseudo = Column(Text, nullable=False)
    family_id = Column(Integer, ForeignKey("pseudopotential_family.id"),
                       nullable=False)
    family = relationship("PseudopotentialFamily",
                          lazy='joined', innerjoin=True,
                          backref=backref('pseudos', lazy='dynamic'))

    format = Column(String(255), nullable=False)
    converted_from_id = Column(UUID(as_uuid=True), ForeignKey('pseudopotential.id'))
    converted_from = relationship("Pseudopotential", remote_side=[id])

    def __repr__(self):
        return "<Pseudopotential(id='{}', element='{}', family_id={}, format={})>".format(
                self.id, self.element, self.family_id, self.format)

    def __str__(self):
        return "{} ({})".format(self.element, self.id)


class PseudopotentialFamily(Base):
    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True, nullable=False)

    def __repr__(self):
        return "<PseudopotentialFamily(name='{}')>".format(self.name)

    def __str__(self):
        return self.name


class Method(Base):
    # proper database design would demand introduction of
    # separate entity for the following field, but let's make it easy
    id = UUIDPKColumn()
    code = Column(String(255), nullable=False)
    pseudopotential_id = Column(Integer,
                                ForeignKey('pseudopotential_family.id'),
                                nullable=False)
    pseudopotential = relationship("PseudopotentialFamily",
                                   lazy='joined', innerjoin=True)
    basis_set_id = Column(Integer, ForeignKey('basis_set_family.id'),
                          nullable=False)
    basis_set = relationship("BasisSetFamily", lazy='joined', innerjoin=True)
    settings = Column(JSONB)

    __mapper_args__ = {
        "order_by": code,
    }

    def __repr__(self):
        return "<Method(id='{}')>".format(self.id)

    def __str__(self):
        settings_info = ""

        if self.settings:
            settingskeys = self.settings.keys()

            if "cutoff_rho" in settingskeys:
                settings_info = ", cutoff: {:.1f}".format(self.settings['cutoff_rho']/13.605692)

            if "rel_settings" in settingskeys:
                settings_info += ", relativistic"

            if "qs_settings" in settingskeys and "epsiso" in self.settings['qs_settings'].keys():
                settings_info += ", epsiso={}".format(self.settings['qs_settings']['epsiso'])

        return "ID: {}, code: {}, pseudopotential: {}, basis set: {} {}".format(
                self.id,
                self.code,
                self.pseudopotential,
                self.basis_set,
                settings_info
                )


TestStructure = Table('test_structure', Base.metadata,
                      Column('test_id', Integer, ForeignKey('test.id'),
                             nullable=False),
                      Column('structure_id', UUID(as_uuid=True), ForeignKey('structure.id'),
                             nullable=False))


class Test(Base):
    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True, nullable=False)
    description = Column(Text)
    structures = relationship('Structure', secondary='test_structure',
                              backref=backref('tests', lazy='dynamic'))

    __mapper_args__ = {
        "order_by": name,
    }

    def __repr__(self):
        return "<Test(name={})>".format(self.name)

    def __str__(self):
        return self.name


class TaskStatus(Base):
    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True, nullable=False)

    __mapper_args__ = {
        "order_by": name,
    }

    def __repr__(self):
        return "<TaskStatus(name={})>".format(self.name)

    def __str__(self):
        return self.name


class Task(Base):
    id = UUIDPKColumn()
    structure_id = Column(UUID(as_uuid=True), ForeignKey('structure.id'), nullable=False)
    structure = relationship("Structure")
    method_id = Column(UUID(as_uuid=True), ForeignKey('method.id'), nullable=False)
    method = relationship("Method")
    status_id = Column(Integer, ForeignKey('task_status.id'), nullable=False)
    status = relationship("TaskStatus", lazy='joined', innerjoin=True,
                          backref=backref('tasks', lazy='dynamic'))
    test_id = Column(Integer, ForeignKey('test.id'), nullable=False)
    test = relationship("Test")
    ctime = Column(DateTime, nullable=False, default=dt.now)
    mtime = Column(DateTime, nullable=False, default=dt.now, onupdate=dt.now)
    machine = Column(String(255), nullable=False)
    priority = Column(Integer, default=0)

    def __repr__(self):
        return "<Task(id='{}', status='{}')>".format(self.id, self.status)


class Result(Base):
    id = UUIDPKColumn()
    task_id = Column(UUID(as_uuid=True), ForeignKey('task.id'), nullable=False)
    task = relationship("Task", backref=backref('results', lazy='dynamic'))
    filename = Column(String(255))
    data = Column(JSONB)


class TestResult(Base):
    id = UUIDPKColumn()
    ctime = Column(DateTime, nullable=False, default=dt.now)
    test_id = Column(Integer, ForeignKey('test.id'), nullable=False)
    test = relationship("Test", backref=backref("results", lazy='dynamic'))
    method_id = Column(UUID(as_uuid=True), ForeignKey('method.id'),
                       nullable=False)
    method = relationship("Method", backref=backref("results", lazy='dynamic'))
    result_data = Column(JSONB)
