
# Design rationale wrt foreign key relationship specifications:
# Whenever referencing by integer keys (internal object references)
# we specify the ORM loading directly in the model using lazy='joined'
# and innerjoin=True
# All other references shall be loaded when during the query and as-needed

import uuid
from datetime import datetime as dt
try:
    from urllib.parse import urlsplit
except ImportError:
    from urlparse import urlsplit
from os import path

from flask_security import UserMixin, RoleMixin

from sqlalchemy import text, or_
from sqlalchemy import Column, ForeignKey, UniqueConstraint, CheckConstraint, Index
from sqlalchemy import Integer, String, Boolean, DateTime, Text, Enum
from sqlalchemy.dialects.postgresql import UUID, JSONB, ExcludeConstraint
from sqlalchemy.orm import column_property
from sqlalchemy.sql.expression import null
from sqlalchemy.sql.functions import coalesce

from werkzeug.datastructures import FileStorage

from . import resultfiles

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
                        primary_key=True),
                 Column('role_id', Integer, ForeignKey('role.id'),
                        primary_key=True))


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
    name = Column(String, nullable=False)
    # ase_structure = Column(JSONB, nullable=False)
    ase_structure = Column(Text, nullable=False)

    replaced_by_id = Column(UUID(as_uuid=True), ForeignKey('structure.id'))
    replaced_by = relationship("Structure", remote_side=[id], lazy='joined', join_depth=2)

    # this is also required to make the cascade work when deleting a structure
    replaced = relationship("Structure", remote_side=[replaced_by_id], lazy='noload')

    def __repr__(self):
        return "<Structure(name='{}')>".format(self.name)

    def __str__(self):
        return self.name

    __table_args__ = (
        # ensure that the name is unique amongst non-replaced structures,
        # and defer constraint to end of transaction to be able to replace
        # structures within one transaction
        ExcludeConstraint(
            ('name', '='),
            using='btree',
            where=replaced_by_id == null(),
            deferrable=True, initially='DEFERRED'),
        )


StructureSetStructure = Table('structure_set_structure', Base.metadata,
                              Column('structure_id', UUID(as_uuid=True),
                                     ForeignKey('structure.id'),
                                     primary_key=True),
                              Column('set_id', Integer,
                                     ForeignKey('structure_set.id'),
                                     primary_key=True))


class StructureSet(Base):
    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True, nullable=False)
    description = Column(Text)
    structures = relationship('Structure', secondary="structure_set_structure",
                              backref=backref('sets'))

    superset_id = Column(Integer, ForeignKey('structure_set.id'))
    subsets = relationship('StructureSet',
                           backref=backref('superset', remote_side=[id]))

    def __repr__(self):
        return "<StructureSet(name='{}')>".format(self.name)

    def __str__(self):
        return self.name


class BasisSet(Base):
    id = UUIDPKColumn()
    element = Column(String(255), nullable=False)
    family_id = Column(Integer, ForeignKey("basis_set_family.id"))
    family = relationship("BasisSetFamily", lazy='joined',
                          backref=backref('basissets', lazy='dynamic'))
    basis = Column(Text, nullable=False)

    def __repr__(self):
        return "<BasisSet(id='{}', element='{}', family_id={})>".format(
                self.id, self.element, self.family_id)

    def __str__(self):
        return "{} ({})".format(self.element, self.id)


class BasisSetFamily(Base):
    id = Column(Integer, primary_key=True)
    name = Column(String(255), unique=True, nullable=False)

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
    family_id = Column(Integer, ForeignKey("pseudopotential_family.id"))
    family = relationship("PseudopotentialFamily", lazy='joined',
                          backref=backref('pseudos', lazy='dynamic'))

    format = Column(String(255), nullable=False)
    core_electrons = Column(Integer, nullable=False)
    description = Column(Text)

    converted_from_id = Column(UUID(as_uuid=True), ForeignKey('pseudopotential.id'))
    converted_from = relationship("Pseudopotential", remote_side=[id])

    __table_args__ = (
        UniqueConstraint('element', 'family_id', 'core_electrons', 'format'),
    )

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
                             primary_key=True),
                      Column('structure_id', UUID(as_uuid=True), ForeignKey('structure.id'),
                             primary_key=True))


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


class Machine(Base):
    id = UUIDPKColumn()
    shortname = Column(String(255), nullable=False, unique=True)
    name = Column(String(255), nullable=False, unique=True)
    settings = Column(JSONB)

    def __repr__(self):
        return "<Machine(id='{}', shortname='{}')>".format(self.id, self.name)

    def __str__(self):
        return "{} ({})".format(self.name, self.shortname)


class Code(Base):
    id = UUIDPKColumn()
    name = Column(String(255), nullable=False, unique=True)
    pseudo_format = Column(String(255), nullable=False)

    def __repr__(self):
        return "<Code(id='{}', name='{}')>".format(self.id, self.name)

    def __str__(self):
        return self.name


class Command(Base):
    """Per-machine and per-code commands to run codes.

    Since we want exactly one command per machine and code
    we use a combined PK (again) instead of introducing a new one.
    """
    code_id = Column(UUID(as_uuid=True),
                     ForeignKey('code.id'),
                     primary_key=True)
    code = relationship("Code", backref="commands")
    machine_id = Column(UUID(as_uuid=True),
                        ForeignKey('machine.id'),
                        primary_key=True)
    machine = relationship("Machine", backref="commands")
    environment = Column(JSONB)  # for setting env vars & loading modules
    commands = Column(JSONB, nullable=False)  # the actual commands to run


class Calculation(Base):
    id = UUIDPKColumn()

    collection_id = Column(UUID(as_uuid=True),
                           ForeignKey('calculation_collection.id'),
                           nullable=False)
    collection = relationship("CalculationCollection",
                              lazy='joined', innerjoin=True,
                              backref=backref("calculations", lazy='dynamic'))

    test_id = Column(Integer, ForeignKey('test.id'))
    test = relationship("Test", backref="calculations")

    structure_id = Column(UUID(as_uuid=True), ForeignKey('structure.id'),
                          nullable=False)
    structure = relationship("Structure", backref="calculations")

    code_id = Column(UUID(as_uuid=True),
                     ForeignKey('code.id'),
                     nullable=False)
    code = relationship("Code", backref="calculations")

    pseudos = relationship("Pseudopotential",
                           secondary="calculation_pseudopotential",
                           backref="calculations",
                           cascade="all", passive_deletes=True)
    basis_sets = relationship("BasisSet",
                              secondary="calculation_basis_set",
                              backref="calculations")
    testresults_query = relationship("TestResult2", secondary="test_result2_calculation", lazy='dynamic')
    tasks_query = relationship("Task2", lazy='dynamic', cascade="all", passive_deletes=True)

    settings = Column(JSONB)
    restrictions = Column(JSONB)
    results = Column(JSONB)

    results_available = column_property(results.isnot(None))

    @property
    def current_task(self):
        try:
            return self.tasks[0]
        except IndexError:
            return None

    # for later: link together already defined calculations (but full copy existing ones)
    # parent_id = Column(UUID(as_uuid=True), ForeignKey('calculation.id'))
    # parent = relationship("Calculation", remote_side=[id])

    def __repr__(self):
        return "<Calculation(id='{}')>".format(self.id)


CalculationPseudopotential = Table(
    'calculation_pseudopotential', Base.metadata,
    Column('calculation_id', UUID(as_uuid=True), ForeignKey('calculation.id', ondelete='CASCADE'),
           primary_key=True),
    Column('pseudo_id', UUID(as_uuid=True), ForeignKey('pseudopotential.id'),
           primary_key=True))


class CalculationBasisSet(Base):
    """Association table for Calculation <-> BasisSet

    See http://docs.sqlalchemy.org/en/rel_1_1/orm/basic_relationships.html#association-object
    but I left out the relationship definitions since I don't need them (?)
    """

    calculation_id = Column(UUID(as_uuid=True),
                            ForeignKey('calculation.id', ondelete='CASCADE'),
                            primary_key=True)
    basis_set_id = Column(UUID(as_uuid=True), ForeignKey('basis_set.id'),
                          primary_key=True)
    btype = Column(Enum("default", "aux_fit", "aux", "lri", "ri_aux",
                        name="basis_set_type"),
                   primary_key=True, default="default")

    calculation = relationship("Calculation",
                               backref=backref("basis_set_associations",
                                               cascade="all", passive_deletes=True))
    basis_set = relationship("BasisSet", backref="calculation_associations")

    def __repr__(self):
        return ("<CalculationBasisSet("
                "calculation='{}', basis_set='{}', type='{}'"
                ")>").format(self.calculation_id,
                             self.basis_set_id,
                             self.btype)


class CalculationDefaultSettings(Base):
    """Default settings for calculations.

    The settings are hard-linked to code and test and
    optionally linked to a structure or a structure set.

    Since the structure and the structure set are optional,
    we can not form a PK out of (code, test, structure, structure set).

    Since NULL is not comparable (NULL == NULL = False), we coalesce NULL values in
    the structure and structure set columns to a 0 UUID, a 0 respectively before the
    uniqueness constraint is applied. Even if improbable, we further ensure that the
    structure and structure set columns can not contain that 0 UUID/0 integer.
    """

    id = UUIDPKColumn()

    code_id = Column(UUID(as_uuid=True), ForeignKey('code.id'), nullable=False)
    code = relationship("Code", backref="default_settings")

    test_id = Column(Integer, ForeignKey('test.id'), nullable=False)
    test = relationship("Test", backref="default_settings")

    structure_id = Column(UUID(as_uuid=True), ForeignKey('structure.id'))
    structure = relationship("Structure", backref="default_settings")

    structure_set_id = Column(Integer, ForeignKey('structure_set.id'))
    structure_set = relationship("StructureSet", backref="default_settings")

    settings = Column(JSONB)

    __table_args__ = (
        Index('unique_calculation_default_settings_idx',
              code_id,
              test_id,
              coalesce(structure_id, str(uuid.UUID(int=0))),
              coalesce(structure_set_id, 0),
              unique=True),
        CheckConstraint(structure_id != str(uuid.UUID(int=0))),
        CheckConstraint(structure_set_id != 0),
        # ensure that a setting is not restricted to both a structure and a structure set:
        CheckConstraint(or_(structure_set_id == None, structure_id == None)),
    )

    def __repr__(self):
        return ("<CalculationDefaultSettings("
                "code='{}', test='{}'"
                ")>").format(self.code_id, self.test_id)


class CalculationCollection(Base):
    id = UUIDPKColumn()
    name = Column(String(255), nullable=False, unique=True)
    desc = Column(Text)

    def __repr__(self):
        return "<CalculationCollection(id='{}', name='{}')>".format(self.id,
                                                                    self.name)

    def __str__(self):
        return self.name


class Task2(Base):
    id = UUIDPKColumn()
    calculation_id = Column(UUID(as_uuid=True), ForeignKey('calculation.id', ondelete='CASCADE'),
                            nullable=False)
    calculation = relationship("Calculation",
                               backref=backref("tasks", order_by=lambda: Task2.ctime.desc(),
                                               cascade="all", passive_deletes=True))
    status_id = Column(Integer, ForeignKey('task_status.id'), nullable=False)
    status = relationship("TaskStatus", lazy='joined', innerjoin=True)

    ctime = Column(DateTime, nullable=False, default=dt.now)
    mtime = Column(DateTime, nullable=False, default=dt.now, onupdate=dt.now)
    machine_id = Column(UUID(as_uuid=True), ForeignKey('machine.id'))
    machine = relationship("Machine", backref="tasks")
    priority = Column(Integer, default=0)
    data = Column(JSONB)  # task-related data like runtime, #(MPI nodes), etc.
    restrictions = Column(JSONB)
    settings = Column(JSONB)
    infiles = relationship("Artifact", secondary="task2_artifact",
                           lazy='dynamic', cascade="all", passive_deletes=True,
                           primaryjoin=("(Task2.id==Task2Artifact.task_id) &"
                                        "(Task2Artifact.linktype=='input')"))
    outfiles = relationship("Artifact", secondary="task2_artifact",
                            lazy='dynamic', cascade="all", passive_deletes=True,
                            primaryjoin=("(Task2.id==Task2Artifact.task_id) & "
                                         "(Task2Artifact.linktype=='output')"))

    def __repr__(self):
        return "<Task2(id='{}', status='{}')>".format(self.id, self.status)

    def __init__(self, calculation_id, priority=0):
        self.calculation_id = calculation_id
        self.status = TaskStatus.query.filter_by(name='new').one()


class Task2Artifact(Base):
    """Association table for Task2 <-> Artifact

    See http://docs.sqlalchemy.org/en/rel_1_1/orm/basic_relationships.html#association-object
    but I left out the task and artifact relationship definitions since I don't need them (?)
    """

    task_id = Column(UUID(as_uuid=True), ForeignKey('task2.id', ondelete='CASCADE'),
                     primary_key=True)
    artifact_id = Column(UUID(as_uuid=True), ForeignKey('artifact.id'),
                         primary_key=True)
    linktype = Column(Enum("input", "output", name="artifact_link_type"),
                      nullable=False)

    artifact = relationship("Artifact")
    task = relationship("Task2")


class Artifact(Base):
    id = UUIDPKColumn()
    name = Column(String(255), nullable=False)
    path = Column(String(255), nullable=False)
    mdata = Column('metadata', JSONB, default={}, nullable=False)
    # ^^ "metadata" is an SQLAlchemy reserved word

    def __repr__(self):
        return "<Artifact(id='{}', name='{}')>".format(self.id, self.name)

    def __init__(self, name, path, metadata={'compressed': None}):
        self.id = uuid.uuid4()
        self.name = name
        # replace {id} and {name} occurrances
        self.path = path.format(id=self.id, name=self.name)
        self.mdata = metadata

    def save(self, buf):
        """Save a byte-stream to the storage specified in path"""

        scheme, nwloc, fullpath, _, _ = urlsplit(self.path)

        if scheme == 'fkup' and nwloc == 'results':
            folder, filename = path.split(fullpath)
            storage = FileStorage(buf, filename=filename)
            resultfiles.save(storage, folder=folder[1:])
        else:
            raise RuntimeError("unknown scheme '{}' or location '{}'".format(
                scheme, nwloc))


class TestResult2(Base):
    id = UUIDPKColumn()
    test_id = Column(Integer, ForeignKey('test.id'), nullable=False)
    test = relationship("Test", lazy='joined')
    calculations = relationship("Calculation", secondary="test_result2_calculation",
                                backref='testresults')
    data = Column(JSONB)
    collections = relationship('TestResult2Collection',
                               secondary="test_result2_test_result2_collection",
                               backref='testresults', lazy='joined')

    def __repr__(self):
        return "<TestResult(id='{}')>".format(self.id)


TestResult2Calculation = Table(
    'test_result2_calculation',
    Base.metadata,
    Column('test_id', UUID(as_uuid=True), ForeignKey('test_result2.id'),
           primary_key=True),
    Column('calculation_id', UUID(as_uuid=True),
           ForeignKey('calculation.id'), primary_key=True))


TestResult2TestResult2Collection = Table(
    'test_result2_test_result2_collection',
    Base.metadata,
    Column('test_id', UUID(as_uuid=True), ForeignKey('test_result2.id'),
           primary_key=True),
    Column('collection_id', UUID(as_uuid=True),
           ForeignKey('test_result2_collection.id'), primary_key=True))


class TestResult2Collection(Base):
    id = UUIDPKColumn()
    name = Column(String(255), nullable=False, unique=True)
    desc = Column(Text)

    def __repr__(self):
        return "<TestResultCollection(id='{}', name='{}')>".format(
            self.id, self.name)

    def __str__(self):
        return self.name


class TaskRuntimeSettings(Base):
    """
    Settings which are applied as soon as we are going to run a task
    and at which point the machine and the command is finally known.

    Currently hard-linked to Command and Test.

    These settings can currently target the input, the runner arguments
    or the command environment based on Machine, Code and Test.
    """

    id = UUIDPKColumn()
    settings = Column(JSONB, nullable=False)

    machine_id = Column(UUID(as_uuid=True),
                        ForeignKey('machine.id'),
                        nullable=False)
    machine = relationship("Machine", backref="task_runtime_settings")
    code_id = Column(UUID(as_uuid=True),
                     ForeignKey('code.id'),
                     nullable=False)
    code = relationship("Code", backref="task_runtime_settings")
    test_id = Column(Integer,
                     ForeignKey('test.id'),
                     nullable=False)
    test = relationship("Test", backref="task_runtime_settings")
