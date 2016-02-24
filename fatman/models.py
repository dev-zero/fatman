
from peewee import *
from playhouse.postgres_ext import BinaryJSONField

from fatman import db

class BaseModel(Model):
    class Meta:
        database = db

class Structure(BaseModel):
    name = CharField(unique=True, null=False)
    ase_structure = TextField(null=False)
    class Meta:
        order_by = ('name',)

    def __str__(self):
        return self.name

class BasissetFamily(BaseModel):
    name = CharField(unique=True, null=False)

    class Meta:
        order_by = ('name',)

    def __str__(self):
        return self.name

class PseudopotentialFamily(BaseModel):
    name = CharField(unique=True, null=False)

    class Meta:
        order_by = ('name',)

    def __str__(self):
        return self.name


class BasisSet(BaseModel):
    family = ForeignKeyField(BasissetFamily, related_name='basisset')
    element = CharField(null=False)
    basis = TextField(null=False)

    class Meta:
        order_by = ('family',)

    def __str__(self):
        return self.family.name

class Pseudopotential(BaseModel):
    family = ForeignKeyField(PseudopotentialFamily, related_name='pseudopotential')
    element = CharField(null=False)
    pseudo = TextField(null=False)

    class Meta:
        order_by = ('family',)

    def __str__(self):
        return self.family.name

class Method(BaseModel):
    # proper database design would demand introduction of
    # separate entities for the following fields, but let's make it easy
    code = CharField()
    pseudopotential = ForeignKeyField(Pseudopotential, related_name = 'method_pp')
    basis_set = ForeignKeyField(BasisSet, related_name = 'method_bs')
    settings = BinaryJSONField(null=True)

    class Meta:
        order_by = ('code',)

    def __str__(self):
        return "code: {}, pseudopotential: {}, basis set: {}".format(
                self.code,
                self.pseudopotential,
                self.basis_set)

class Test(BaseModel):
    name = CharField(unique=True, null=False)
    description = TextField(null=True)

    class Meta:
        order_by = ('name',)

    def __str__(self):
        return self.name

class TestStructure(BaseModel):
    structure = ForeignKeyField(Structure, related_name='tests')
    test = ForeignKeyField(Test, related_name='structures')

class TaskStatus(BaseModel):
    name = CharField(unique=True, null=False)

    class Meta:
        order_by = ('name',)

    def __str__(self):
        return self.name

class Task(BaseModel):
    structure = ForeignKeyField(Structure, related_name='tasks')
    method = ForeignKeyField(Method, related_name='tasks')
    status = ForeignKeyField(TaskStatus)
    ctime = DateTimeField()
    mtime = DateTimeField()
    machine = CharField()

    def __str__(self):
        return "id: {}, structure: {}, status: {}".format(
                self.id,
                self.structure,
                self.status)

class Result(BaseModel):
    energy = DoubleField()
    task = ForeignKeyField(Task, related_name='results')
    filename = CharField()



