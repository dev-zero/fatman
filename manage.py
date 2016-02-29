#!/usr/bin/env python

from flask.ext.script import Manager
from flask.ext.script.commands import ShowUrls, Clean

from fatman import app, db, models

manager = Manager(app)

@manager.command
def initdb():
    """Initialize the database (structure and initial data)"""

    db.create_tables([models.Role,
                      models.User,
                      models.UserRole,
                      models.Structure,
                      models.Method,
                      models.Test,
                      models.TestStructure,
                      models.TaskStatus,
                      models.Task,
                      models.Result,
                      models.BasissetFamily,
                      models.PseudopotentialFamily,
                      models.BasisSet,
                      models.Pseudopotential,
                      models.TestResult,
                     ], safe=True)

    for name in ['new', 'pending', 'running', 'done', 'error', 'resting']:
        models.TaskStatus.create_or_get(name=name)

    for name in ['SZV-GTH',
                 'DZV-GTH',
                 'DZVP-GTH',
                 'TZVP-GTH',
                 'TZV2P-GTH',
                 'QZV2P-GTH',
                 'QZV3P-GTH',
                 'aug-DZVP-GTH',
                 'aug-TZVP-GTH',
                 'aug-TZV2P-GTH',
                 'aug-QZV2P-GTH',
                 'aug-QZV3P-GTH',
                 '6-31G*',
                 '6-311ppG3f2d',
                 '6-31ppG3f2d',
                 'TZVP-pob',
                 'DZVP-MOLOPT-SR-GTH',
                 'DZVP-MOLOPT-GTH',
                ]:
        models.BasissetFamily.create_or_get(name=name)

    for name in ['GTH-PBE', 'GTH-NLCC-PBE', 'GTH-NLCC2015-PBE', 'ALL']:
        models.PseudopotentialFamily.create_or_get(name=name)


@manager.command
def cleardb():
    """Delete all the data in the Structure, Test, TestStructure, and Task databases.
       The Result and Method tables remain intact."""

    for table in [models.Structure,
                  models.Test,
                  models.TestStructure,
                  models.Task,
                  models.BasisSet,
                  models.BasissetFamily,
                  models.Pseudopotential,
                  models.PseudopotentialFamily,
                 ]:
        table.delete().execute()

@manager.command
def createconfig():
    """Create initial configuration file"""

    from os import urandom

    with open('fatman.cfg', 'w') as cfg:
        cfg.write("SECRET_KEY = {}\n".format(urandom(24)))

# add predefined commands from flask.script
manager.add_command("show_urls", ShowUrls())
manager.add_command("clean", Clean())

@manager.shell
def make_shell_context():
    """Automatically load our app, db and models in the shell started by `manage.py shell`"""
    return dict(app=app, db=db, models=models)

if __name__ == '__main__':
    manager.run()
