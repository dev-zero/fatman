#!/usr/bin/env python

from flask.ext.script import Manager
from flask.ext.script.commands import ShowUrls, Clean
from flask.ext.security.script import \
        AddRoleCommand, RemoveRoleCommand, \
        CreateUserCommand, ActivateUserCommand, DeactivateUserCommand

from fatman import app, db, models, user_datastore

manager = Manager(app)

@manager.command
def initdb():
    """Initialize the database (structure and initial data)
       A VIEW  is also created by manually executing a SQL statement,
       as there doesn't seem to be any other way to create one.
    """

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
    
    db.execute_sql("""
     DROP VIEW resultwithouttestresult;
     CREATE VIEW resultwithouttestresult AS 
     SELECT result.* FROM result JOIN task a  ON a.id = result.task_id 
         WHERE NOT EXISTS 
             (SELECT 1 FROM testresult b 
                  WHERE a.method_id = b.method_id AND 
                        a.test_id = b.test_id
             );
    """)

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
                 'TZVP-MOLOPT-GTH',
                 'TZV2PX-MOLOPT-GTH',
                 'pc-1',
                 'pc-2',
                 'pc-3',
                 'pc-4',
                ]:
        models.BasissetFamily.create_or_get(name=name)

    for name in ['GTH-PBE', 'GTH-NLCC-PBE', 'GTH-NLCC2015-PBE', 'ALL']:
        models.PseudopotentialFamily.create_or_get(name=name)

    user_datastore.find_or_create_role(name='admin', description='Administrator')


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

# add commands from flask.security
manager.add_command("create_user", CreateUserCommand())
manager.add_command("add_role", AddRoleCommand())
manager.add_command("remove_role", RemoveRoleCommand())
manager.add_command("deactivate_user", DeactivateUserCommand())
manager.add_command("activate_user", ActivateUserCommand())

# add predefined commands from flask.script
manager.add_command("show_urls", ShowUrls())
manager.add_command("clean", Clean())

@manager.shell
def make_shell_context():
    """Automatically load our app, db and models in the shell started by `manage.py shell`"""
    return dict(app=app, db=db, models=models)

if __name__ == '__main__':
    manager.run()
