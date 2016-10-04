#!/usr/bin/env python

from flask_script import Manager
from flask_script.commands import ShowUrls, Clean
from flask_security.script import \
        AddRoleCommand, RemoveRoleCommand, \
        CreateUserCommand, ActivateUserCommand, DeactivateUserCommand
from flask_migrate import MigrateCommand

from fatman import app, db, models, tasks, migrate

manager = Manager(app)

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

manager.add_command("db", MigrateCommand)

@manager.shell
def make_shell_context():
    """Automatically load our app, db and models in the shell started by `manage.py shell`"""
    return dict(app=app, db=db, models=models, tasks=tasks)

if __name__ == '__main__':
    manager.debug=True
    manager.run()
