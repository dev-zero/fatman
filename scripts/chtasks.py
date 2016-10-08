#!/usr/bin/env python

import requests
import click
from tabulate import tabulate

TASKS_URL = '{}/api/v1/tasks'
TASK_URL = '{}/api/v1/tasks/{}'

@click.command()
@click.option('--status', '-s', type=str, help='only tasks with given status (default: all)')
@click.option('--machine', '-m', type=str, help='machine used for a task')
@click.option('--pseudofamily', '-p', type=str,
              help='pseudopotential family used in the method for a task')
@click.option('--structure', '-t', type=str, help='structure calculated (also partially matched)')
@click.option('--id', '-i', 'ids', type=click.UUID, multiple=True,
              help='task ID to change (can be specified multiple times)')
@click.option('--url', type=str, default='https://tctdb.chem.uzh.ch/fatman',
              help='The URL where FATMAN is running (default: https://tctdb.chem.uzh.ch/fatman)')
@click.argument('newstatus', type=str)
def chtasks(url, newstatus, ids, **params):
    """Change status of a task"""

    if not any(params.values()) and not ids:
        raise click.UsageError("at least one filter option must be used")

    sess = requests.Session()
    # the certificate is signed by the inofficial TC-Chem CA
    sess.verify = False

    tasks = []

    # we may have gotten only '-i' arguments
    if any(params.values()):
        params['timeorder'] = True

        req = sess.get(TASKS_URL.format(url), data=params)
        req.raise_for_status()

        tasks += req.json()

    if ids:
        for tid in ids:
            req = sess.get(TASK_URL.format(url, tid))
            req.raise_for_status()

            tasks += [req.json()]

    if not tasks:
        click.echo('No tasks found matching the criteria')
        return

    headers = ['id', 'status', 'machine', 'mtime', 'pseudofamily', 'structure']
    table = [[t['id'], t['status'], t['machine'], t['mtime'],
              t['method']['pseudopotential'], t['structure']['name']]
             for t in tasks]

    click.echo(tabulate(table, headers=headers))

    click.confirm("The status of all listed tasks will be set to '{}', continue?".format(newstatus),
                  abort=True)

    for task in tasks:
        req = sess.patch(TASKS_URL.format(url) + "/{}".format(task['id']),
                         data={'status': newstatus})
        req.raise_for_status()

if __name__ == "__main__":
    chtasks()
