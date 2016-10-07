#!/usr/bin/env python

import requests
import click
import collections
from tabulate import tabulate

TASKS_URL = '{}/api/v1/tasks'

# from http://stackoverflow.com/questions/6027558/flatten-nested-python-dictionaries-compressing-keys
def flatten(d, parent_key='', sep='/'):
    items = []
    for k, v in d.items():
        new_key = parent_key + sep + k if parent_key else k
        if isinstance(v, collections.MutableMapping):
            items.extend(flatten(v, new_key, sep=sep).items())
        else:
            items.append((new_key, v))
    return dict(items)

@click.command()
@click.option('--status', '-s', type=str, help='Only  list tasks with given status (default: all)')
@click.option('--limit', '-l', type=int, default=100, help='Limit the number of tasks to the given value (default: 100)')
@click.option('--machine', '-m', type=str, help='Machine used for a task')
@click.option('--pseudofamily', '-p', type=str, help='Pseudopotential family used in the method for a task')
@click.option('--structure', '-t', type=str, help='Structure calculated')
@click.option('--format', '-f', 'formatstring', type=str, help='Format string for each output line')
@click.option('--url', type=str, default='https://tctdb.chem.uzh.ch/fatman',
              help='The URL where FATMAN is running (default: https://tctdb.chem.uzh.ch/fatman)')
def lstasks(url, formatstring, **params):
    """Change state of a task"""

    sess = requests.Session()
    # the certificate is signed by the inofficial TC-Chem CA
    sess.verify = False

    # filter None values
    params = {k: params[k] for k in params.keys() if params[k]}
    params['timeorder'] = True

    req = sess.get(TASKS_URL.format(url), data=params)
    req.raise_for_status()

    tasks = req.json()

    if formatstring:
        for task in tasks:
            click.echo(formatstring.format(**flatten(task)))

    else:
        headers = ['id', 'status', 'machine', 'mtime', 'pseudofamily', 'structure']
        table = [[t['id'], t['status'], t['machine'], t['mtime'],
                  t['method']['pseudopotential'], t['structure']['name']]
                 for t in tasks]
        click.echo(tabulate(table, headers=headers))

if __name__ == "__main__":
    lstasks()
