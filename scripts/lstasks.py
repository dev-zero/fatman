#!/usr/bin/env python

from __future__ import print_function

import requests
import click
from tabulate import tabulate

TASKS_URL = '{}/tasks'
TASKS_OPTS = ['status', 'limit', 'method', 'pseudofamily', 'structure']

@click.command()
@click.option('--status', '-s', type=str, help='Only  list tasks with given status (default: all)')
@click.option('--limit', '-l', type=int, default=100, help='Limit the number of tasks to the given value (default: 100)')
@click.option('--machine', '-m', type=str, help='Machine used for a task')
@click.option('--pseudofamily', '-p', type=str, help='Pseudopotential family used in the method for a task')
@click.option('--structure', '-t', type=str, help='Structure calculated')
@click.option('--url', type=str, default='https://tctdb.chem.uzh.ch/fatman',
              help='The URL where FATMAN is running (default: https://tctdb.chem.uzh.ch/fatman)')
def lstasks(url, **kwargs):
    """Change state of a task"""

    sess = requests.Session()
    # the certificate is signed by the inofficial TC-Chem CA
    sess.verify = False

    params = {o: kwargs[o] for o in kwargs.keys() if o in TASKS_OPTS and kwargs[o]}
    params['timeorder'] = True

    req = sess.get(TASKS_URL.format(url), data=params)
    req.raise_for_status()

    tasks = req.json()

    headers = ['id', 'status', 'machine', 'mtime', 'pseudofamily', 'structure']
    table = [[t['id'], t['status'], t['machine'], t['mtime'],
              t['method']['pseudopotential'], t['structure']['name']]
             for t in tasks]

    print(tabulate(table, headers=headers))

if __name__ == "__main__":
    lstasks()
