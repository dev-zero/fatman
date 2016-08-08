#!/usr/bin/env python

from __future__ import print_function

import requests
import click

TASKS_URL = '{}/tasks'

@click.command()
@click.argument('state', type=str)
@click.argument('tids', type=int, nargs=-1, metavar='[task ids ...]')
@click.option('--url', type=str, default='http://localhost:5000',
              help='The URL where FATMAN is running (default: http://localhost:5000/)')
def chtask(state, tids, url):
    """Change state of a task"""

    sess = requests.Session()
    # the certificate is signed by the inofficial TC-Chem CA
    sess.verify = False

    for tid in tids:
        req = sess.patch(TASKS_URL.format(url) + "/{}".format(tid), data={'status': state})
        req.raise_for_status()

if __name__ == "__main__":
    chtask()
