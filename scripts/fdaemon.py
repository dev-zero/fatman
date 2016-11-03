#!/usr/bin/env python

import socket
import logging
import subprocess
import os
from os import path
from time import sleep
from urllib.parse import urlparse

import click
import click_log
import requests

TASKS_URL = '{}/api/v2/tasks'

logger = logging.getLogger(__name__)


@click.command()
@click.option('--url', type=str, default='http://localhost:5000',
              show_default=True, help="The URL where FATMAN is running")
@click.option('--hostname', type=str, default=socket.gethostname,
              help="Override hostname-detection")
@click.option('--nap-time', type=str, default=5*60,
              show_default=True,
              help="Time to sleep if no new tasks are available")
@click.option('--data-dir', type=click.Path(exists=True, resolve_path=True),
              default='./fdaemon-data', show_default=True,
              help="Data directory")
@click_log.simple_verbosity_option()
@click_log.init(__name__)
def main(url, hostname, nap_time, data_dir):
    """FATMAN Calculation Runner Daemon"""

    os.chdir(data_dir)

    sess = requests.Session()
    sess.verify = False

    parsed_uri = urlparse(url)
    server = '{uri.scheme}://{uri.netloc}'.format(uri=parsed_uri)

    while True:
        logger.info("fetching new task")
        req = sess.get(TASKS_URL.format(url),
                       params={'limit': 1, 'status': 'new'})
        req.raise_for_status()
        tasks = req.json()

        if len(tasks) == 0:
            logger.info("no new tasks available, taking a nap")
            sleep(nap_time)
            continue

        req = sess.patch(server + tasks[0]['_links']['self'],
                         json={'status': 'pending', 'machine': hostname})
        req.raise_for_status()
        task = req.json()

        logger.info("aquired new task %s", task['id'])

        task_dir = path.join(data_dir, task['id'])
        os.mkdir(task_dir)

        # download each input file by streaming
        for infile in task['infiles']:
            req = sess.get(server + infile['_links']['download'], stream=True)
            req.raise_for_status()
            with open(path.join(task_dir, infile['name']), 'wb') as fhandle:
                for chunk in req.iter_content(1024):
                    fhandle.write(chunk)

        runner = task['settings']['machine']['runner']
        if runner == 'slurm':
            subprocess.check_call(['sbatch', task['settings']['cmd']],
                                  cwd=task_dir)
            req = sess.patch(server + task['_links']['self'],
                             json={'status': 'running'})
            req.raise_for_status()
            task = req.json()

        elif runner == 'download-only':
            req = sess.patch(server + task['_links']['self'],
                             json={'status': 'running'})
            req.raise_for_status()
            task = req.json()

        else:
            raise NotImplementedError("runner {} not yet implemented".format(
                runner))

if __name__ == '__main__':
    main()
