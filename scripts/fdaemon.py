#!/usr/bin/env python

import socket
import logging
import subprocess
import os
from os import path
from time import sleep
from urllib.parse import urlparse
from glob import glob

import click
import click_log
import requests

TASKS_URL = '{}/api/v2/tasks'

logger = logging.getLogger(__name__)


def run_with_slurm(sess, server, task, task_dir):
    subprocess.check_call(['sbatch', task['settings']['cmd']],
                          cwd=task_dir)


def download_only(sess, server, task, task_dir):
    pass


def run_direct(sess, server, task, task_dir):
    task['data']['commands'] = {}
    commands = task['data']['commands']

    # since we block below, we set the task to running right away
    req = sess.patch(server + task['_links']['self'],
                     json={'status': 'running'})
    req.raise_for_status()

    for entry in task['settings']['commands']:
        name = entry['name']
        stdout_fn = path.join(task_dir, "{}.out".format(name))
        stderr_fn = path.join(task_dir, "{}.err".format(name))

        commands[name] = {'return': -1}

        logger.info("task %s: running command %s", task['id'], name)

        try:
            with open(stdout_fn, 'w') as stdout, \
                 open(stderr_fn, 'w') as stderr:
                cproc = subprocess.run(
                    [entry['cmd']] + entry['args'],
                    stdout=stdout, stderr=stderr,
                    cwd=task_dir)

                commands[name]['return'] = cproc.returncode

                if not entry.get('ignore_failure', False):
                    cproc.check_returncode()

        except EnvironmentError as error:
            commands[name]['error'] = (
                "error when opening stdout/stderr files: {}"
                .format(error))
            raise
        except subprocess.SubprocessError as error:
            commands[name]['error'] = (
                "error when running: {}"
                .format(error))
            raise


RUNNERS = {
    'slurm': run_with_slurm,
    'download-only': download_only,
    'direct': run_direct,
    }


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

        try:
            runner = RUNNERS[task['settings']['machine']['runner']]
        except KeyError:
            raise NotImplementedError(
                "runner {} is not (yet) implemented".format(runner))

        if not task['data']:
            task['data'] = {}

        try:
            runner(sess, server, task, task_dir)

            for a_name in task['settings']['output_artifacts']:
                full_name = path.join(task_dir, a_name)

                filepaths = glob(full_name)

                if not filepaths:
                    logger.warning("task %s: glob for '%s' returned 0 files",
                                   task['id'], a_name)
                    # TODO: record this failure in the task
                    continue

                for filepath in filepaths:
                    data = {
                        'name': path.relpath(filepath, task_dir),
                        }
                    with open(filepath, 'rb') as data_fh:
                        req = sess.post(server + task['_links']['uploads'],
                                        data=data, files={'data': data_fh})
                        req.raise_for_status()

            req = sess.patch(server + task['_links']['self'],
                             json={'status': 'done', 'data': task['data']})
            req.raise_for_status()

        except requests.exceptions.HTTPError as error:
            logger.error("task %s: HTTP error occurred: %s\n%s",
                         task['id'], error, error.response.text)

        except Exception as error:
            logger.error("task %s: error occurred during run: %s",
                         task['id'], error)
            req = sess.patch(server + task['_links']['self'],
                             json={'status': 'error', 'data': task['data']})
            req.raise_for_status()


if __name__ == '__main__':
    main()
