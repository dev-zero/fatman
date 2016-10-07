#!/usr/bin/env python

import time
import requests
import click


# define an empty command group
@click.group()
@click.pass_context
@click.option('--url', type=str, default='https://tctdb.chem.uzh.ch/fatman',
              show_default=True, help="The URL where FATMAN is running")
@click.option('--watch/--no-watch', '-w', default=False, show_default=True,
              help="Watch the actions")
def cli(ctx, url, watch):
    """Trigger actions on FATMAN"""
    ctx.obj['base'] = url
    ctx.obj['watch'] = watch


@cli.group()
@click.pass_context
@click.option('--update/--no-update', default=False, show_default=True,
              help="Update already existing results")
@click.option('--id', '-i', 'rid', type=click.UUID,
              help="result ID to change (or all if empty)")
def results(ctx, update, rid):
    """Call actions for results"""
    ctx.obj['update'] = update
    ctx.obj['rid'] = rid


@results.command()
@click.pass_context
def postproc(ctx):
    """Trigger postprocessing on given result or all results"""

    sess = requests.Session()
    # the certificate is signed by the inofficial TC-Chem CA
    sess.verify = False

    action_data = {'doPostprocessing': {'update': ctx.obj['update'], }}

    if ctx.obj['rid']:
        req_url = '{base}/api/v1/results/{rid}/action'.format(
            base=ctx.obj['base'], rid=ctx.obj['rid'])
    else:
        req_url = '{base}/api/v1/results/action'.format(
            base=ctx.obj['base'])

    req = sess.post(req_url, json=action_data)
    req.raise_for_status()

    response_url = req.headers['location']
    click.echo("Result URL:\n  {}\n".format(response_url))

    if ctx.obj['watch']:
        status = 202
        while status != 200:
            req = sess.get(response_url)
            req.raise_for_status()

            status = req.status_code

            if req.text:
                click.echo(req.json())

            time.sleep(1)

if __name__ == "__main__":
    cli(obj={})
