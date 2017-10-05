
from . import app, calculation_finished

try:
    from slackclient import SlackClient
    SLACK_TOKEN = app.config['SLACK_API_TOKEN']
except (ImportError, KeyError):
    pass
else:
    sc = SlackClient(SLACK_TOKEN)

    def slack_calc_finished(sender, task, calculation, **extras):
        text = "task {} for calculation {} on structure {} finished: {}".format(
                    task.id, calculation.id,
                    calculation.structure.name,
                    task.status.name.upper())

        sc.api_call(
            "chat.postMessage",
            channel=app.config.get('SLACK_CHANNEL', "#fatman"),
            text=text
        )

    calculation_finished.connect(slack_calc_finished)
