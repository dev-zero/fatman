
import bz2

from celery.utils.log import get_task_logger
from celery import group

from . import capp, resultfiles, tools, db
from .models import Result, Task, Method

logger = get_task_logger(__name__)


@capp.task
def postprocess_result_files(update=False):
    """Postprocess all results files.
    Returns a tuple of rids and a group task."""

    query = (Result.select(Result.id)
             .where(~(Result.filename >> None)))

    rids = [r.id for r in query]

    return (rids, group(postprocess_result_file.s(r, update)
                        for r in rids)())


@capp.task
@db.atomic()
def postprocess_result_file(rid, update=False):
    """Parse the file uploaded for a result.
    Set update to True to re-parse the file and update already present data"""

    result = (Result.select(Result, Task, Method)
              .where(Result.id == rid)
              .join(Task).join(Method)
              .get())

    if result.data and not update:
        logger.warning(("Parsed data available and update=False, "
                        "skipping file parsing for result id %d"),
                       rid)
        return False

    if not result.filename:
        logger.error("Missing filename for result id %d", rid)
        return False

    filepath = resultfiles.path(result.filename)
    code = result.task.method.code

    with bz2.open(filepath, mode='rt') as infile:
        try:
            data = tools.get_data_from_output(infile, code)
        except tools.OutputParseError as exc:
            logger.error("Parsing %s for result id %d failed: %s",
                         filepath, rid, exc)
            return False

    logger.info("Parsing %s for result id %d succeeded", filepath, rid)

    if result.data is None:
        result.data = data
    else:
        # merging works only for the top-level dict,
        # nested dicts/lists get overwritten
        result.data.update(data)

    result.save()

    return True
