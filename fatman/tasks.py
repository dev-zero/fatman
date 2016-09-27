
import bz2
from datetime import datetime
from collections import OrderedDict

from celery.utils.log import get_task_logger
from celery import group

from ase.units import kcal, mol

from . import capp, resultfiles, tools, db, cache
from .models import Result, Task, Method, Test, TestResult
from .tools import Json2Atoms, deltatest_ev_curve, gmtkn_coefficients

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
        logger.warning(("Parsed data available and update=False"
                        ", skipping file parsing for result id %s"),
                       rid)
        return False

    if not result.filename:
        logger.error("Missing filename for result id %s", rid)
        return False

    filepath = resultfiles.path(result.filename)
    code = result.task.method.code

    with bz2.open(filepath, mode='rt') as infile:
        try:
            data = tools.get_data_from_output(infile, code)
        except tools.OutputParseError as exc:
            logger.error("Parsing %s for result id %s failed: %s",
                         filepath, rid, exc)
            return False

    logger.info("Parsing %s for result id %s succeeded", filepath, rid)

    if result.data is None:
        result.data = data
    else:
        # merging works only for the top-level dict,
        # nested dicts/lists get overwritten
        result.data.update(data)

    result.save()

    return True


@capp.task(bind=True)
def calculate_deltatest(self, tid, rid, trid=None):
    """Calculate the Δ-value for a given test ID `tid`,
    including the result with ID `rid`.
    If `trid` is given, update the respective TestResult instead of
    creating a new entry"""

    test = Test.get(Test.id == tid)
    method = (Method.select()
              .join(Task)
              .join(Result)
              .where(Result.id == rid)
              .get())

    # select all results which share the same method
    # and are from the same test (deltatest_<El>_<Vol>)
    results = (Result.select()
               .join(Task).join(Method)
               .switch(Task).join(Test)
               .where((Method.id == method.id) & (Test.id == tid)))

    if results.count() < 5:
        logger.info(("Not enough datapoints to calculate %s value"
                     " for Method %s using result %s, skipping"),
                    test.name, method.id, rid)
        return False

    if results.count() > 5:
        logger.critical(("More than 5 elements found to calculate %s value"
                         " for Method %s using result %s"),
                        test.name, method.id, rid)
        return False

    lock_id = '{}-lock-{}-method-{}'.format(self.name, test.name, method.id)
    if not cache.cache.add(lock_id, True, 60*60):
        # another task is already calculating the deltavalue for this result
        logger.info(("Skipping deltatest calculation for %s using Method %s"
                     ", another calculation task is already running"),
                    test.name, method.id)
        return False

    try:
        energies = []
        volumes = []
        natom = 0

        for result in results:
            struct = Json2Atoms(result.task.structure.ase_structure)
            natom = len(struct.get_atomic_numbers())

            energies.append(result.energy/natom)
            volumes.append(struct.get_volume()/natom)

        v, e, B0, B1, R = deltatest_ev_curve(volumes, energies)

        if isinstance(v, complex) or (v == "fail"):
            result_data = {"_status": "unfittable"}
        else:
            result_data = {"_status": "fitted",
                           "V": v, "E0": e, "B0": B0, "B1": B1, "R": R}

        if trid:
            testresult = TestResult.get(TestResult.id == trid)
            testresult.ctime = datetime.now()
            testresult.result_data = result_data
            testresult.save()
        else:
            TestResult.create(test=test,
                              method=method,
                              ctime=datetime.now(),
                              result_data=result_data)

        logger.info("Calculated deltatest value for %s using Method %s",
                    test.name, method.id)

    finally:
        cache.cache.delete(lock_id)

    # only if no exception was thrown somewhere above
    return True


@capp.task(bind=True)
def calculate_GMTKN(self, tid, rid, trid=None):
    """Calculate the GMTKN-value for a given test ID `tid`,
    including the result with ID `rid`.
    If `trid` is given, update the respective TestResult instead of
    creating a new entry"""

    test = Test.get(Test.id == tid)
    method = (Method.select()
              .join(Task)
              .join(Result)
              .where(Result.id == rid)
              .get())

    # select all results which share the same method
    # and are from the same test (GMTKN_<sub-set>)
    results = (Result.select()
               .join(Task).join(Method)
               .switch(Task).join(Test)
               .where((Method.id == method.id) & (Test.id == tid)))

    sub_db = test.name[6:]
    all_structures = set([item[0]
                          for sublist in gmtkn_coefficients[sub_db]
                          for item in sublist])

    results_count = results.count()
    required_count = len(all_structures)
    if results_count < required_count:
        logger.info(("Only %d of %d required results to calculate %s available"
                     ", skipping"), results_count, required_count, test.name)
        return False

    lock_id = '{}-lock-{}-method-{}'.format(self.name, test.name, method.id)
    if not cache.cache.add(lock_id, True, 60*60):
        # another task is already calculating the deltavalue for this result
        logger.info(("Skipping GMTKN calculation for %s using Method %s"
                     ", another calculation task is already running"),
                    test.name, method.id)
        return False

    try:
        result_data = {"_status": "incomplete", "energies": []}

        prefix = "gmtkn30_{}_".format(sub_db)
        struct2energies = OrderedDict([(x.task.structure.name, x.energy)
                                       for x in results])
        for rxn in gmtkn_coefficients[sub_db]:
            etot = 0.
            for struct, coeff in rxn:
                try:
                    etot += coeff*struct2energies[prefix+struct]
                except KeyError:
                    logger.warning(("Result for structure %s%s missing"
                                    ", ignored to calculate %s"),
                                   prefix, struct, test.name)

            logger.debug("%s: %f/%f", rxn, etot/kcal*mol, etot)
            result_data["energies"].append(etot)

        result_data["_status"] = "done"

        if trid:
            testresult = TestResult.get(TestResult.id == trid)
            testresult.ctime = datetime.now()
            testresult.result_data = result_data
            testresult.save()
        else:
            TestResult.create(test=test,
                              method=method,
                              ctime=datetime.now(),
                              result_data=result_data)

        logger.info("Calculated GMTKN values for %s using Method %s",
                    test.name, method.id)

    finally:
        cache.cache.delete(lock_id)

    return True


@capp.task
def postprocess_result(rid, update=False):
    testresult = (TestResult.select()
                  .join(Test).join(Task).join(Result)
                  .where((Task.method == TestResult.method) &
                         (Result.id == rid)))

    trid = None

    if testresult.count():
        trid = testresult.get().id
        if update:
            logger.info("Updating TestResult %s for result %s",
                        trid, rid)
        else:
            logger.warning(("Found TestResult %s for result %s"
                            ", skipping postprocessing"),
                           trid, rid)
            return False

    test = (Test.select(Test)
            .join(Task).join(Result)
            .where(Result.id == rid)
            .get())

    if test.name.startswith('deltatest'):
        return calculate_deltatest(test.id, rid, trid)
    elif test.name.startswith('GMTKN'):
        return calculate_GMTKN(test.id, rid, trid)
    else:
        logger.error("No TestResult process function found for Test %s",
                     test.name)
        return False

#  vim: set ts=4 sw=4 tw=0 :
