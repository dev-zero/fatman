
import bz2
from collections import OrderedDict
from urllib.parse import urlsplit

from ase.units import kcal, mol

from celery.utils.log import get_task_logger
from celery import group

from sqlalchemy.orm import joinedload, contains_eager

from . import capp, resultfiles, tools, db, cache
from .models import (
    Result,
    Task,
    Method,
    Test,
    TestResult,
    Task2,
    Calculation,
    Artifact,
    TaskStatus,
    TestResult2,
    TestResult2Collection,
    Pseudopotential,
    )

from .tools import Json2Atoms, deltatest_ev_curve, gmtkn_coefficients

logger = get_task_logger(__name__)


@capp.task
def postprocess_result_files(update=False):
    """Postprocess all results files.
    Returns a tuple of rids and a group task."""

    query = (db.session.query(Result.id)
             .filter(Result.filename != None)
             .all())

    rids = [r.id for r in query]

    return (rids, group(postprocess_result_file.s(r, update)
                        for r in rids)())


@capp.task
def postprocess_result_file(rid, update=False):
    """Parse the file uploaded for a result.
    Set update to True to re-parse the file and update already present data"""

    result = (Result.query
              .options(joinedload("task")
                       .joinedload("method"))
              .get(rid))

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

    db.session.commit()

    return True


@capp.task(bind=True)
def calculate_deltatest(self, tid, rid, trid=None):
    """Calculate the Δ-value for a given test ID `tid`,
    including the result with ID `rid`.
    If `trid` is given, update the respective TestResult instead of
    creating a new entry"""

    test = Test.query.get(tid)
    method = (Method.query
              .join(Task)
              .join(Result)
              .filter(Result.id == rid)
              .one())

    # select all results which share the same method
    # and are from the same test (deltatest_<El>_<Vol>)
    results = (Result.query
               .join(Task)
               .join(Method)
               .join(Test)
               .options(contains_eager(Result.task))
               .filter(Method.id == method.id, Test.id == test.id))

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

    for result in results:
        if 'total_energy' not in result.data.keys():
            logger.info(("total_energy missing for %s to calculate %s value"
                         " for Method %s from result %s, skipping"),
                        result.id, test.name, method.id, rid)
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

            energies.append(result.data['total_energy']/natom)
            volumes.append(struct.get_volume()/natom)

        v, e, B0, B1, R = deltatest_ev_curve(volumes, energies)

        if isinstance(v, complex) or (v == "fail"):
            result_data = {"_status": "unfittable"}
        else:
            result_data = {"_status": "fitted",
                           "V": v, "E0": e, "B0": B0, "B1": B1, "R": R}

        if trid:
            testresult = TestResult.query.get(trid)
            testresult.result_data = result_data
        else:
            testresult = TestResult(test=test,
                                    method=method,
                                    result_data=result_data)
        db.session.commit()

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

    test = Test.query.get(tid)
    method = (Method.query
              .join(Task)
              .join(Result)
              .filter(Result.id == rid)
              .one())

    # select all results which share the same method
    # and are from the same test (GMTKN_<sub-set>)
    results = (Result.query
               .join(Task)
               .join(Method)
               .join(Test)
               .options(contains_eager(Result.task))
               .filter(Method.id == method.id, Test.id == test.id))

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

    for result in results:
        if 'total_energy' not in result.data.keys():
            logger.info(("total_energy missing for %s to calculate %s value"
                         " for Method %s from result %s, skipping"),
                        result.id, test.name, method.id, rid)
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
        struct2energies = OrderedDict([(x.task.structure.name,
                                        x.data['total_energy'])
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
            testresult = TestResult.query.get(trid)
            testresult.result_data = result_data
        else:
            testresult = TestResult(test=test,
                                    method=method,
                                    result_data=result_data)
        db.session.commit()

        logger.info("Calculated GMTKN values for %s using Method %s",
                    test.name, method.id)

    finally:
        cache.cache.delete(lock_id)

    return True


@capp.task
def postprocess_result(rid, update=False):
    # check whether there is already testresult using this result
    # and the same method
    trid = (db.session.query(TestResult.id)
            .join(Test)
            .join(Task)
            .join(Result)
            .filter(Task.method_id == TestResult.method_id,
                    Result.id == rid)
            .scalar())

    if trid:
        if update:
            logger.info("Updating TestResult %s for result %s",
                        trid, rid)
        else:
            logger.warning(("Found TestResult %s for result %s"
                            ", skipping postprocessing"),
                           trid, rid)
            return False

    test = (Test.query
            .join(Task)
            .join(Result)
            .filter(Result.id == rid)
            .one())

    if test.name.startswith('deltatest'):
        return calculate_deltatest(test.id, rid, trid)
    elif test.name.startswith('GMTKN'):
        return calculate_GMTKN(test.id, rid, trid)
    else:
        logger.error("No TestResult process function found for Test %s",
                     test.name)
        return False


@capp.task
def generate_calculation_results(calc_id, update=False):
    """
    Parse the outputs attached to the first succeeded task for the given calculation,
    if the calculation does not have any results yet.

    Args:
        update: update the results if they already exist, careful: test results
                based on this calculation are not automatically updated!
    """
    calc = (Calculation.query
            .options(joinedload('code'))
            .get(calc_id))

    if calc is None:
        logger.error("calculation %s: not found")
        return False

    if calc.results and not update:
        logger.info("calculation %s: has already a result and update=False, skipping", calc.id)
        return False

    task = (calc.tasks_query
            .filter(TaskStatus.name == 'done')
            .order_by(Task2.mtime.desc())
            .first())

    if task is None:
        logger.info("calculation %s: no successfully finished task found, skipping", calc.id)
        return False

    logger.info("calculation %s: taking results from task %s", calc.id, task.id)

    if calc.code.name == 'CP2K':
        resultfile_name = 'calc.out'
    else:
        logger.error("calculation %s: no parser available (yet) for code %s",
                     calc.id, task.calculation.code.name)
        return False

    artifact = (task.outfiles
                .filter(Artifact.name == resultfile_name)
                .one_or_none())

    if artifact is None:
        logger.info("calculation %s: no artifact with name '%s' linked to task %s",
                    calc.id, resultfile_name, task.id)
        return False

    scheme, nwloc, path, _, _ = urlsplit(artifact.path)

    results = None

    if scheme == 'fkup' and nwloc == 'results':
        filepath = resultfiles.path(path[1:])

        # Artifact.name contains a full path, including possible subdirs
        compressed = artifact.mdata.get('compressed', None)

        if compressed is None:
            with open(filepath, 'r') as fhandle:
                results = tools.get_data_from_output(fhandle, calc.code.name)
        elif compressed == 'bz2':
            with bz2.open(filepath, mode='rt') as fhandle:
                results = tools.get_data_from_output(fhandle, calc.code.name)
        else:
            logger.error("unsupported compression scheme found: %s", compressed)
            return False

    else:
        logger.error("unsupported storage scheme found: %s", scheme)
        return False

    if not results:
        logger.error("calculation %s: no parseable data found in artifact %s",
                     calc.id, artifact.id)
        return False

    task.calculation.results = results
    db.session.commit()

    return True


@capp.task(bind=True)
def generate_test_result_deltatest(self, calc_id, update=False, force_new=False):
    """
    Calculate the Δ-test value using a given calculation id if one does not exist already (see arguments).
    Creates a new Test Result in a Test Result Collection named
    after the Calculation Collection of the given calculation.

    Args:
        calc_id: the UUID of the calculation to be included in the generated test result
        update: whether or not to update an already existing test result for the same collection and test
        force_new: always create a new test result, irrespective of already existing ones
    """

    calc = (Calculation.query
            .options(joinedload('test'))
            .options(joinedload('collection'))
            .get(calc_id))

    tresult = (calc.testresults_query
               .filter(TestResult2.test == calc.test)
               .join(TestResult2.collections)
               .filter(TestResult2Collection.name == calc.collection.name)
               .first())

    if tresult is not None and not (force_new or update):
        logger.info("calculation %s: a test result for this test and collection already exists, skipping", calc.id)
        return False

    if force_new:
        tresult = None

    # getting the element from the first (and only) pseudo
    element = calc.pseudos[0].element

    result_data = {'element': element}

    # from the same calculation collection as the given calculation,
    # get the calculations for the same elements (selected via pseudos)
    calcs = (Calculation.query
             .filter_by(collection=calc.collection)
             .options(joinedload('structure'))
             .join(Calculation.pseudos)
             .filter(Pseudopotential.element == element)
             .limit(6)  # limit to one more than we can actually use, to catch errors
             .all())

    if len(calcs) < 5:
        logger.info("Not enough datapoints to calculate deltatest value using calculation %s", calc.id)
        return False

    if len(calcs) > 5:
        logger.critical("More than 5 elements found to calculate deltatest value using calculation %s", calc.id)
        return False

    for calc in calcs:
        if not calc.results or 'total_energy' not in calc.results.keys():
            logger.error("total_energy missing for calculation %s to calculate deltatest value", calc.id)
            return False

    lock_id = '{}-lock-{}-collection-{}'.format(self.name, calc.test.name, calc.collection.id)
    if not cache.cache.add(lock_id, True, 60*60):
        # another task is already calculating the deltavalue for this result
        logger.info("Skipping deltatest calculation for %s, another calculation task is already running", calc.id)
        return False

    try:
        energies = []
        volumes = []
        natom = 0

        for calc in calcs:
            struct = Json2Atoms(calc.structure.ase_structure)
            natom = len(struct.get_atomic_numbers())

            energies.append(calc.results['total_energy']/natom)
            volumes.append(struct.get_volume()/natom)

        result_data.update({
            'energies': energies,
            'volumes': volumes,
            })

        coeffs = deltatest_ev_curve(volumes, energies)

        if isinstance(coeffs[0], complex) or (coeffs[0] == "fail"):
            result_data['status'] = "unfittable"
        else:
            result_data.update({
                'status': "fitted",
                'coefficients': dict(zip(('v', 'e', 'B0', 'B1', 'R'), coeffs)),
                })

        if not tresult:
            # get or create a new test result collection with the same name as the calc collection
            trcollection = (TestResult2Collection.query
                            .filter_by(name=calc.collection.name)
                            .first())

            if not trcollection:
                trcollection = TestResult2Collection(name=calc.collection.name)

            tresult = TestResult2(test=calc.test, data=result_data)
            tresult.collections.append(trcollection)
            tresult.calculations.extend(calcs)

            db.session.add(tresult)

        else:
            tresult.data = result_data
            tresult.calculations.clear()
            tresult.calculations.extend(calcs)

        db.session.commit()

        logger.info("Calculated deltatest value for element %s in collection %s",
                    element, calc.collection.id)

    finally:
        cache.cache.delete(lock_id)

    # only if no exception was thrown somewhere above
    return True

#  vim: set ts=4 sw=4 tw=0 :
