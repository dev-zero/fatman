
import json
import jsonschema
from os import path

from flask_testing import TestCase as BaseTestCase

from fatman import app

TASK = "c0735f0b-78c2-4deb-88ef-b307904d6c8c"
STRUCTURE = "2f56a08f-2e13-478c-94e6-b9430a99a890"
METHOD = "3fd3acad-0b20-4038-a0d5-ce49380712ce"
METHOD2 = "2ad699d5-1878-4a32-908e-dee3d7a32cd6"
RESULT = "a1b6e8ae-b7d9-4e1e-a8ae-f6ba0eb17c5e"
PSEUDO = "0d87c358-28f8-4e6e-931a-85d9aa3cbb26"


def setUpModule():
    app.config['TESTING'] = True
    app.config['PRESERVE_CONTEXT_ON_EXCEPTION'] = False


def tearDownModule():
    pass


class TestCase(BaseTestCase):
    schema_base = path.join(
        path.dirname(path.abspath(__file__)),
        'schemas')

    def create_app(self):
        return app

    def assertJSONSchema(self, data, schema_file=None):
        schema_path = path.join(self.schema_base,
                                schema_file if schema_file
                                else self.schema_file)

        with open(schema_path, 'r') as schema_fh:
            schema = json.load(schema_fh)
            resolver_path = 'file://{}/'.format(self.schema_base)
            resolver = jsonschema.RefResolver(resolver_path, schema)
            (jsonschema
             .Draft4Validator(schema, resolver=resolver)
             .validate(data))


class TestStructures(TestCase):
    """Tests for the /structure endpoint"""

    def test_get_by_id(self):
        """single structure by id"""
        resp = self.client.get('/structure?id=%s' % STRUCTURE)
        self.assertEqual(resp.status_code, 200)

    def test_get_by_test(self):
        """single structure by test"""
        resp = self.client.get('/structure?test=deltatest_H')
        self.assertEqual(resp.status_code, 200)


class TestStats(TestCase):
    """Tests for the /stats endpoint"""

    def test_get_tasks(self):
        """task stats"""
        resp = self.client.get('/stats/tasks')
        self.assertEqual(resp.status_code, 200)
        self.assertJSONSchema(resp.json, "stats_tasks.json")


class TestMachinestatus(TestCase):
    """Tests for the /machinestatus endpoint"""

    def test_get(self):
        """machine stats"""
        resp = self.client.get('/machinestatus')
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, list)
        self.assertIsInstance(resp.json[0], list)


class TestTasks(TestCase):
    """Tests for the /tasks endpoint"""

    schema_file = "tasks.json"

    def test_get(self):
        """task list"""
        resp = self.client.get('/tasks?limit=3')
        self.assertEqual(resp.status_code, 200)
        self.assertJSONSchema(resp.json)


class TestTask(TestCase):
    """Tests for the /tasks endpoint"""

    schema_file = "task.json"

    def test_get(self):
        """single task"""
        resp = self.client.get('/tasks/%s' % TASK)
        self.assertEqual(resp.status_code, 200)
        self.assertJSONSchema(resp.json)

    def test_patch(self):
        """change task"""
        resp = self.client.patch('/tasks/%s?status=done&machine=-&priority=100'
                                 % TASK)
        self.assertEqual(resp.status_code, 200)
        self.assertJSONSchema(resp.json)

    def test_post_multiple(self):
        """create tasks for method"""
        resp = self.client.post('/tasks?method=%s&test=deltatest_H' % METHOD)
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, list)

    def test_post_single(self):
        """create single task for method and structure"""
        resp = self.client.post('/tasks?method=%s&test=deltatest_H&structure=deltatest_H_1.00' % METHOD)
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, list)
        self.assertEqual(len(resp.json), 1)


class TestPlot(TestCase):
    """Tests for the /plot endpoint"""

    def test_get(self):
        """single plot"""
        resp = self.client.get('/plot?test=deltatest_H&method=%s' % METHOD)
        self.assertEqual(resp.status_code, 200)


class TestResult(TestCase):
    def test_get(self):
        """single plot"""
        resp = self.client.get('/plot?test=deltatest_H&method=%s' % METHOD)
        self.assertEqual(resp.status_code, 200)


class TestResult(TestCase):
    def test_get(self):
        """single result"""
        resp = self.client.get('/results/%s' % RESULT)
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, dict)

    def test_patch(self):
        """change result"""
        resp = self.client.patch('/results/%s' % RESULT)
        self.assertEqual(resp.status_code, 200)

    def test_post(self):
        """create result for task"""
        data = json.dumps({
            'task_id': TASK,
            'energy': 0.,
            'data': json.dumps({'foo': 'bar'}),
        })
        resp = self.client.post('/results', data=data,
                                content_type='application/json')
        self.assertEqual(resp.status_code, 201)

        resp = self.client.get(resp.location)
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, dict)


    def test_get_list(self):
        """result list"""
        resp = self.client.get('/results?code=ref-wien2k')
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, list)


class TestBasis(TestCase):
    def test_get(self):
        """single basis set"""
        resp = self.client.get('/basis?element=H&family=TZVP-GTH')
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, dict)


class TestPseudos(TestCase):
    def test_get(self):
        """single pseudo"""
        resp = self.client.get('/pseudos/%s' % PSEUDO)
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, dict)

    def test_get_list(self):
        """get pseudo list"""
        resp = self.client.get('/pseudos?format=CP2K&family=GTH-PBE&element=H')
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, list)

    # TODO: implement pseudo post

class TestPseudoFamilies(TestCase):
    def test_get_list(self):
        """pseudofamily list"""
        resp = self.client.get('/pseudofamilies')
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, list)


class TestResults(TestCase):
    def test_get(self):
        """testresult list by method and test"""
        resp = self.client.get('/testresults?method=%s&test=deltatest_H' % METHOD)
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, list)


class TestCompare(TestCase):
    def test_get(self):
        """comparison (2 methods)"""
        resp = self.client.get('/compare?method1=%s&method2=%s' % (METHOD, METHOD2))
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, dict)

    # TODO: implement 1 method comparison and methodbytest

class TestMethods(TestCase):
    def test_get(self):
        """single method"""
        resp = self.client.get('/methods/%s' % METHOD)
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, dict)

    def test_get_list(self):
        """method list"""
        resp = self.client.get('/methods')
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, list)

    def test_post(self):
        """create method task"""
        data = json.dumps({
            'code': "CP2K",
            'pseudopotential': "GTH-PBE",
            'basis_set': "TZVP-GTH",
            'settings': json.dumps({'foo': 'bar'}),
        })
        resp = self.client.post('/methods', data=data,
                                content_type='application/json')
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, dict)


class TestTests(TestCase):
    def test_get(self):
        """single test"""
        resp = self.client.get('/tests/deltatest_H')
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, list)

    def test_get_list(self):
        """test list"""
        resp = self.client.get('/tests')
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, list)
