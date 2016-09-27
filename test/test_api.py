
import json

from flask_testing import TestCase as BaseTestCase

from fatman import app
from fatman.models import Structure

TASK = "ee408098-ff00-436c-b407-40d94055b84e"
STRUCTURE = "339a5953-5f49-4425-aa25-a85532f92e70"
METHOD = "a60ead16-e4d9-4fc2-b1f7-52dda03d61ac"
RESULT = "4e2857df-84bb-4c7d-95dd-f79f2e197fb2"
PSEUDO = "e0f435b6-da72-4560-aa19-6b2d5eb73f4d"


def setUpModule():
    app.config['TESTING'] = True
    app.config['PRESERVE_CONTEXT_ON_EXCEPTION'] = False


def tearDownModule():
    pass


class TestCase(BaseTestCase):
    def create_app(self):
        return app


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
        self.assertIsInstance(resp.json, list)
        self.assertIsInstance(resp.json[0], dict)


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

    def test_get_list(self):
        """task list"""
        resp = self.client.get('/tasks?limit=3')
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, list)

    def test_get(self):
        """single task"""
        resp = self.client.get('/tasks/%s' % TASK)
        self.assertEqual(resp.status_code, 200)
        self.assertIsInstance(resp.json, dict)

    def test_patch(self):
        """change task"""
        resp = self.client.patch('/tasks/%s?status=done&machine=-&priority=100'
                                 % TASK)
        self.assertEqual(resp.status_code, 200)

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
        resp = self.client.get('/compare?method1=%s&method2=2735a114-5031-4dbc-928b-98feb630df74' % METHOD)
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


class TestPseudos(TestCase):
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
