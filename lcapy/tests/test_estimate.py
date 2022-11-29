from lcapy import *
from numpy import arange
import unittest


class LcapyTester(unittest.TestCase):
    """Unit tests for lcapy"""

    def test_estimate1(self):

        e = expr('a * exp(-t  / tau) * u(t)')

        tv = arange(10)

        params = {'a': 1, 'tau': 10}

        vv = e.subs(params).evaluate(tv)

        eparams = e.estimate(tv, vv, method='trf',
                             ranges={'a': (0, 10), 'tau': (1, 20)}).params

        assert abs(params['a'] - eparams['a']) < 1e-6 and \
            abs(params['tau'] - eparams['tau']) < 1e-6

        eparams = e.estimate(tv, vv, method='Nelder-Mead',
                             ranges={'a': (0, 10), 'tau': (1, 20)}).params

        assert abs(params['a'] - eparams['a']) < 1e-5 and \
            abs(params['tau'] - eparams['tau']) < 1e-4

        eparams = e.estimate(tv, vv, method='brute',
                             ranges={'a': (0, 10), 'tau': (1, 20)}).params

        assert abs(params['a'] - eparams['a']) < 1e-5 and \
            abs(params['tau'] - eparams['tau']) < 1e-4
