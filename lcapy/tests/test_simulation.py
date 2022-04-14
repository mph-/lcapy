from lcapy import *
import numpy as np
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_sim(self):

        cct = Circuit("""
        V1 1 0 step 10; down
        R1 1 2 5; right
        L1 2 0_2 2; down
        W 0 0_2; right""")

        results = cct.sim(np.linspace(0, 10))

        self.assertEqual(results.R1.v[0], 0, 'vR[0]')
        self.assertEqual(np.round(results.R1.v[-1], 3), 10, 'vR[-1]')
        self.assertEqual(np.round(results.R1.i[-1], 3), 2, 'iR[-1]')
        self.assertEqual(np.round(results.L1.v[-1], 3), 0, 'vL[0]')
        self.assertEqual(np.round(results.L1.i[-1], 3), 2, 'iL[-1]')
