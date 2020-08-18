from lcapy import *
from lcapy.polyphase import *
from lcapy.discretetime import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_polyphase(self):

        A = symbol('A')
        alpha = polyphase_alpha(3)
        
        V = PhaseVoltageVector((A, A * alpha, A * alpha**2))
