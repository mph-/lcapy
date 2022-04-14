from lcapy import *
import unittest


class LcapyTester(unittest.TestCase):

    """Unit tests for lcapy

    """

    def test_pole_zero_plot(self):

        H = s / (s + 2)
        ax = H.plot(xlabel='real', ylabel='imag', title='test')
        self.assertEqual(ax.get_xlabel(), 'real', "x_label incorrect.")
        self.assertEqual(ax.get_ylabel(), 'imag', "y_label incorrect.")
        self.assertEqual(ax.get_title(), 'test', "title incorrect.")

    def test_time_plot(self):

        h = cos(2 * pi * t)
        ax = h.plot(xlabel='time', ylabel='voltage', title='test')
        self.assertEqual(ax.get_xlabel(), 'time', "x_label incorrect.")
        self.assertEqual(ax.get_ylabel(), 'voltage', "y_label incorrect.")
        self.assertEqual(ax.get_title(), 'test', "title incorrect.")

    def test_frequency_plot(self):

        h = cos(2 * pi * t)
        H = h(f)
        ax = h.plot(xlabel='frequency',
                    ylabel='voltage spectral density', title='test')
        self.assertEqual(ax.get_xlabel(), 'frequency', "x_label incorrect.")
        self.assertEqual(
            ax.get_ylabel(), 'voltage spectral density', "y_label incorrect.")
        self.assertEqual(ax.get_title(), 'test', "title incorrect.")

    def test_discrete_time_plot(self):

        h = cos(2 * pi * n * 5)
        ax = h.plot(xlabel='time', ylabel='voltage', title='test')
        self.assertEqual(ax.get_xlabel(), 'time', "x_label incorrect.")
        self.assertEqual(ax.get_ylabel(), 'voltage', "y_label incorrect.")
        self.assertEqual(ax.get_title(), 'test', "title incorrect.")
