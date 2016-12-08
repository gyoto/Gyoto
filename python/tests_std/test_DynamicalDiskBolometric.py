import numpy
import unittest
import gyoto
gyoto.loadPlugin('stdplug')
import gyoto_std

class TestDynamicalDiskBolometric(unittest.TestCase):

    def test_DynamicalDiskBolometric(self):
        ao=gyoto.Astrobj("DynamicalDiskBolometric")
        ao=gyoto_std.DynamicalDiskBolometric()
        self.assertTrue(True)
