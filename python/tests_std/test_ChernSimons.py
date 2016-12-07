import numpy
import unittest
import gyoto
gyoto.loadPlugin('stdplug')
import gyoto_std

class TestChernSimons(unittest.TestCase):

    def test_setInitCoord(self):
        gg=gyoto_std.ChernSimons()
        zeta=0.5
        gg.dzetaCS(zeta)
        self.assertTrue((gg.dzetaCS() == zeta))
