import numpy
import unittest
import gyoto
gyoto.loadPlugin('stdplug')
import gyoto_std

class TestInflateStar(unittest.TestCase):

    def test_InflateStar(self):
        ao=gyoto.Astrobj("InflateStar")
        ao=gyoto_std.InflateStar()
        ao.radius(0.)
        ao.radiusStop(2.)
        ao.timeInflateInit(0.)
        ao.timeInflateStop(2.)
        self.assertTrue(abs(ao.radiusAt(1.)-1.) < 1e-6)
