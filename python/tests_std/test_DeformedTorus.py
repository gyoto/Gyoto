import numpy
import unittest
import gyoto
gyoto.loadPlugin('stdplug')
import gyoto_std

class TestDeformedTorus(unittest.TestCase):

    def test_DeformedTorus(self):
        ao=gyoto.Astrobj("DeformedTorus")
        ao=gyoto_std.DeformedTorus()
        sp=gyoto_std.BlackBody()
        ao.spectrum(sp)
        r=4.
        ao.largeRadius(r)
        self.assertTrue((ao.largeRadius() == r))
