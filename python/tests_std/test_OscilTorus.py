import numpy
import unittest
import gyoto
import gyoto_std

class TestOscilTorus(unittest.TestCase):

    def test_OscilTorus(self):
        ao=gyoto.Astrobj("OscilTorus")
        ao=gyoto_std.OscilTorus()
        r=4.
        ao.largeRadius(r)
        self.assertTrue((ao.largeRadius() == r))
        perturbkind="Vertical"
        ao.perturbKind(perturbkind)
        self.assertTrue((ao.perturbKind() == perturbkind))
