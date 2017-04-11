import numpy
import unittest
import gyoto
import gyoto_std

class TestEquatorialHotSpot(unittest.TestCase):

    def test_EquatorialHotSpot(self):
        ao=gyoto.Astrobj("EquatorialHotSpot")
        ao=gyoto_std.EquatorialHotSpot()
        r=4.
        ao.spotRadSize(r)
        self.assertTrue((ao.spotRadSize() == r))
        kind="RadialBeaming"
        ao.beaming(kind)
        self.assertTrue((ao.beaming() == kind))
        self.assertRaises(gyoto.Error, ao.beaming, "foo")
