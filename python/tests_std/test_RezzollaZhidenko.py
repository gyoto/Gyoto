import numpy
import unittest
import gyoto
gyoto.loadPlugin('stdplug')
import gyoto_std

class TestRezzollaZhidenko(unittest.TestCase):

    def test_mass(self):
        gg=gyoto_std.RezzollaZhidenko()
        m=2.
        gg.mass(m)
        self.assertTrue((gg.mass() == m))

    def test_aparam(self):
        gg=gyoto_std.RezzollaZhidenko()
        a=(1., 2., 3., 4.)
        gg.aparam(a)
        self.assertTrue((gg.aparam() == a))
