import unittest
import gyoto

class TestUnit(unittest.TestCase):

    def test___str__(self):
        self.assertEqual(str(gyoto.Unit("km")), "km")

    def test_throw(self):
        self.assertRaises(gyoto.Error, lambda: gyoto.Unit("deadbeaf"))

    def test_to(self):
        self .assertEqual(gyoto.Unit("m").To(5, gyoto.Unit("km")), 5000.)

    def test_from(self):
        self .assertEqual(gyoto.Unit("m").From(5000, gyoto.Unit("km")), 5.)

class TestConverter(unittest.TestCase):

    def test_trivial(self):
        a=gyoto.Converter()
        self.assertEqual(a(2.), 2.)

    def test_converter(self):
        a=gyoto.Converter(gyoto.Unit("m"), gyoto.Unit("km"))
        self.assertEqual(5., a(5000))

    def test_reset_trivial(self):
        a=gyoto.Converter(gyoto.Unit("m"), gyoto.Unit("km"))
        a.reset()
        self.assertEqual(a(2.), 2.)

    def test_reset(self):
        a=gyoto.Converter()
        a.reset(gyoto.Unit("m"), gyoto.Unit("km"))
        self.assertEqual(5., a(5000))
