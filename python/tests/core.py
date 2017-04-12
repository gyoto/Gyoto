import unittest
import gyoto.core
import numpy

gyoto.core.requirePlugin('stdplug')

class TestUnit(unittest.TestCase):

    def test___str__(self):
        self.assertEqual(str(gyoto.core.Unit("km")), "km")

    def test_throw(self):
        self.assertRaises(gyoto.core.Error, lambda: gyoto.core.Unit("deadbeaf"))

    def test_to(self):
        self .assertEqual(gyoto.core.Unit("m").To(5, gyoto.core.Unit("km")), 5000.)

    def test_from(self):
        self .assertEqual(gyoto.core.Unit("m").From(5000, gyoto.core.Unit("km")), 5.)

class TestConverter(unittest.TestCase):

    def test_trivial(self):
        a=gyoto.core.Converter()
        self.assertEqual(a(2.), 2.)

    def test_converter(self):
        a=gyoto.core.Converter(gyoto.core.Unit("m"), gyoto.core.Unit("km"))
        self.assertEqual(5., a(5000))

    def test_reset_trivial(self):
        a=gyoto.core.Converter(gyoto.core.Unit("m"), gyoto.core.Unit("km"))
        a.reset()
        self.assertEqual(a(2.), 2.)

    def test_reset(self):
        a=gyoto.core.Converter()
        a.reset(gyoto.core.Unit("m"), gyoto.core.Unit("km"))
        self.assertEqual(5., a(5000))

class TestError(unittest.TestCase):

    def test_get_message(self):
        a=gyoto.core.Error("toto")
        self.assertEqual(a.get_message(), "toto")

    def test_throwError(self):
        self.assertRaises(gyoto.core.Error, lambda: gyoto.core.throwError("msg"))

    def test_getErrcode(self):
        a=gyoto.core.Error("toto")
        self.assertEqual(a.getErrcode(), 1)

class TestMetric(unittest.TestCase):

    def test_circularVelocity(self):
        gg=gyoto.core.Metric('KerrBL')
        vel=numpy.zeros(4, float)
        vel2=numpy.zeros(4, float)
        pos_list=[0, 6, numpy.pi*0.5, 0]
        pos_numpy=numpy.asarray(pos_list)
        gg.circularVelocity(pos_numpy, vel)
        gg.circularVelocity(pos_list, vel2)
        vel3=gg.circularVelocity(pos_list)
        self.assertTrue((vel == vel2).all())
        self.assertTrue((vel == vel3).all())

    def test_SysPrimeToTdot(self):
        gg=gyoto.core.Metric('KerrBL')
        self.assertAlmostEqual(gg.SysPrimeToTdot([0, 6, numpy.pi/2, 0], (0, 0, 0.1)), 1.8057877962865378)

    def test_isStopCondition(self):
        gg=gyoto.core.Metric('KerrBL')
        self.assertFalse(gg.isStopCondition((0, 6, 3.14, 0, 0, 0, 0, 0)))
        self.assertTrue(gg.isStopCondition((0, 0, 3.14, 0, 0, 0, 0, 0)))

    def test_gmunu(self):
        gg=gyoto.core.Metric('KerrBL')
        tt=gg.gmunu((0, 6, 3.14, 0), 0, 0)
        self.assertAlmostEqual(tt, -0.6666666666666667)
        dst=numpy.zeros((4, 4), float)
        gg.gmunu(dst, (0, 6, 3.14, 0))
        self.assertEqual(tt, dst[0, 0])
        dst2=gg.gmunu((0, 6, 3.14, 0))
        self.assertEqual(tt, dst2[0, 0])

    def test_christoffel(self):
        gg=gyoto.core.Metric('KerrBL')
        tt=gg.christoffel((0, 6, 3.14, 0), 0, 0, 0)
        self.assertAlmostEqual(tt, 0)
        dst=numpy.zeros((4, 4, 4), float)
        gg.christoffel(dst, (0, 6, 3.14, 0))
        self.assertEqual(tt, dst[0, 0, 0])
        dst2=gg.christoffel((0, 6, 3.14, 0))
        self.assertEqual(tt, dst2[0, 0, 0])

class TestValue(unittest.TestCase):

    def test_good(self):
        a=gyoto.core.Value(5)
        self.assertEqual(a.toLong(), 5)
        self.assertEqual(a.toULong(), 5)
        self.assertEqual(a.toSizeT(), 5)
        a=gyoto.core.Value((1, 2, 3.))
        self.assertEqual(a.toVDouble(), (1., 2., 3.))
        a=gyoto.core.Value((1, 2, 3))
        self.assertEqual(a.toVULong(), (1, 2, 3))
        a=gyoto.core.Value('a')
        self.assertEqual(a.toString(), 'a')
        s=gyoto.core.Screen()
        a=gyoto.core.Value(s)
        self.assertEqual(a.toScreen().this, s.this)

    def test_assign(self):
        a=gyoto.core.Value(0)
        a.assign(gyoto.core.Value(5))
        self.assertEqual(a.toLong(), 5)
        a.assign(gyoto.core.Value(5.))
        self.assertEqual(a.toDouble(), 5.)

    def test_bad(self):
        self.assertRaises(gyoto.core.Error, lambda: gyoto.core.Value(5).toDouble())
        self.assertRaises(gyoto.core.Error, lambda: gyoto.core.Value(5.).toLong())
        self.assertRaises(gyoto.core.Error, lambda: gyoto.core.Value((1,)).toVDouble())
        self.assertRaises(gyoto.core.Error, lambda: gyoto.core.Value('a').toVULong())

class TestProperty(unittest.TestCase):

    def test_get(self):
        s=gyoto.core.Screen()
        self.assertEqual(s.get('Distance'), 1.0)
        self.assertRaises(gyoto.core.Error, lambda: s.get('NonExistentProperty'))

    def test_set(self):
        s=gyoto.core.Screen()
        s.set('Distance', 8., 'kpc')
        self.assertAlmostEqual(s.get('Distance'), 8.*gyoto.core.GYOTO_KPC, -15)

    def test_describe(self):
        s=gyoto.core.Screen()
        p=s.property('Distance')
        self.assertIn('Distance: double with unit', s.describeProperty(p))
