import unittest
import gyoto
import numpy

gyoto.loadPlugin('stdplug')

class TestMetric(unittest.TestCase):

    def test_circularVelocity(self):
        gg=gyoto.Metric('KerrBL')
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
        gg=gyoto.Metric('KerrBL')
        self.assertAlmostEqual(gg.SysPrimeToTdot([0, 6, numpy.pi/2, 0], (0, 0, 0.1)), 1.8057877962865378)

    def test_isStopCondition(self):
        gg=gyoto.Metric('KerrBL')
        self.assertFalse(gg.isStopCondition((0, 6, 3.14, 0, 0, 0, 0, 0)))
        self.assertTrue(gg.isStopCondition((0, 0, 3.14, 0, 0, 0, 0, 0)))

    def test_gmunu(self):
        gg=gyoto.Metric('KerrBL')
        tt=gg.gmunu((0, 6, 3.14, 0), 0, 0)
        self.assertAlmostEqual(tt, -0.6666666666666667)
        dst=numpy.zeros((4, 4), float)
        gg.gmunu(dst, (0, 6, 3.14, 0))
        self.assertEqual(tt, dst[0, 0])
        dst2=gg.gmunu((0, 6, 3.14, 0))
        self.assertEqual(tt, dst2[0, 0])

    def test_christoffel(self):
        gg=gyoto.Metric('KerrBL')
        tt=gg.christoffel((0, 6, 3.14, 0), 0, 0, 0)
        self.assertAlmostEqual(tt, 0)
        dst=numpy.zeros((4, 4, 4), float)
        gg.christoffel(dst, (0, 6, 3.14, 0))
        self.assertEqual(tt, dst[0, 0, 0])
        dst2=gg.christoffel((0, 6, 3.14, 0))
        self.assertEqual(tt, dst2[0, 0, 0])
