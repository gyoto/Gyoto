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
