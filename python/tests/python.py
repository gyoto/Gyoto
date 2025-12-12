import unittest, numpy
import gyoto.python

class TestPythonSpectrum(unittest.TestCase):

    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sp = gyoto.python.PythonSpectrum()
        self.sp.Module = "gyoto_sample_spectra"

    def test_BlackBody6000(self):
        self.sp.Class = "BlackBody6000"
        stdsp=gyoto.std.BlackBody()
        stdsp.Temperature = 6000.
        stdsp.Scaling = 1.
        stdsp.ColorCorrection = 1.
        self.assertAlmostEqual(self.sp(1.)/stdsp(1.), 1., places=2)
        self.assertAlmostEqual(self.sp(1.e3)/stdsp(1.e3), 1., places=5)
        self.assertAlmostEqual(self.sp(1.e6)/stdsp(1.e6), 1., places=7)

    def test_PowerLaw(self):
        self.sp.Class = "PowerLaw"
        stdsp=gyoto.std.PowerLaw()
        # compare python integrator with Gyoto integrator
        self.sp.Parameters = (2., 0.)
        stdsp.Exponent = 0.
        stdsp.Constant = 2.
        self.assertAlmostEqual(self.sp.integrate(1., 2.)
                               /stdsp.integrate(1., 2.), 1., places=7)
        self.sp.Parameters = (2., 1.)
        stdsp.Exponent = 1.
        stdsp.Constant = 2.
        self.assertAlmostEqual(self.sp.integrate(1., 2.)
                               /stdsp.integrate(1., 2.), 1., places=7)
        self.sp.Parameters = (5., -1.)
        stdsp.Exponent = -1.
        stdsp.Constant = 5.
        self.assertAlmostEqual(self.sp.integrate(1., 2.)
                               /stdsp.integrate(1., 2.), 1., places=2)


class TestPythonMetric(unittest.TestCase):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.gg = gyoto.python.PythonMetric()
        self.gg.Module = "gyoto_sample_metrics"

    def pos(self, metric):
        if metric.coordKind() is gyoto.core.GYOTO_COORDKIND_SPHERICAL:
            pos=(10, 6., numpy.pi/4, numpy.pi/3)
        else:
            pos=(10, 6, 2, 4.)
        return pos

    def test_Minkowski(self):
        self.gg.Class = "Minkowski"
        # check christoffel consistency in Cartesian and spherical coordinates
        for self.gg.Spherical in True, False:
            try:
                gyoto.metric.check_christoffel(self.gg,
                                               poslist=(self.pos(self.gg),),
                                               epsilon=1e-7)
            except AssertionError as e:
                self.fail(e.__str__())

    def test_KerrBL(self):
        self.gg.Class = "KerrBL"
        stdgg=gyoto.std.KerrBL()
        # set spin
        self.gg.Parameters = (0.5,)
        stdgg.Spin = 0.5
        # test christoffel consistency
        try:
            gyoto.metric.check_christoffel(self.gg,
                                           poslist=(self.pos(self.gg),),
                                           epsilon=1e-7)
        except AssertionError as e:
            self.fail(e.__str__())
        # Test getRms
        self.assertEqual(self.gg.getRms(), stdgg.getRms(), "Rms is different")
        # Test getRmb
        self.assertEqual(self.gg.getRmb(), stdgg.getRmb(), "Rmb is different")
        # Test getSpecificAngularMomentum
        self.assertEqual(self.gg.getSpecificAngularMomentum(10.),
                         stdgg.getSpecificAngularMomentum(10.),
                         "Specific angular momentum is different")
        # Test gmunu
        pos=[0, 10, 1, 1.]
        self.assertEqual(self.gg.gmunu(pos).tolist(),
                         stdgg.gmunu(pos).tolist(),
                         "gmunu is different")
        # Test christoffel
        pos=[0, 10, 1, 1.]
        self.assertEqual(self.gg.christoffel(pos).tolist(),
                         stdgg.christoffel(pos).tolist(),
                         "christoffel is different")
        # Test getPotential
        self.assertEqual(self.gg.getPotential(pos, 10.),
                         stdgg.getPotential(pos, 10.),
                         "getPotential is different")
        # Test isStopCondition
        coord=[0, 10, 1, 1., 0., 0., 0., 0]
        self.assertEqual(self.gg.isStopCondition(coord),
                         stdgg.isStopCondition(coord),
                         "isStopCondition is different")
        # Test circularVelocity
        self.gg.Keplerian = False
        stdgg.Keplerian = False
        self.assertEqual(self.gg.circularVelocity(pos).tolist(),
                         stdgg.circularVelocity(pos).tolist(),
                         "circularVelocity is different")
        vp=numpy.ndarray(4, dtype=float)
        vc=numpy.ndarray(4, dtype=float)
        self.gg.circularVelocity(pos, vp)
        stdgg . circularVelocity(pos, vc)
        self.assertEqual(vp.tolist(), vc.tolist(),
                         "circularVelocity is different")
        self.gg.circularVelocity(pos, vp, -1)
        stdgg . circularVelocity(pos, vc, -1)
        self.assertEqual(vp.tolist(), vc.tolist(),
                         "circularVelocity is different")
