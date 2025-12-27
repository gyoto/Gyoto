import unittest, numpy, sys
import gyoto.python

class IntrospectingThinDisk(gyoto.python.ThinDiskBase):
    properties = {"ThisPointer": "long",
                  "SelfPointer": "long",
                  "MyDouble": "double",
                  "MyLong": "long",
                  "MyUnsignedLong": "unsigned_long",
                  "MySize_t": "size_t",
                  "MyBool": "bool",
                  "MyString": "string",
                  "MyFilename": "filename",
                  "MyVectorDouble": "vector_double",
                  "MyVectorUnsignedLong": "vector_unsigned_long",
                  "MyMetric": "metric",
                  "MyScreen": "screen",
                  "MyAstrobj": "astrobj",
                  "MySpectrum": "spectrum",
                  "MySpectrometer": "spectrometer",
                  "MyEmpty": "empty",
                  }
    def set(self, key, *args):
        if key in ("ThisPointer", "SelfPointer"):
            raise ValueError(f"read-only property: {key}")
        else: super().set(key, *args)
    def get(self, key, *args):
        if key == "ThisPointer":
            return self.this.getPointer()
        if key == "SelfPointer":
            return self.getPointer()
        return super().get(key, *args)

class TestIntrospectingThinDisk(unittest.TestCase):
    '''Test a ThinDisk implemented in Python
    '''
    def test_pointers(self):
        '''Test consistency of pointers of classes implemented in Python
        '''
        td = IntrospectingThinDisk()
        self.assertEqual(td.ThisPointer, td.SelfPointer,
                         "instance.ThisPointer==instance.SelfPointer with self instance")
        self.assertEqual(td.ThisPointer, td.getPointer(),
                         "instance.ThisPointer==instance.getPointer() with self instance")
        self.assertEqual(gyoto.core.gyotoid(td), td.Instance,
                         "gyotoid(instance)==instance.Instance with self instance")
        ao=gyoto.python.PythonThinDisk()
        ao.Module="python"
        ao.Class="IntrospectingThinDisk"
        self.assertEqual(ao.ThisPointer, ao.SelfPointer,
                         "thindisk.ThisPointer==thindisk.SelfPointer with internal instance")
        self.assertEqual(ao.ThisPointer, ao.getPointer(),
                         "thindisk.ThisPointer==thindisk.getPointer() with internal instance")
        ao=gyoto.python.PythonThinDisk()
        ao.Instance=gyoto.core.gyotoid(td)
        self.assertEqual(ao.ThisPointer, ao.SelfPointer,
                         "thindisk.ThisPointer==thindisk.SelfPointer with external instance")
        self.assertEqual(ao.ThisPointer, ao.getPointer(),
                         "thindisk.ThisPointer==thindisk.getPointer() with external instance")
        self.assertEqual(ao.InnerRadius, td.InnerRadius,
                         "thindisk.InnerRadius==instance.InnerRadius for default value")
        gen=gyoto.astrobj.Generic(ao)
        default=gen.InnerRadius
        gen.InnerRadius=default+1
        self.assertEqual(ao.InnerRadius, td.InnerRadius,
                         "thindisk.InnerRadius==instance.InnerRadius with external instance")
        self.assertEqual(ao.InnerRadius, gen.InnerRadius,
                         "thindisk.InnerRadius==generic.InnerRadius with external instance")
        self.assertEqual(ao.InnerRadius, default+1,
                         "thindisk.InnerRadius==<set value> with external instance")

    def test_properties(self):
        '''Test setting and getting properties of all types in a PythonThinDisk
        '''
        td=IntrospectingThinDisk()
        gen=gyoto.core.Astrobj(td)

        # scalar properties
        for key, val in {"MyDouble": 2.0,
                         "MyLong": -1023,
                         "MyUnsignedLong": 1023,
                         "MySize_t": 123456,
                         "MyBool": True,
                         "MyBool": False,
                         "MyString": "Some , interresting: text ?",
                         "MyFilename": "zqlieksdh filserkjd glskerjh"
                         }.items():
            try:
                setattr(td, key, val)
            except:
                self.fail(f"failed setting thindisk.{key} to {val}")
            self.assertEqual(getattr(gen, key), val,
                             f"retrieved value for generic.{key} is not the same as was set")
            try:
                setattr(gen, key, val)
            except:
                self.fail(f"failed setting generic.{key} to {val}")
            self.assertEqual(getattr(td, key), val,
                             f"retrieved value for thindisk.{key} is not the same as was set")

        # vector properties
        for key, val in {"MyVectorDouble": (1.0, 2.0, 3.4),
                         "MyVectorUnsignedLong": (1, 0, 2, 3)
                         }.items():
            try:
                setattr(td, key, val)
            except:
                self.fail(f"failed setting thindisk.{key} to {val}")
            try:
                out=getattr(gen, key)
            except:
                self.fail(f"failed getting generic.{key}")
            self.assertEqual(len(out), len(val),
                             f"retrieved value for generic.{key} is not the same length as was set")
            for i in range(len(out)):
                self.assertEqual(out[i], val[i],
                                 f"retrieved value for generic{key}[{i}] is not the same length as was set")
            try:
                setattr(gen, key, val)
            except:
                self.fail(f"failed setting generic.{key} to {val}")
            try:
                out=getattr(td, key)
            except:
                self.fail(f"failed getting thindisk.{key}")
            self.assertEqual(len(out), len(val),
                             f"retrieved value for thindisk.{key} is not the same length as was set")
            for i in range(len(out)):
                self.assertEqual(out[i], val[i],
                                 f"retrieved value for thindisk{key}[{i}] is not the same length as was set")

        # object properties
        for key, val in {"MyMetric": gyoto.metric.KerrBL(),
                         "MyScreen": gyoto.core.Screen(),
                         "MyAstrobj": gyoto.astrobj.Star(),
                         "MySpectrum": gyoto.spectrum.PowerLaw(),
                         "MySpectrometer": gyoto.spectrometer.Uniform()
                         }.items():
            try:
                setattr(td, key, val)
            except:
                self.fail(f"failed setting thindisk.{key} to {val}")
            try:
                out=getattr(gen, key)
            except:
                self.fail(f"failed getting generic.{key}")
            self.assertEqual(out.getPointer(), val.getPointer(),
                             f"retrieved value for generic.{key} is not the same as was set")
            try:
                setattr(gen, key, val)
            except:
                self.fail(f"failed setting generic.{key} to {val}")
            try:
                out=getattr(td, key)
            except:
                self.fail(f"failed getting thindisk.{key}")
            self.assertEqual(out.getPointer(), val.getPointer(),
                             f"retrieved value for thindisk.{key} is not the same as was set")

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
