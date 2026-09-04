import numpy
import unittest
import gyoto.core
import gyoto.std
import inspect
from . import helpers

try:
    import gyoto.lorene
    class TestNeutronStar(unittest.TestCase):

        def test_NeutronStar(self):
            ao=gyoto.core.Astrobj("NeutronStar")
            ao=gyoto.lorene.NeutronStar()
            gg=gyoto.std.KerrBL()
            with self.assertRaises(gyoto.Error):
                ao.metric(gg)
            self.assertIsNone(ao.metric())
            gg=gyoto.lorene.NumericalMetricLorene()
            ao.metric(gg)
            self.assertIsNotNone(ao.metric())
            ao.metric(None)
            self.assertIsNone(ao.metric())
            
        def test_NeutronStarAnalyticEmission(self):
            ao=gyoto.core.Astrobj("NeutronStarAnalyticEmission")
            ao=gyoto.lorene.NeutronStarAnalyticEmission()
            gg=gyoto.std.KerrBL()
            with self.assertRaises(gyoto.Error):
                ao.metric(gg)
            self.assertIsNone(ao.metric())
            gg=gyoto.lorene.NumericalMetricLorene()
            ao.metric(gg)
            self.assertIsNotNone(ao.metric())
            ao.metric(None)
            self.assertIsNone(ao.metric())
            self.assertIsNone(ao.spectrum())
            sp=gyoto.std.BlackBody()
            ao.spectrum(sp)
            self.assertIsNotNone(ao.spectrum())
            
        def test_NeutronStarModelAtmosphere(self):
            ao=gyoto.core.Astrobj("NeutronStarModelAtmosphere")
            ao=gyoto.lorene.NeutronStarModelAtmosphere()
            gg=gyoto.std.KerrBL()
            with self.assertRaises(gyoto.Error):
                ao.metric(gg)
            self.assertIsNone(ao.metric())
            gg=gyoto.lorene.NumericalMetricLorene()
            ao.metric(gg)
            self.assertIsNotNone(ao.metric())
            ao.metric(None)
            self.assertIsNone(ao.metric())

    class TestRotStar3_1(unittest.TestCase):
        def test_christoffel(self):
            metric=gyoto.lorene.RotStar3_1()
            try:
                metric.file('../bin/.check-lorene/resu.d')
            except gyoto.core.Error as e:
                self.skipTest('RotStar3_1::christoffel (metric needs to be precomputed)')
            try:
                gyoto.metric.check_christoffel(metric)
            except AssertionError as e:
                self.fail(e.__str__())

    class TestLoreneXMLio(unittest.TestCase):

        default_verbosity=gyoto.core.verbose()

        def setUp(self):
            gyoto.core.verbose(0)

        def tearDown(self):
            gyoto.core.verbose(self.default_verbosity)

        def invalid(self, classname, cls):
            return (not inspect.isclass(cls)
                    or not issubclass(cls, gyoto.core.Object))

        def test_xmlio(self):
            helpers.test_xmlio(self, gyoto.lorene, self.invalid)

except ImportError:            
    import warnings
    warnings.warn('Could not load plug-in "lorene"')

