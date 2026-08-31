import numpy
import unittest
import gyoto.core
import gyoto.std
import inspect

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

except ImportError:            
    import warnings
    warnings.warn('Could not load plug-in "lorene"')

