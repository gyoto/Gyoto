import numpy
import unittest
import gyoto.core
import gyoto.std

gyoto.core.requirePlugin("lorene", 1)

if gyoto.core.havePlugin("lorene"):
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
else:
    import warnings
    warnings.warn('Could not load plug-in "lorene"')

