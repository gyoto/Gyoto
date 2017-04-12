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
            
else:
    import warnings
    warnings.warn('Could not load plug-in "lorene"')

