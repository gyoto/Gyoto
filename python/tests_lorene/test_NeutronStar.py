import numpy
import unittest
import gyoto
import gyoto_std
import gyoto_lorene

class TestNeutronStar(unittest.TestCase):

    def test_NeutronStar(self):
        ao=gyoto.Astrobj("NeutronStar")
        ao=gyoto_lorene.NeutronStar()
        gg=gyoto_std.KerrBL()
        with self.assertRaises(gyoto.Error):
            ao.metric(gg)
        self.assertIsNone(ao.metric())
        gg=gyoto_lorene.NumericalMetricLorene()
        ao.metric(gg)
        self.assertIsNotNone(ao.metric())
        ao.metric(None)
        self.assertIsNone(ao.metric())

