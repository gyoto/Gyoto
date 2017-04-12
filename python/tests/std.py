import numpy
import unittest
import gyoto.core
import gyoto.std

class TestChernSimons(unittest.TestCase):

    def test_setInitCoord(self):
        gg=gyoto.std.ChernSimons()
        zeta=0.5
        gg.dzetaCS(zeta)
        self.assertTrue((gg.dzetaCS() == zeta))

class TestDeformedTorus(unittest.TestCase):

    def test_DeformedTorus(self):
        ao=gyoto.core.Astrobj("DeformedTorus")
        ao=gyoto.std.DeformedTorus()
        sp=gyoto.std.BlackBody()
        ao.spectrum(sp)
        r=4.
        ao.largeRadius(r)
        self.assertTrue((ao.largeRadius() == r))

class TestDynamicalDiskBolometric(unittest.TestCase):

    def test_DynamicalDiskBolometric(self):
        ao=gyoto.core.Astrobj("DynamicalDiskBolometric")
        ao=gyoto.std.DynamicalDiskBolometric()
        self.assertTrue(True)

class TestEquatorialHotSpot(unittest.TestCase):

    def test_EquatorialHotSpot(self):
        ao=gyoto.core.Astrobj("EquatorialHotSpot")
        ao=gyoto.std.EquatorialHotSpot()
        r=4.
        ao.spotRadSize(r)
        self.assertTrue((ao.spotRadSize() == r))
        kind="RadialBeaming"
        ao.beaming(kind)
        self.assertTrue((ao.beaming() == kind))
        self.assertRaises(gyoto.core.Error, ao.beaming, "foo")

class TestInflateStar(unittest.TestCase):

    def test_InflateStar(self):
        ao=gyoto.core.Astrobj("InflateStar")
        ao=gyoto.std.InflateStar()
        ao.radius(0.)
        ao.radiusStop(2.)
        ao.timeInflateInit(0.)
        ao.timeInflateStop(2.)
        self.assertTrue(abs(ao.radiusAt(1.)-1.) < 1e-6)

class TestOscilTorus(unittest.TestCase):

    def test_OscilTorus(self):
        ao=gyoto.core.Astrobj("OscilTorus")
        ao=gyoto.std.OscilTorus()
        r=4.
        ao.largeRadius(r)
        self.assertTrue((ao.largeRadius() == r))
        perturbkind="Vertical"
        ao.perturbKind(perturbkind)
        self.assertTrue((ao.perturbKind() == perturbkind))

class TestRezzollaZhidenko(unittest.TestCase):

    def test_mass(self):
        gg=gyoto.std.RezzollaZhidenko()
        m=2.
        gg.mass(m)
        self.assertTrue((gg.mass() == m))

    def test_aparam(self):
        gg=gyoto.std.RezzollaZhidenko()
        a=(1., 2., 3., 4.)
        gg.aparam(a)
        self.assertTrue((gg.aparam() == a))

class TestStar(unittest.TestCase):

    def test_setInitCoord(self):
        st=gyoto.std.Star()
        gg=gyoto.std.KerrBL()
        st.metric(gg)
        pos_list=(600., 9., 1.5707999999999999741, 0)
        vel_list=(0., 0., 0.037037)
        st.setInitCoord(pos_list, vel_list)
        dst=numpy.zeros(8, float)
        dst2=numpy.zeros(8, float)
        st.getInitialCoord(dst)
        st.setPosition(pos_list)
        st.setVelocity(vel_list)
        st.getInitialCoord(dst2)
        self.assertTrue((dst == dst2).all())
