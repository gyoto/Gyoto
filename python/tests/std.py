import numpy
import unittest
import gyoto.core
import gyoto.std
import inspect

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

class TestHayward(unittest.TestCase):

    def test_mass(self):
        gg=gyoto.std.Hayward()
        m=2.
        gg.mass(m)
        self.assertTrue((gg.mass() == m))

    def test_charge(self):
        gg=gyoto.std.Hayward()
        b=0.5
        gg.charge(b)
        self.assertTrue((gg.charge() == b))

class TestStar(unittest.TestCase):

    def test_setInitCoord(self):
        st=gyoto.std.Star()
        gg=gyoto.std.KerrBL()
        st.metric(gg)
        pos_list=(600., 9., 1.5707999999999999741, 0)
        vel_list=(0., 0., 0.037037)
        st.setInitCoord(pos_list, vel_list)
        dst=gyoto.core.vector_double()
        dst2=gyoto.core.vector_double()
        st.getInitialCoord(dst)
        st.setPosition(pos_list)
        st.setVelocity(vel_list)
        st.getInitialCoord(dst2)
        self.assertTrue((numpy.asarray(dst) == numpy.asarray(dst2)).all())

class TestMinkowski(unittest.TestCase):

    def _compute_r_norm(self, met, st, pos, v, tmax=1e6):
        t, x, y, z, tdot, xdot, ydot, zdot=self._compute_orbit(met, st, pos, v, tmax)
        if met.spherical():
            r=x
        else:
            r=numpy.sqrt(x**2+y**2+z**2)
        n=t.size
        norm=numpy.asarray([met.norm([t[i], x[i], y[i], z[i]],
                                     [tdot[i], xdot[i], ydot[i], zdot[i]])
                            for  i in range(n)])
        return r, norm

    def _compute_orbit(self, met, st, pos, v, tmax=1e6):
        st.initCoord(numpy.append(pos, v))
        st.xFill(tmax)
        n=st.get_nelements()
        t=numpy.ndarray(n)
        x=numpy.ndarray(n)
        y=numpy.ndarray(n)
        z=numpy.ndarray(n)
        tdot=numpy.ndarray(n)
        xdot=numpy.ndarray(n)
        ydot=numpy.ndarray(n)
        zdot=numpy.ndarray(n)
        st.get_t(t)
        st.getCoord(t, x, y, z, tdot, xdot, ydot, zdot)
        return t, x, y, z, tdot, xdot, ydot, zdot

    def test_Keplerian(self):
        met=gyoto.std.Minkowski()
        met.keplerian(True)
        st=gyoto.std.Star()
        st.metric(met)

        # Cartesian coordinates
        # particlae is at x=1000:
        pos=[0., 1000., 0., 0.]
        v=met.circularVelocity(pos)
        r, norm=self._compute_r_norm(met, st, pos, v)
        self.assertLess( numpy.abs(r-1000.).max(), 1e-6)
        self.assertLess( numpy.abs(norm+1.).max(), 1e-6 )

        # same velocity at z=1000:
        pos=[0., 0., 0., 1000.]
        r, norm=self._compute_r_norm(met, st, pos, v)
        self.assertLess( numpy.abs(r-1000.).max(), 1e-6)
        self.assertLess( numpy.abs(norm+1.).max(), 1e-6 )

        # same velocity at x=z=1000/sqrt(2)
        pos=[0., 1000./numpy.sqrt(2.), 0., 1000./numpy.sqrt(2.)]
        r, norm=self._compute_r_norm(met, st, pos, v)
        self.assertLess( numpy.abs(r-1000.).max(), 1e-6)
        self.assertLess( numpy.abs(norm+1.).max(), 1e-6 )

        # a elliptical orbit
        pos=[0., 1000., 0., 0.]
        v3=[0, 0.015, 0]
        tdot=met.SysPrimeToTdot(pos, v3)
        v=[tdot, v3[0]*tdot, v3[1]*tdot, v3[2]/tdot]
        r, norm=self._compute_r_norm(met, st, pos, v)
        self.assertLess( numpy.abs(norm+1.).max(), 1e-6 )

        # circular orbit starting at x=3
        pos=[0., 3., 0., 0.]
        v=met.circularVelocity(pos)
        st.deltaMaxOverR(0.1)
        r, norm=self._compute_r_norm(met, st, pos, v, tmax=50.)
        self.assertLess( numpy.abs(r-3.).max(), 1e-6)
        self.assertLess( numpy.abs(norm+1.).max(), 1e-6 )

        # Spherical coordinates
        st=gyoto.std.Star()
        st.metric(met)
        met.spherical(True)

        # again starting at x=1000
        pos=[0., 1000., numpy.pi*0.5, 0.]
        v=met.circularVelocity(pos)
        r, norm=self._compute_r_norm(met, st, pos, v)
        self.assertLess( numpy.abs(r-1000.).max(), 1e-6)
        self.assertLess( numpy.abs(norm+1.).max(), 1e-6 )

        # elliptical orbit
        v[3] = v[3]/2.
        r, norm=self._compute_r_norm(met, st, pos, v)
        self.assertLess( numpy.abs(norm+1.).max(), 1e-3 )

        # circular orbit starting at x=z=1000/sqrt(2)
        pos=[0., 1000., numpy.pi/2., 0.]
        v=met.circularVelocity(pos)
        pos[2]=numpy.pi/4.
        v[3]=v[3]/numpy.sin(numpy.pi/4.)
        r, norm=self._compute_r_norm(met, st, pos, v)
        self.assertLess( numpy.abs(r-1000.).max(), 1e-6)
        self.assertLess( numpy.abs(norm+1.).max(), 1e-6 )

        # circular orbit starting at x=3
        pos=[0., 3., numpy.pi*0.5, 0.]
        st.deltaMaxOverR(0.1)
        st.initCoord(numpy.append(pos, v))
        v=met.circularVelocity(pos)
        r, norm=self._compute_r_norm(met, st, pos, v, tmax=50.)
        self.assertLess( numpy.abs(r-3.).max(), 1e-6)
        self.assertLess( numpy.abs(norm+1.).max(), 1e-6 )

class TestStdMetric(unittest.TestCase):
    def test_christoffel(self):
        nspace=gyoto.std
        for classname, cls in inspect.getmembers(nspace):
            if (not inspect.isclass(cls)
                or not issubclass(cls, gyoto.core.Metric)):
                continue
            try:
                gyoto.metric.check_christoffel(cls)
            except AssertionError as e:
                self.fail(e.__str__())
