import numpy
import unittest
import gyoto.core
import gyoto.std
import gyoto.metric
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

    def test_gmunu(self):
        metric=gyoto.std.Hayward()
        I=numpy.zeros((4,4))
        for i in range(4):
            I[i,i]=1.
        for r in (-2, 0.5, 6):
            pos=(10, r, numpy.pi/4, numpy.pi/3)
            g=metric.gmunu(pos)
            gup=metric.gmunu_up(pos)
            ggup=numpy.linalg.multi_dot((g, gup))
            for mu in range(4):
                for nu in range(4):
                    self.assertAlmostEqual(metric.gmunu(pos, mu, nu), g[mu, nu])
                    self.assertAlmostEqual(metric.gmunu_up(pos, mu, nu), gup[mu, nu], 7,
                                           "mu: {}, nu: {}".format(mu, nu))
                    self.assertAlmostEqual(ggup[mu, nu],I[mu,nu])


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
    def pos(self, metric):
        if metric.coordKind() is gyoto.core.GYOTO_COORDKIND_SPHERICAL:
            pos=(10, 6., numpy.pi/4, numpy.pi/3)
        else:
            pos=(10, 6, 2, 4.)
        return pos

    def metric(self, cls):
        metric=cls()
        # All Kerr-like: use non-zero spin
        try:
            metric.spin(0.5)
        except AttributeError:
            pass
        # Chern-Simons: non-zero dzeta
        try:
            metric.dzetaCS(0.5)
        except AttributeError:
            pass
        # Complex: add two metrics
        try:
            metric.append(gyoto.std.KerrBL())
            metric.append(gyoto.std.KerrBL())
        except AttributeError:
            pass
        # Shift: add a submetric
        if (cls == gyoto.std.Shift):
            metric.subMetric(gyoto.std.KerrKS())
            metric.offset((1., 1., 1., 1.))
        return metric

    def invalid(self, classname, cls):
        return (not inspect.isclass(cls)
                or not issubclass(cls, gyoto.core.Metric))

    def test_christoffel(self):
        nspace=gyoto.std
        for classname, cls in inspect.getmembers(nspace):
            if (self.invalid(classname, cls)):
                continue
            metric=self.metric(cls)
            try:
                gyoto.metric.check_christoffel(metric, poslist=(self.pos(metric),), epsilon=1e-7)
            except AssertionError as e:
                self.fail(e.__str__())
            pos=self.pos(metric)
            G=metric.christoffel(pos)
            G2=numpy.ones((4,4,4))
            retval=metric.christoffel(G2, pos)
            self.assertEqual(retval, 0, 'christoffel errors out')
            for a in range(4):
                for mu in range(4):
                    for nu in range(4):
                        self.assertAlmostEqual(metric.christoffel(pos, a, mu, nu), G[a, mu, nu], 7, classname)
                        self.assertAlmostEqual(metric.christoffel(pos, a, mu, nu), G2[a, mu, nu], 7, classname)

    def test_jacobian(self):
        nspace=gyoto.std
        for classname, cls in inspect.getmembers(nspace):
            if (self.invalid(classname, cls)):
                continue
            metric=self.metric(cls)
            pos=self.pos(metric)
            jac=metric.jacobian(pos)
            jacn=gyoto.metric.jacobian_numerical(metric, pos, epsilon=1e-7)
            for a in range(4):
                for m in range(4):
                    for n in range(4):
                        self.assertAlmostEqual(jac[a,m,n],
                                               jacn[a,m,n], 7, classname)

    def test_gmunu(self):
        nspace=gyoto.std
        for classname, cls in inspect.getmembers(nspace):
            if (self.invalid(classname, cls)):
                continue
            metric=self.metric(cls)
            pos=self.pos(metric)
            g=metric.gmunu(pos)
            for mu in range(4):
                for nu in range(4):
                    self.assertAlmostEqual(metric.gmunu(pos, mu, nu), g[mu, nu], 7, classname)

    def test_gmunu_up(self):
        nspace=gyoto.std
        for classname, cls in inspect.getmembers(nspace):
            if (self.invalid(classname, cls)):
                continue
            metric=self.metric(cls)
            pos=self.pos(metric)
            gup=metric.gmunu_up(pos)
            gup2, jac=metric.gmunu_up_and_jacobian(pos)
            g=metric.gmunu(pos)
            gup3=numpy.linalg.inv(g)
            for mu in range(4):
                for nu in range(4):
                    self.assertAlmostEqual(metric.gmunu_up(pos, mu, nu), gup[mu, nu], 7,
                                           "class: {}, mu: {}, nu: {} ({})".format(classname, mu, nu, 'coef/matrix'))
                    self.assertAlmostEqual(gup3[mu, nu], gup2[mu, nu], 7,
                                           "class: {}, mu: {}, nu: {} ({})".format(classname, mu, nu, 'inv(g)/*_and_jacobian'))
                    self.assertAlmostEqual(gup2[mu, nu], gup[mu, nu], 7,
                                           "class: {}, mu: {}, nu: {} ({})".format(classname, mu, nu, 'matrix/*_and_jacobian'))

    def test_gmunu_gmunu_up(self):
        nspace=gyoto.std
        I=numpy.zeros((4,4))
        for i in range(4):
            I[i,i]=1.
        for classname, cls in inspect.getmembers(nspace):
            if (self.invalid(classname, cls)):
                continue
            metric=self.metric(cls)
            pos=self.pos(metric)
            g=metric.gmunu(pos)
            gup=metric.gmunu_up(pos)
            self.assertAlmostEqual(numpy.abs(numpy.linalg.multi_dot((g, gup))-I).max(), 0.)
