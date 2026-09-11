#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Unit test file for the BalasinGrumiller metric.

Three kinds of tests are performed:

1. Parameter setting and retrieval.

2. The metric components and all 64 Christoffel symbols are compared, at
   three chosen points, against independently computed values.

3. A consistency check over a wider set of points, comparing the Christoffel
   symbols analytically evaluated using christoffel() against the values
   obtained numerically from gmunu() alone, by inversion and finite
   differences.
"""

import unittest
import numpy
import gyoto.core
import gyoto.std
import gyoto.metric


class TestBalasinGrumiller(unittest.TestCase):
    # Parameters used throughout.
    #
    V0 = 0.000733333
    R = 100.0
    r0 = 1.0
    
    def setUp(self):
        """Set up a BalasinGrumiller metric instance before each test"""
        self.metric = gyoto.core.Metric("BalasinGrumiller")
        self.metric.set("V0", self.V0)
        self.metric.set("R", self.R)
        self.metric.set("r0", self.r0)

    # ---------------------------------------------------------------
    # helpers
    # ---------------------------------------------------------------

    @staticmethod
    def _lookup_g(table, mu, nu):
        return table.get((min(mu, nu), max(mu, nu)), 0.0)

    @staticmethod
    def _lookup_chr(table, a, mu, nu):
        return table.get((a, min(mu, nu), max(mu, nu)), 0.0)

    def assertClose(self, actual, expected, label, rtol=1e-9, atol=1e-15):
        """Compare to a fixed relative precision.

        A relative tolerance keeps the precision uniform across the very
        different magnitudes involved (from 1e-9 to 2500), which an absolute
        criterion such as assertAlmostEqual(places=N) cannot do.  The additive
        atol covers the components that should vanish exactly, where a
        relative tolerance would demand bit-exact equality; round-off in the
        trigonometric cancellations leaves residues of order 1e-16 there.
        """
        tol = atol + rtol * abs(expected)
        self.assertLessEqual(
            abs(actual - expected), tol,
            "{}: got {!r}, expected {!r} (difference {:.3e}, tolerance {:.3e})"
            .format(label, actual, expected, abs(actual - expected), tol))

    # ---------------------------------------------------------------
    # Test 1 - parameter setting and retrieval
    # ---------------------------------------------------------------

    def test_parameters(self):
        """Test parameter getting/setting"""
        self.assertClose(self.metric.get("V0"), self.V0, "V0")
        self.assertClose(self.metric.get("R"), self.R, "R")
        self.assertClose(self.metric.get("r0"), self.r0, "r0")

    # ---------------------------------------------------------------
    # Test 2 - comparison of metric components and Christoffel symbols
    #          against independently computed reference values
    # ---------------------------------------------------------------
    #
    # Only independent components are listed.  The lookup helpers apply the
    # symmetry of the lower indices, and treat any component absent from a
    # table as zero -- so the loops below also assert that the components
    # which should vanish do vanish.

    # 2.1 Equatorial point
    EQUATORIAL = {
        'pos': [0., 50., numpy.pi / 2, 0.],
        'gmunu': {
            (0, 0): -1.0,
            (0, 3): 0.02728482768962097,
            (1, 1): 1.0,
            (2, 2): 2500.0,
            (3, 3): 2499.9992555381777,
        },
        'christoffel': {
            (0, 0, 1): 2.211325543338302e-9,
            (0, 1, 3): 0.000343081548578509,
            (1, 0, 3): -0.00020261494487827393,
            (1, 2, 2): -50.0,
            (1, 3, 3): -49.99998894337228,
            (2, 1, 2): 1.0 / 50.0,
            (3, 0, 1): 8.104597795130958e-8,
            (3, 1, 3): 0.019999997788674458,
        },
    }

    # 2.2 Off-equatorial point
    OFF_EQUATORIAL = {
        'pos': [0., 50., numpy.pi / 4, 0.],
        'gmunu': {
            (0, 0): -1.0,
            (0, 3): 0.030958415625370057,
            (1, 1): 1.0,
            (2, 2): 2500.0,
            (3, 3): 1249.999041576502,
        },
        'christoffel': {
            (0, 0, 1): 6.43865063749218e-9,
            (0, 0, 2): -1.0303305100929405e-7,
            (0, 1, 3): 0.0003591963693300576,
            (0, 2, 3): 0.035118557864864106,
            (1, 0, 3): -0.0002599717438469211,
            (1, 2, 2): -50.0,
            (1, 3, 3): -24.999983903373405,
            (2, 0, 3): 1.6640556199016153e-6,
            (2, 1, 2): 1.0 / 50.0,
            (2, 3, 3): -0.500000103033051,
            (3, 0, 1): 2.079773950775369e-7,
            (3, 0, 2): -3.3281112398032305e-6,
            (3, 1, 3): 0.019999993561349363,
            (3, 2, 3): 1.000000103033051,
        },
    }

    # 2.3 Third point (r=150, \theta=2\pi/3)
    THIRD = {
        'pos': [0., 150., 2 * numpy.pi / 3, 0.],
        'gmunu': {
            (0, 0): -1.0,
            (0, 3): 0.05418322070858099,
            (1, 1): 1.0,
            (2, 2): 22500.0,
            (3, 3): 16874.997064178595,
        },
        'christoffel': {
            (0, 0, 1): 1.9212278381854763e-10,
            (0, 0, 2): 2.271528330815674e-8,
            (0, 1, 3): 0.0003013861109490975,
            (0, 2, 3): -0.038357220641105715,
            (1, 0, 3): -0.00005983535003161125,
            (1, 2, 2): -150.0,
            (1, 3, 3): -112.49999351585605,
            (2, 0, 3): -3.1442321549592736e-7,
            (2, 1, 2): 1.0 / 150.0,
            (2, 3, 3): 0.4330127359651443,
            (3, 0, 1): 3.545798520391778e-9,
            (3, 0, 2): 4.1923095399456977e-7,
            (3, 1, 3): 0.006666666474543883,
            (3, 2, 3): -0.5773502919049092,
        },
    }

    def _check_gmunu(self, case, name):
        if not case['gmunu']:
            self.skipTest(
                "reference values for the {} point not yet filled in".format(name))
        pos = case['pos']
        for mu in range(4):
            for nu in range(4):
                with self.subTest(mu=mu, nu=nu):
                    expected = self._lookup_g(case['gmunu'], mu, nu)
                    self.assertClose(self.metric.gmunu(pos, mu, nu), expected,
                                     "{} g[{}][{}]".format(name, mu, nu))

    def _check_christoffel(self, case, name):
        if not case['christoffel']:
            self.skipTest(
                "reference values for the {} point not yet filled in".format(name))
        pos = case['pos']
        for a in range(4):
            for mu in range(4):
                for nu in range(4):
                    with self.subTest(alpha=a, mu=mu, nu=nu):
                        expected = self._lookup_chr(case['christoffel'], a, mu, nu)
                        self.assertClose(
                            self.metric.christoffel(pos, a, mu, nu), expected,
                            "{} Gamma^{}_[{}][{}]".format(name, a, mu, nu))

    def test_gmunu_equatorial(self):
        """Metric components at r=50, theta=pi/2, against Mathematica."""
        self._check_gmunu(self.EQUATORIAL, "equatorial")

    def test_gmunu_off_equatorial(self):
        """Metric components at r=50, theta=pi/4, against Mathematica."""
        self._check_gmunu(self.OFF_EQUATORIAL, "off-equatorial")

    def test_gmunu_third(self):
        """Metric components at r=150, theta=2pi/3, against Mathematica."""
        self._check_gmunu(self.THIRD, "third")

    def test_christoffel_equatorial(self):
        """All 64 Christoffel components at theta=pi/2, against Mathematica."""
        self._check_christoffel(self.EQUATORIAL, "equatorial")

    def test_christoffel_off_equatorial(self):
        """All 64 Christoffel components at theta=pi/4, against Mathematica.

        This is the test that exercises the theta-derivative of W, which the
        equatorial test cannot: every such term vanishes at theta = pi/2.
        """
        self._check_christoffel(self.OFF_EQUATORIAL, "off-equatorial")

    def test_christoffel_third(self):
        """All 64 Christoffel components at theta=2pi/3, against Mathematica.

        Outside the rods and in the opposite hemisphere, so that cos(theta)
        changes sign relative to the other two reference points.
        """
        self._check_christoffel(self.THIRD, "third")

    def test_christoffel_symmetry(self):
        """Gamma^a_{mn} is symmetric in its lower indices."""
        for case, name in ((self.EQUATORIAL, "equatorial"),
                           (self.OFF_EQUATORIAL, "off-equatorial"),
                           (self.THIRD, "third")):
            pos = case['pos']
            for a in range(4):
                for mu in range(4):
                    for nu in range(mu + 1, 4):
                        with self.subTest(point=name, alpha=a, mu=mu, nu=nu):
                            self.assertClose(
                                self.metric.christoffel(pos, a, mu, nu),
                                self.metric.christoffel(pos, a, nu, mu),
                                "{} symmetry of Gamma^{}_[{}][{}]"
                                .format(name, a, mu, nu))

    # ---------------------------------------------------------------
    # Test 3 - Consistency of analytical and numerical values
    # ---------------------------------------------------------------
    #
    # These are a free choice. Some examples are given below
    CONSISTENCY_POSITIONS = [
        [0., 3.0,   numpy.pi / 4,     0.5],
        [0., 3.0,   1.2,              0.5],
        [0., 50.0,  numpy.pi / 4,     0.5],
        [0., 50.0,  numpy.pi / 3,     2.0],
        [0., 50.0,  numpy.pi / 2,     0.5],
        [0., 99.0,  0.6,              1.0],
        [0., 150.0, numpy.pi / 4,     0.5],
        [0., 150.0, 2 * numpy.pi / 3, 1.0],
    ]

    def test_christoffel_consistency(self):
        """christoffel() agrees with the connection derived from gmunu().

        Compares the given analytical expressions for the Christoffel symbols 
				against a finite-difference estimate obtained by inverting gmunu() 
				and differentiating it.
        """
        gyoto.metric.check_christoffel(self.metric,
                                       poslist=self.CONSISTENCY_POSITIONS)


if __name__ == '__main__':
    unittest.main()
