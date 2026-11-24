"""
Unit tests for Gyoto's Kerr2PN metric.
"""
import unittest

import numpy as np
import gyoto.core, gyoto.std, gyoto.metric

CHI = 0.5                                      # same spin as the Gyoto unit test
DIRECTION = np.array([6., 2., 4.])
DIRECTION /= np.linalg.norm(DIRECTION)         # same direction as the unit-test point (6,2,4)
TEST_POINT = (10., 6., 2., 4.)                 # the unit-test point

NAMES = ["Gamma^0_0i", "Gamma^0_ij", "Gamma^i_00", "Gamma^i_0j", "Gamma^i_jk"]


def make(cls, spin=CHI):
    m = cls()
    m.spin(spin)
    return m


def at_radius(r):
    """Point at coordinate radius r along DIRECTION (t = 10)."""
    return (10., *(r * DIRECTION))


def per_class(G, H=None):
    """Largest |value| (or |difference|) of each family of Christoffel symbols."""
    D = G if H is None else G - H
    D = np.abs(D)
    return [D[0, 0, 1:].max(),      # Gamma^0_0i
            D[0, 1:, 1:].max(),     # Gamma^0_ij
            D[1:, 0, 0].max(),      # Gamma^i_00
            D[1:, 0, 1:].max(),     # Gamma^i_0j
            D[1:, 1:, 1:].max()]    # Gamma^i_jk


def accel(G, v):
    """Coordinate acceleration of a test particle with u = (1, v)."""
    u = np.array([1., *v])
    num = np.einsum('amn,m,n->a', G, u, u)
    return -(num[1:] - num[0] * v)


def log_slope(x0, x1, y0, y1):
    """Exponent p if y goes like x^p between two points."""
    return np.log(y1 / y0) / np.log(x1 / x0)


class TestKerr2PN(unittest.TestCase):
    """Sanity checks for gyoto.std.Kerr2PN (Kerr truncated at 2PN order)."""

    @classmethod
    def setUpClass(cls):
        cls.kerr2pn = make(gyoto.std.Kerr2PN)
        cls.kerrharm = make(gyoto.std.KerrHarmonic)   # exact Kerr in harmonic coordinates

    # ------------------------------------------------------------------
    # CHECK 1.  Is gmunu_up really the inverse of gmunu?
    def test_gmunu_up_inverts_gmunu(self):
        fresh = make(gyoto.std.Kerr2PN)
        g = fresh.gmunu(TEST_POINT)
        # Kerr2PN has no gmunu_up of its own, so Gyoto's generic one is used:
        # it inverts gmunu using a pivot algorithm.
        for label in ("1st call", "2nd call"):   # generic gmunu_up
                                                 # used to not behave
                                                 # the same on the 1st
                                                 # and later calls
            with self.subTest(call=label):
                gup = fresh.gmunu_up(TEST_POINT)
                self.assertLess(np.abs(g @ gup - np.eye(4)).max(), 1e-15,
                                "|g.g_up - 1|")
                self.assertLess(np.abs(gup - np.linalg.inv(g)).max(), 1e-15,
                                "|g_up - numpy inv(g)|")
                self.assertLess(np.abs(gup - gup.T).max(), 1e-16,
                                "|g_up - g_up transposed|")
        # For comparison, KerrHarmonic has its own gmunu_up.
        g2 = self.kerrharm.gmunu(TEST_POINT)
        gup2 = self.kerrharm.gmunu_up(TEST_POINT)
        self.assertLess(np.abs(g2 @ gup2 - np.eye(4)).max(), 1e-15)
        self.assertLess(np.abs(gup2 - gup2.T).max(), 1e-16)

    # ------------------------------------------------------------------
    # CHECK 2.  Christoffels written in the code versus Christoffels
    #           estimated numerically from gmunu itself.
    def test_christoffel_matches_numerical(self):
        G = self.kerr2pn.christoffel(TEST_POINT)
        Gn = gyoto.metric.christoffel_numerical(self.kerr2pn, TEST_POINT, 1e-7)
        # The existing Gyoto unit test allows 1e-6 on this component.
        self.assertLess(abs(G[0, 0, 1] - Gn[0, 0, 1]), 4e-4,
                        "Gamma^0_01: code vs numerical at the test point")
        # Away from the test point, |code - numerical| must decrease with r
        # (it is made of higher-order PN terms), until the finite-difference
        # noise (~1e-9) takes over beyond r ~ 60.
        radii = [7.5, 15., 30., 60.]
        rows = [per_class(self.kerr2pn.christoffel(at_radius(r)),
                          gyoto.metric.christoffel_numerical(self.kerr2pn, at_radius(r), 1e-7))
                for r in radii]
        for k, name in enumerate(NAMES):
            with self.subTest(family=name):
                self.assertLess(rows[-1][k], rows[0][k],
                                name + ": difference should decrease with r")

    # ------------------------------------------------------------------
    # CHECK 3.  Kerr2PN versus the exact harmonic-coordinate Kerr.
    def test_gmunu_close_to_kerr_harmonic(self):
        radii = [7.5, 15., 30., 60., 120., 240.]
        gd = [np.abs(self.kerr2pn.gmunu(at_radius(r))
                     - self.kerrharm.gmunu(at_radius(r))).max() for r in radii]
        # The 2PN truncation error must decrease with r; the first missing term
        # is 3PN (or the 2.5PN spin term), i.e. slope p <= -3 for the metric.
        self.assertLess(gd[-1], gd[0])
        p = log_slope(radii[0], radii[-1], gd[0], gd[-1])
        self.assertLess(p, -2., "metric difference should fall faster than r^-2 (p = %.2f)" % p)

    def test_christoffel_close_to_kerr_harmonic(self):
        radii = [7.5, 15., 30., 60., 120., 240.]
        rows = [per_class(self.kerr2pn.christoffel(at_radius(r)),
                          self.kerrharm.christoffel(at_radius(r))) for r in radii]
        for k, name in enumerate(NAMES):
            with self.subTest(family=name):
                self.assertLess(rows[-1][k], rows[0][k],
                                name + ": |Kerr2PN - KerrHarmonic| should decrease with r")

    # ------------------------------------------------------------------
    # CHECK 4.  Does the code give the right motion?  Coordinate
    #           acceleration of a test particle (u = (1, v), v^2 ~ 1/r),
    #           Kerr2PN versus the exact harmonic-coordinate Kerr geodesic,
    #           scaled by the Newtonian 1/r^2.
    def test_acceleration_close_to_kerr_harmonic(self):
        w = np.array([0.35, 0.8, 0.45])
        radii = [20., 40., 80., 160., 320., 640., 1280.]
        for spin in (0., 0.5):
            with self.subTest(spin=spin):
                k2 = make(gyoto.std.Kerr2PN, spin)
                kh = make(gyoto.std.KerrHarmonic, spin)
                err = []
                for r in radii:
                    p, v = at_radius(r), np.sqrt(1. / r) * w
                    err.append(np.linalg.norm(accel(k2.christoffel(p), v)
                                              - accel(kh.christoffel(p), v)) * r**2)
                # The scaled error must decrease with r...
                for i in range(len(radii) - 1):
                    self.assertLess(err[i + 1], err[i])
                # ...with an asymptotic slope p = -3 without spin (the first
                # missing term is 3PN) and p ~ -2.5 with spin (the 2.5PN spin
                # term that the 2PN truncation leaves out on purpose).
                p = log_slope(radii[-2], radii[-1], err[-2], err[-1])
                self.assertLess(p, -2., "asymptotic slope p = %.2f" % p)


if __name__ == "__main__":
    unittest.main()
