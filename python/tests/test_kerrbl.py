#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright 2011, 2014 Frederic Vincent, Thibaut Paumard
#    Copyright 2026 Julien Brulé, Thibaut Paumard
#
#    This file is part of Gyoto.
#
#    Gyoto is free software: you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    Gyoto is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.
#

import unittest
import gyoto.core
import numpy as np
import gyoto.metric
import gyoto.astrobj
import gyoto.spectrum
import gyoto.spectrometer
import gyoto.std

gyoto.core.requirePlugin('stdplug')

class TestKerrBL(unittest.TestCase):
    def construct_kerrbl(self):
        gg = gyoto.std.KerrBL()
        self.assertIsInstance(gg, gyoto.std.KerrBL)

    def test_MakeCoord_kerrbl(self):
        gg = gyoto.std.KerrBL()
        aa = 0.995
        gg.set("Spin", aa)
        # initial conditions :
        ri = 10.791
        thetai = 1.5708
        phii = 0.
        ti = 0.
        pri = 0.  # canonical momentum
        pthetai = 0.
        cst = [1, 0.921103, 2., 0.]  # 4 Kerr cst of motion mu, E, L, Q

        if cst[3] == 0.:
            cst = np.append(cst, 1.)
        else:
            cst = np.append(cst, 1. / cst[3])

        yinit = [ti, ri, thetai, phii, -cst[1], pri, pthetai, cst[2]]

        coordout = np.zeros(8)
        gg.MakeCoord(yinit, cst, coordout)

        np.testing.assert_array_almost_equal(
            coordout,
            [0, 10.791, 1.5708, 0, 1.12641, 0, 0, 0.0187701]
        )

    def test_setter_getter_kerrbl(self):
        gg = gyoto.std.KerrBL()
        gg.set("Spin", 0.7)
        self.assertEqual(gg.get("Spin"), 0.7)
        gg.set("Spin", 0.9)
        self.assertEqual(gg.get("Spin"), 0.9)

        self.assertNotEqual(gg.getPointer(), 0)

        gg.deltaMin(40)
        self.assertEqual(gg.get("DeltaMin"), 40)

        gg.deltaMax(400)
        self.assertEqual(gg.get("DeltaMax"), 400)

        gg.difftol(1e-3)
        self.assertEqual(gg.get("DiffTol"), 1e-3)

        gg.deltaMaxOverR(2e-3)
        self.assertEqual(gg.get("DeltaMaxOverR"), 2e-3)

        gg2 = gg.clone()
        self.assertEqual(gg2.get("DeltaMin"), 40)
        self.assertEqual(gg2.get("DeltaMax"), 400)
        self.assertEqual(gg2.get("DiffTol"), 1e-3)
        self.assertEqual(gg2.get("DeltaMaxOverR"), 2e-3)

    def test_integrators(self):
        gg = gyoto.std.KerrBL()
        gg.set("Spin", 0.995)
        gg2 = gg.clone()
        gg2.set("GenericIntegrator", True)
        st = gyoto.std.Star()
        st.metric(gg)
        st.setInitCoord([0., 10.791, np.pi / 2., 0],
                        [0., 0., 0.016664])
        st2 = gyoto.std.Star()
        st2.metric(gg2)
        st2.setInitCoord([0., 10.791, np.pi / 2., 0],
                         [0., 0., 0.016664])

        st.integrator("Legacy")
        st2.integrator("Legacy")

        dates = np.arange(100, dtype=np.float64) + 1.
        n = len(dates)

        r = np.ndarray(n)
        theta = np.ndarray(n)
        phi = np.ndarray(n)

        r2 = np.ndarray(n)
        theta2 = np.ndarray(n)
        phi2 = np.ndarray(n)

        # récupérer les 3 autres coordonnées
        st.getCoord(dates, r, theta, phi)
        st2.getCoord(dates, r2, theta2, phi2)

        np.testing.assert_array_almost_equal(
            [r, theta, phi],
            [r2, theta2, phi2],
            decimal=3
        )

        for integrator in ["runge_kutta_cash_karp54",
                           "runge_kutta_fehlberg78",
                           "runge_kutta_dopri5",
                           "runge_kutta_cash_karp54_classic"]:
            st2.Integrator = integrator
            st2.getCoord(dates, r2, theta2, phi2)
            np.testing.assert_array_almost_equal(
                [r, theta],
                [r2, theta2],
                decimal=3
            )
            np.testing.assert_array_almost_equal(
                np.mod(phi - phi2 + np.pi, 2. * np.pi) - np.pi,
                np.zeros(n),
                decimal=3
            )

if __name__ == '__main__':
    unittest.main()
