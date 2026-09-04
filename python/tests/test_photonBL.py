#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Copyright 2011, 2013 Thibaut Paumard
#    Copyright 2026 Julien Brulé & Thibaut Paumard
#
#    This file is part of Gyoto.
#
#    Gyoto is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
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
# import os
import unittest
import gyoto.core
import numpy as np
import gyoto.metric
import gyoto.astrobj
import gyoto.spectrum
import gyoto.spectrometer
import matplotlib.pyplot as plt

gyoto.core.requirePlugin('stdplug')


class TestPhotonBL(unittest.TestCase):
    def test_create_obj_PhotonBL(self):
        aa = 0.
        gg = gyoto.std.KerrBL()
        gg.spin(aa)
        gg2 = gg.clone()
        ph = gyoto.core.Photon()

        # Photon mass should be zero
        self.assertEqual(ph.getMass(), 0.0)

        # Attach metric
        ph.Metric = gg

        # initial conditions :
        ri = 6
        thetai = np.pi / 2.
        phii = 0.

        ti = 0.
        pri = 0.  # canonical momentum
        pthetai = 0.

        # yinit[5] and [7] are ignored by MakeCoord
        yinit = [ti, ri, thetai, phii, 0., pri, pthetai, 0.]
        cst = [1, 0.921103, 2., 0., 0]  # 4 Kerr cst of motion a, E, L, Q and 1/Q or 0.
        coord = np.zeros(8)
        gg.MakeCoord(yinit, cst, coord)

        self.assertLess(np.max(abs((coord
                                    - [0, 6, 1.5708, 0, 1.38165, 0, 0, 0.0555556]))),
                        1e-5,
                        "computed coordinate is not as expected")

        gg.spin(0.95)
        gg.MakeCoord(yinit, cst, coord)

        pos = coord[:4]
        v = coord[5:] / coord[4]

        st = gyoto.std.Star()
        st.metric(gg)
        st.radius(1.)
        st.setInitCoord(pos, v)
        st.xFill(770.)

        # Ray tracing

        print("Checking gyoto_Metric_setObserverPos: ")
        screen = gyoto.core.Screen()
        screen.metric(gg)
        screen.setObserverPos([1000., 100., 0.05, 0.])
        gg.spin(0.)
        gg.deltaMin(40)

        orbit = gyoto.std.Star()
        orbit.metric(gg)
        orbit.radius(2)
        orbit.setInitCoord((600., 6., 1.57, 0.), (0., 0., 0.068041381745))

        N = 21
        delta = np.pi / (10. * N)
        screen.FieldOfView = np.pi / 10.
        screen.Resolution = 21
        gg.GenericIntegrator = True
        i = 15
        j = 9
        xscr = delta * (i - (N + 1) / 2.)
        yscr = delta * (j - (N + 1) / 2.)

        # Checking gyoto_Photon_setInitialCondition
        ph2 = gyoto.core.Photon()
        coordout = np.empty(8)
        screen.getRayTriad(i, j, coordout)
        ph2.setInitialCondition(gg, orbit, coordout)

        ph1 = gyoto.core.Photon()
        ph1.setInitialCondition(gg, orbit, screen, -xscr, yscr)

        self.assertEqual(abs(np.asarray(ph1.InitCoord) - ph2.InitCoord).max(), 0,
                         "the two versions of setInitialCondition don't agree")

        hit1 = ph1.hit()
        hit2 = ph2.hit()
        self.assertEqual(hit1, hit2)
        self.assertEqual(hit1, True)

        n1 = ph1.get_nelements()
        n2 = ph2.get_nelements()
        self.assertEqual(n1, n2)

        hitmap1 = np.zeros((N, N))
        hitmap2 = np.zeros((N, N))

        for i in range(1, N):
            xscr = delta * (i - (N + 1) / 2.)
            for j in range(1, N):
                yscr = delta * (j - (N + 1) / 2.)
                screen.getRayTriad(i, j, coordout)
                ph2.InitCoord = coordout
                ph1.setInitialCondition(ph.metric(), ph.astrobj(),
                                        screen, -xscr, yscr)

                ph2.delta(1.)
                ph1.delta(1.)
                hitmap1[i, j] = ph1.hit()
                hitmap2[i, j] = ph2.hit()

        self.assertTrue( (hitmap1 == hitmap2).all() )

if __name__ == '__main__':
    unittest.main()
