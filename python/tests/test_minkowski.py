#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Copyright 2014 Thibaut Paumard
#    Copyright 2026 Julien Brulé & Thibaut Paumard
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
# import os
import unittest
import gyoto.core
import gyoto.std
import numpy as np
import gyoto.metric
import gyoto.astrobj
import gyoto.spectrum
import gyoto.spectrometer

gyoto.core.requirePlugin('stdplug')

class TestMinkowski(unittest.TestCase):
    def test_Minkowski(self):
        gg = gyoto.std.Minkowski()
        self.assertIsInstance(gg, gyoto.std.Minkowski)

    def test_Minkowski_cartesian(self):

        positions = np.array([[0, 10., 12., 5.],
                              [0, 5., 2., 7.],
                              [0, -10., 0., 50.],
                              [0, 0., 0., 10000.]])

        gg = gyoto.std.Minkowski()

        self.assertEqual(gg.Cartesian, True,
                         "Cartesian should be True by default")

        gyoto.metric.check_christoffel(gg, positions)

        st = gyoto.std.Star()
        st.metric(gg)
        pos = (0., 0., 0., 0.)
        st.DeltaMax = 100.
        st.DeltaMaxOverR = np.inf
        vel = (0., 0., 0.)
        st.setInitCoord(pos, vel)
        st.xFill(1000.)

        n = st.get_nelements()
        t = np.ndarray(n)
        x = np.ndarray(n)
        y = np.ndarray(n)
        z = np.ndarray(n)
        tdot = np.ndarray(n)
        xdot = np.ndarray(n)
        ydot = np.ndarray(n)
        zdot = np.ndarray(n)

        st.get_t(t)
        self.assertTrue(max(t) > 1000 and min(t) < 1000)
        st.getCoord(t, x, y, z, tdot, xdot, ydot, zdot)

        self.assertEqual(np.max(abs(x)), 0.)
        self.assertEqual(np.max(abs(y)), 0.)
        self.assertEqual(np.max(abs(z)), 0.)
        self.assertEqual(np.max(abs(tdot-1.)), 0.)
        self.assertEqual(np.max(abs(xdot)), 0.)
        self.assertEqual(np.max(abs(ydot)), 0.)
        self.assertEqual(np.max(abs(zdot)), 0.)

    def test_Minkowski_cartesian_other_star(self):
        gg = gyoto.std.Minkowski()

        pos = (0., 4., 2., 6.)
        vel = (0.5, 0.2, 0.4)
        st = gyoto.std.Star()
        st.metric(gg)
        st.setInitCoord(pos, vel)
        st.xFill(1000.)
        n = st.get_nelements()

        t = np.ndarray(n)
        x = np.ndarray(n)
        y = np.ndarray(n)
        z = np.ndarray(n)
        tdot = np.ndarray(n)
        xdot = np.ndarray(n)
        ydot = np.ndarray(n)
        zdot = np.ndarray(n)

        st.get_t(t)
        st.getCoord(t, x, y, z, tdot, xdot, ydot, zdot)

        data = np.column_stack((t, x, y, z))
        maxerr = (np.max(
            np.abs(minmax(
                data[:, 1:4]
                - np.expand_dims(data[:, 0], 1)
                * vel[:] - pos[1:4]))))
        self.assertTrue(maxerr < 1e-10,
                        "integration produced wrong results")

    def test_Minkowski_cartesian_photon(self):
        gg = gyoto.std.Minkowski()

        pos = (0., 4., 2., 6.)
        vel = (0.5, 0.2, 0.4)
        photon = gyoto.core.Photon()
        photon.metric(gg)

        # use nullifyCoord to compute dt/dtau
        coord = np.concatenate((pos, (0,), vel))
        gg.nullifyCoord(coord)
        photon.initCoord(coord)

        photon.xFill(1000.)
        n = photon.get_nelements()

        t = np.ndarray(n)
        x = np.ndarray(n)
        y = np.ndarray(n)
        z = np.ndarray(n)
        tdot = np.ndarray(n)
        xdot = np.ndarray(n)
        ydot = np.ndarray(n)
        zdot = np.ndarray(n)

        photon.get_t(t)
        photon.getCoord(t, x, y, z, tdot, xdot, ydot, zdot)

        data = np.column_stack((t, x, y, z, tdot, xdot, ydot, zdot))
        nrows = data.shape[0]
        norm = np.zeros(nrows)
        for i in range(nrows):
            norm[i] = gg.ScalarProd(data[i, 0:4],
                                    data[i, 4:8], data[i, 4:8])
        self.assertTrue(max(norm) < 1e-10, "norm was not conserved")

        maxerr = (np.max(
            np.abs(
                minmax(data[:, 1:4] -
                       np.expand_dims(data[:, 0], 1)
                       * np.expand_dims(vel, 0)
                       / np.expand_dims(data[:, 4], 1)
                       - np.expand_dims(pos[1:], 0)))))
        self.assertTrue(maxerr < 1e-10,
                        "integration produced wrong results")

    def test_Minkowski_spherical_coordinates(self):
        gg = gyoto.std.Minkowski()
        gg.spherical(True)
        self.assertEqual(gg.Spherical, True,
                         "could not set Spherical to True")
        positions=[[0, 10., np.pi/2., 0.],
                   [100., 10., np.pi/4., 2.],
                   [1000., 100., 1., 1.]]
        gyoto.metric.check_christoffel(gg, positions)


    def test_Minkowski_motionless_star(self):
        gg = gyoto.std.Minkowski()
        gg.spherical(True)

        pos = (0., 4., 2., 6.)
        vel = (0., 0., 0.)
        st = gyoto.std.Star()
        st.metric(gg)
        st.setInitCoord(pos, vel)

        st.xFill(1000.)
        n = st.get_nelements()

        t = np.zeros(n)
        r = np.zeros(n)
        theta = np.zeros(n)
        phi = np.zeros(n)
        tdot = np.zeros(n)
        rdot = np.zeros(n)
        thetadot = np.zeros(n)
        phidot = np.zeros(n)

        st.get_t(t)
        st.getCoord(t, r, theta, phi, tdot, rdot, thetadot, phidot)

        self.assertEqual(np.max(abs(r-4.)), 0.)
        self.assertEqual(np.max(abs(theta-2.)), 0.)
        self.assertEqual(np.max(abs(phi-6.)), 0.)
        self.assertEqual(np.max(abs(tdot-1.)), 0.)
        self.assertEqual(np.max(abs(rdot)), 0.)
        self.assertEqual(np.max(abs(thetadot)), 0.)
        self.assertEqual(np.max(abs(phidot)), 0.)

    def test_Minkowski_moving_star(self):
        gg = gyoto.std.Minkowski()
        gg.spherical(True)

        pos = [0., 10.791, np.pi / 2., 0]
        vel = [0., 0., 0.016664]
        st = gyoto.std.Star()
        st.metric(gg)
        st.setInitCoord(pos, vel)

        st.adaptive(True)
        st.delta(0.1)

        st.xFill(1000.)

        n = st.get_nelements()
        t = np.ndarray(n)
        x = np.ndarray(n)
        y = np.ndarray(n)
        z = np.ndarray(n)

        xp = np.ndarray(n)
        yp = np.ndarray(n)
        zp = np.ndarray(n)

        st.get_t(t)
        st.getCartesian(t, x, y, z, xp, yp, zp)

        self.assertLess(abs(xp - xp[0]).max(), 1e-10,
                        "x velocity varied")
        self.assertLess(abs(yp - yp[0]).max(), 1e-10,
                        "y velocity varied")
        self.assertLess(abs(zp - zp[0]).max(), 1e-10,
                        "z velocity varied")

        self.assertLess(abs((x[1:] - x[:-1]) / (t[1:] - t[:-1])
                            - xp[0]).max(), 1e-10,
                        "x velocity different from x variation")
        self.assertLess(abs((y[1:] - y[:-1]) / (t[1:] - t[:-1])
                            - yp[0]).max(), 1e-10,
                        "y velocity different from y variation")
        self.assertLess(abs((z[1:] - z[:-1]) / (t[1:] - t[:-1])
                            - zp[0]).max(), 1e-10,
                        "z velocity different from z variation")

def minmax(data):
    data = data.flatten()
    max = min = data[0]
    for i in range(1, len(data), 2):
        try:
            a, b = data[i:i + 2]
        except ValueError:
            a = b = data[i]
        if b > a:
            a, b = b, a
        if a > max:
            max = a
        if b < min:
            min = b
    return [min, max]

if __name__ == '__main__':
    unittest.main()
