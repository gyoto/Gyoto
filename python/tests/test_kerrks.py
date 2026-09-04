#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright 2026 Thibaut Paumard, Julien Brulé
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

import unittest
import gyoto.core
import numpy as np
import gyoto.metric
import gyoto.astrobj
import gyoto.spectrum
import gyoto.spectrometer

gyoto.core.requirePlugin('stdplug')

class TestKerrKS(unittest.TestCase):

    def construct_kerrks(self):
        gg = gyoto.std.KerrKS()
        self.assertIsInstance(gg, gyoto.std.KerrKS)

    def test_setter_getter_kerrks(self):
        gg = gyoto.std.KerrKS()
        gg.set("Spin", 0.)
        self.assertEqual(gg.get("Spin"), 0.)

    def test_KerrKS(self):
        gg = gyoto.std.KerrKS()
        gg.spin(0.)

        positions = ((0, 10., 12., 5.),
                     (0, 5., 2., 7.),
                     (0, -10., 0., 50.),
                     (0, 0., 0., 10000.))

        # check_gmunu, gg, positions;
        # check_gmunu_up, gg, positions;
        # already performed in tests/std.py

        pos = [0, 10., 12., 5.]
        dk1 = self.dkfunc(gg, pos)
        dk3 = self.dknum(gg, pos)
        np.testing.assert_array_almost_equal(dk1, dk3)

        pos = np.array([0, 10., 12., 5.])
        df1 = self.dffunc(gg, pos)
        df3 = self.dfnum(gg, pos)
        np.testing.assert_array_almost_equal(df1, df3)

        a = gg.jacobian(pos)
        b = gyoto.metric.jacobian_numerical(gg, pos, epsilon=1e-06)
        np.testing.assert_array_almost_equal(a, b)

        pos = (0, 10., 12., 5.)
        gyoto.metric.christoffel_numerical(gg, pos, epsilon=1e-06)

    def dkfunc(self, gg, pos):
        spin_ = gg.Spin
        a2_ = spin_**2
        x, y, z = pos[1:4]
        x2=x*x; y2=y*y; z2=z*z; a2z2=a2_*z2
        tau=x2+y2+z2-a2_;
        rho2=tau*tau+4.*a2z2; rho=np.sqrt(rho2)
        r2=0.5*(tau+rho)
        r=np.sqrt(r2); r3=r2*r; r4=r2*r2; r2_a2=r2+a2_
        rx_ay=r*x+spin_*y; ry_ax=r*y-spin_*x
        f=2.*r3/(r4+a2_*z2); fr2=f*r2

        frac1=1./(r2_a2*r2_a2*rho)
        frac2=z/(r2_a2*r*rho)
        frac3=-z/(r*rho)

        return np.asarray([
	    # d/dt
	    [0., 0., 0., 0.],
	    # d/dx
	    [
	        0.,
	        (
                    r3*(x2+rho)
                    - rx_ay*x*(x2+y2+z2+rho)+a2_*(rx_ay*x+r*(x2+rho))
                ) * frac1,
	        (
                    x*(r3*y+a2_*(ry_ax+r*y)-ry_ax*(x2+y2+z2))
                    - (spin_*r2_a2+ry_ax*x)*rho
                ) * frac1,
	        x * frac3
	    ],
	    # d/dy
	    [
	        0.,
	        (
                    a2_*(rx_ay+r*x)*y+r2_a2*spin_*rho
                    - y*(-r3*x+rx_ay*(x2+y2+z2+rho))
                ) * frac1,
	        (
                    r3*(y2+rho)
                    - ry_ax*y*(x2+y2+z2+rho)
                    + a2_*(ry_ax*y+r*(y2+rho))
                ) * frac1,
	        y * frac3
	    ],
	    # d/dz
	    [
	        0.,
	        ((a2_-r2)*x-2*spin_*r*y)*frac2,
	        ((a2_-r2)*y+2*spin_*r*x)*frac2,
	        (2.*r2- (z2*(a2_ + x2 + y2 + z2 + rho))/rho)/(2.*r3)
	    ]
        ])

    def dknum(self, gg, pos, eps=1e-10):
        grad=np.zeros((4, 4))

        for i in range(4):
            delta = np.zeros(4)
            delta[i] = eps
            grad[:, i] = (
                self.kfunc(gg, pos + delta)
                - self.kfunc(gg, pos - delta)
            ) / (2. * eps)

        return np.asarray(grad)

    def kfunc(self, gg, pos):
        spin_ = gg.Spin; a2_ = spin_**2
        x, y, z = pos[1:4]
        x2=x*x; y2=y*y; z2=z*z; a2z2=a2_*z2
        tau=x2+y2+z2-a2_
        rho2=tau*tau+4.*a2z2; rho=np.sqrt(rho2)
        r2=0.5*(tau+rho)
        r=np.sqrt(r2); r3=r2*r; r4=r2*r2; r2_a2=r2+a2_
        rx_ay=r*x+spin_*y; ry_ax=r*y-spin_*x
        f=2.*r3/(r4+a2_*z2); fr2=f*r2
        return np.asarray([
            1.,
            rx_ay/r2_a2,
            ry_ax/r2_a2,
            z/r
        ])

    def ffunc(self, gg, pos):
        spin_ = gg.Spin; a2_ = spin_**2
        x, y, z = pos[1:4]
        x2=x*x; y2=y*y; z2=z*z; a2z2=a2_*z2
        tau=x2+y2+z2-a2_
        rho2=tau*tau+4.*a2z2; rho=np.sqrt(rho2)
        r2=0.5*(tau+rho)
        r=np.sqrt(r2); r3=r2*r; r4=r2*r2; r2_a2=r2+a2_
        rx_ay=r*x+spin_*y; ry_ax=r*y-spin_*x
        f=2.*r3/(r4+a2_*z2); fr2=f*r2
        return    f;

    def dffunc(self, gg, pos):
        spin_ = gg.Spin; a2_ = spin_**2
        x, y, z = pos[1:4]
        x2=x*x; y2=y*y; z2=z*z; a2z2=a2_*z2
        tau=x2+y2+z2-a2_
        rho2=tau*tau+4.*a2z2; rho=np.sqrt(rho2)
        r2=0.5*(tau+rho)
        r=np.sqrt(r2); r3=r2*r; r4=r2*r2; r2_a2=r2+a2_
        rx_ay=r*x+spin_*y; ry_ax=r*y-spin_*x
        f=2.*r3/(r4+a2_*z2); fr2=f*r2
        a4=a2_*a2_

        r4_a2z2=r4+a2z2;
        temp=-(2.*r3*(r4-3.*a2z2))/(r4_a2z2*r4_a2z2*rho);
        x2_y2_z2=x2+y2+z2;
        temp2=(a4+2.*r2*x2_y2_z2 - a2_* (x2_y2_z2 - 4.* z2 + rho));
        return np.asarray([
            0.,
	    x*temp,
	    y*temp,
            -((4.*r*z*(2.* a4*a2_ + (a2_ + 2.*r2)*x2_y2_z2*x2_y2_z2
                       + a4*(-3.*x2 - 3.*y2 + z2 - 2.*rho)
                       + a2_*(x2 + y2 - z2)*rho))/(rho*temp2*temp2)
              )
        ])

    def dfnum(self, gg, pos, eps=1e-10):
        grad=np.zeros(4);

        for i in range(4):
            delta = np.zeros(4)
            delta[i] = eps
            grad[i] = (
                self.ffunc(gg, pos + delta)
                - self.ffunc(gg, pos - delta)
            ) / (2. * eps)

        return np.asarray(grad)

if __name__ == '__main__':
    unittest.main()
