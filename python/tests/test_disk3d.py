#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Copyright 2012, 2014-2015 Frederic Vincent, Thibaut Paumard
#    Copyright 2026 Thibaut Paumard, Julien Brulé
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
import os
import sys
import unittest
import gyoto.core
import numpy as np
import gyoto.metric
import gyoto.astrobj
import gyoto.spectrum
import gyoto.spectrometer

class TestDisk3D(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts_dir = "tests/artifacts/"
        os.makedirs(cls.artifacts_dir, exist_ok=True)

    def test_Disk3D(self):

        # TO BE REMOVED:
        # load Yorick-generated objects
        ysc = gyoto.core.Factory("../yorick/check-disk3d.xml").scenery()
        yd3d = gyoto.std.Disk3D(ysc.Astrobj)
        ypemissquant_shape = gyoto.core.array_size_t(4)
        yd3d.getEmissquantNaxes(ypemissquant_shape)
        print([ypemissquant_shape[i] for i in range(4)])
        ypemissquant = gyoto.core.array_double.frompointer(yd3d.getEmissquant())
        ypvelocity = gyoto.core.array_double.frompointer(yd3d.getVelocity())
        # END TO BE REMOVED

        # emiss
        # ! axis ordering is reversed
        emissquant_shape = np.asarray((1, 2, 10, 10), np.uint64)
        pemissquant_shape = gyoto.core.array_size_t.fromnumpy1(
            emissquant_shape)
        emissquant = np.zeros(emissquant_shape[::-1], np.float64)
        emissquant[0:3, :, 0, :] = 100.
        emissquant[3:10, :, 1, :] = 100.
        pemissquant = gyoto.core.array_double.fromnumpy4(emissquant)
        # TO BE REMOVED:
        for k in range(200):
            assert ypemissquant[k] == pemissquant[k]
        # END TO BE REMOVED

        # velocity
        # ! axis ordering is reversed and leading 3 is implicit
        velocity_shape = np.asarray((3, 2, 10, 10), np.uint64)
        pvelocity_shape = gyoto.core.array_size_t.fromnumpy1(velocity_shape[1:])
        velocity = np.ones(velocity_shape[::-1])
        pvelocity = gyoto.core.array_double.fromnumpy4(velocity)
        # TO BE REMOVED:
        for k in range(600):
            assert ypvelocity[k] == pvelocity[k]
        # END TO BE REMOVED

        # metric
        metric = gyoto.std.KerrBL()
        metric.set("Mass", 4e6 * gyoto.core.GYOTO_SUN_MASS)

        # astrobj
        d3d = gyoto.std.Disk3D()
        d3d.copyEmissquant(pemissquant, pemissquant_shape)
        d3d.copyVelocity(pvelocity, pvelocity_shape)
        d3d.rin(3)
        d3d.rout(28)
        d3d.zmin(1.)
        d3d.zmax(10.)
        d3d.phimin(0.)
        d3d.phimax(2. * np.pi)
        d3d.repeatPhi(8)
        d3d.metric(metric)

        # screen
        screen = gyoto.core.Screen()
        screen.metric(metric)
        screen.resolution(64)
        screen.time(1000. * metric.unitLength() / gyoto.core.GYOTO_C)
        screen.distance(100. * metric.unitLength())
        screen.set("FieldOfView", 30. / 100.)
        screen.inclination(110. / 180. * np.pi)
        screen.set("PALN", np.pi)

        # scenery
        sc = gyoto.core.Scenery()
        sc.metric(metric)
        sc.screen(screen)
        sc.astrobj(d3d)

        # Save Scenery
        d3d.fitsWrite("!test_disk3d.fits.gz", self.artifacts_dir)
        self.assertEqual(d3d.File, "test_disk3d.fits.gz")
        gyoto.core.Factory(sc).write(self.artifacts_dir
                                     + "test_disk3d.xml")

        # Read scenery
        sc2 = gyoto.core.Factory(self.artifacts_dir
                                 + "test_disk3d.xml").scenery()

        # Clone
        sc3 = sc.clone()

        # Note: we check both write+read and clone by comparing them
        # to one another

        # compare the three cube emissquant values
        d3d2 = gyoto.std.Disk3D(sc2.astrobj())
        d3d3 = gyoto.std.Disk3D(sc3.astrobj())

        # EmissquantNaxes
        naxes2 = np.zeros(4, np.uint64)
        pnaxes2 = gyoto.core.array_size_t.fromnumpy1(naxes2)
        d3d2.getEmissquantNaxes(pnaxes2)
        naxes3 = np.zeros(4, np.uint64)
        pnaxes3 = gyoto.core.array_size_t.fromnumpy1(naxes3)
        d3d3.getEmissquantNaxes(pnaxes3)
        self.assertTrue((naxes2 == naxes3).all(),
                        "naxes not conserved in i/o or clone")

        buf_emiss_size = naxes2.prod()
        pemissquant2 = gyoto.core.array_double.frompointer(
            d3d2.getEmissquant())
        pemissquant3 = gyoto.core.array_double.frompointer(
            d3d3.getEmissquant())

        for k in range(buf_emiss_size):
            self.assertEqual(pemissquant2[k], pemissquant3[k],
                             'emissquant not conserved')

            # TO BE REMOVED
            self.assertEqual(pemissquant2[k], ypemissquant[k],
                             'emissquant not conserved')
        for k in range(4):
            assert naxes2[k] == ypemissquant_shape[k]
        # END TO BE REMOVED

        # velocity
        buf_velocity_size = velocity_shape.prod()
        pvelocity2 = gyoto.core.array_double.frompointer(
            d3d2.getVelocity())
        pvelocity3 = gyoto.core.array_double.frompointer(
            d3d3.getVelocity())
        for k in range(buf_velocity_size):
            assert pvelocity2[k] == pvelocity3[k], "Velocity changed"
            # TO BE REMOVED
            self.assertEqual(pvelocity2[k], ypvelocity[k],
                             'velocity not conserved')
            # END TO BE REMOVED

if __name__ == '__main__':
    unittest.main()
