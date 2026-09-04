#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Copyright 2014 Frederic Vincent, Thibaut Paumard
#    Copyright 2026 Julien Brulé, Thibaut Paumard
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
import os
import unittest
import gyoto.core
import numpy as np
import gyoto.metric
import gyoto.astrobj
import gyoto.spectrum
import gyoto.spectrometer
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

gyoto.core.requirePlugin('stdplug')

class TestMetric(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts_dir = "tests/artifacts/"
        os.makedirs(cls.artifacts_dir, exist_ok=True)

    def test_DirectionalDisk(self):
        nnu = 2
        ni = 2
        nr = 10
        naxes = (nr, ni, nnu)

        # ! DirectionalDisk convention of axes ordering is reversed
        # ! compared to default Python convention
        gridshape = np.asarray(naxes[::-1], np.uint64)
        pgridshape = gyoto.core.array_size_t.fromnumpy1(gridshape)
        intensity = np.zeros(naxes)
        pintensity = gyoto.core.array_double.fromnumpy3(intensity)

        # Fill intensity array
        intensity[0::2, 0, :] = 10.
        intensity[1::2, 0, :] = 7.
        intensity[2::8, 0, :] = 4.
        intensity[8::10, 0, :] = 1.
        intensity[:, 1, :] = 0.1 * intensity[:, 0, :]

        # Create grid arrays
        freq = np.zeros(nnu)
        pfreq = gyoto.core.array_double.fromnumpy1(freq)
        cosi = np.zeros(ni)
        pcosi = gyoto.core.array_double.fromnumpy1(cosi)
        radius = 5. * (np.arange(nr) + 1)
        pradius = gyoto.core.array_double.fromnumpy1(radius)

        # Specials Values
        freqobs = 1e18
        freq[0] += freqobs * 100.
        freq[1] += freqobs / 100.  # in decreasing order
        cosi += 0.5
        cosi[0] = 0.3
        cosi[1] = 0.6

        # create metric
        metric = gyoto.std.KerrBL()
        metric.mass(4e6 * gyoto.core.GYOTO_SUN_MASS)

        # create object
        directional_disk = gyoto.std.DirectionalDisk()
        directional_disk.copyIntensity(pintensity, pgridshape)
        directional_disk.copyGridFreq(pfreq, nnu)
        directional_disk.copyGridCosi(pcosi, ni)
        directional_disk.copyGridRadius(pradius, nr)
        directional_disk.metric(metric)
        directional_disk.rMax(100)
        directional_disk.innerRadius(6)
        directional_disk.outerRadius(50)

        # create screen
        screen = gyoto.core.Screen()
        screen.resolution(64)
        screen.time(1000. * metric.unitLength() / gyoto.core.GYOTO_C)
        screen.distance(100. * metric.unitLength())
        screen.set("FieldOfView", 1.1)
        screen.inclination(110. / 180. * np.pi)
        screen.set("PALN", np.pi)
        screen.freqObs(freqobs)

        # create scenery
        sc = gyoto.core.Scenery()
        sc.metric(metric)
        sc.screen(screen)
        sc.astrobj(directional_disk)

        # write arrays
        directional_disk.fitsWrite("!test_directionaldisk.fits.gz",
                                   self.artifacts_dir)

        # write XML
        gyoto.core.Factory(sc).write(
            self.artifacts_dir + "test_directionaldisk.xml")

        # Read the scenery we just saved
        sc2 = gyoto.core.Factory(
            self.artifacts_dir + "test_directionaldisk.xml"
        ).scenery()

        # check FieldOfView
        self.assertEqual(sc2.screen().get("FieldOfView"),
                         sc.screen().get("FieldOfView"), "different fov")

        # Create PDF file
        file_output = PdfPages(self.artifacts_dir + 'test_directionaldisk.pdf')

        # Compare DirectionalDisk
        # compare shape
        dd2 = gyoto.std.DirectionalDisk(sc2.astrobj())
        pgridshape2 = gyoto.core.array_size_t(3)
        dd2.getIntensityNaxes(pgridshape2)
        for k in range(3):
            assert pgridshape2[k] == pgridshape[k], "shape of grid changed"
            bufsize = gridshape.prod()

        # compare intensity
        buf = gyoto.core.array_double.frompointer(dd2.getIntensity())
        for k in range(bufsize):
            assert buf[k] == pintensity[k], "Intensity changed"

        # compare freq
        buf = gyoto.core.array_double.frompointer(dd2.getGridFreq())
        for k in range(nnu):
            assert buf[k] == pfreq[k], "gridfreq changed"

        # compare cosi
        buf = gyoto.core.array_double.frompointer(dd2.getGridCosi())
        for k in range(ni):
            assert buf[k] == pcosi[k], "gridcosi changed"

        # compare radius
        buf = gyoto.core.array_double.frompointer(dd2.getGridRadius())
        for k in range(nr):
            assert buf[k] == pradius[k], "gridradius changed"

        # ray-trace
        frame = sc[:,:]["Intensity"]

        # plot and close PDF
        plt.figure()
        plt.imshow(frame, origin='lower')
        file_output.savefig()
        plt.close()
        file_output.close()

if __name__ == '__main__':
    unittest.main()
