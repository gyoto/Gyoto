#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright 2014 Frederic Vincent, Thibaut Paumard
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
import ctypes
import gyoto.metric
import gyoto.astrobj
import gyoto.spectrum
import gyoto.spectrometer
# import inspect
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

gyoto.core.requirePlugin('stdplug')

# from yorick
# include "check-helpers.i"


# create object
# create metric
# create parameters
# associate metric /param to the object
# create scenery

# change pdf output filename to a standard os.path python
# add a delete rule for the Makefile


class TestMetric(unittest.TestCase):
    def test_DirectionalDisk(self):
        GYOTO_PYTHON_DIR = "./"
        # GYOTO_EXAMPLES_DIR = GYOTO_PYTHON_DIR + "../doc/examples/"
        GYOTO_ARTIFACTS_DIR = GYOTO_PYTHON_DIR + "tests/artifacts/"
        GYOTO_YORICK_DIR = GYOTO_PYTHON_DIR + "../yorick/"
        # gyoto.core.verbose(10)
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

        # create object
        directional_disk = gyoto.std.DirectionalDisk()
        # create metric
        metric = gyoto.std.KerrBL()
        metric.mass(4e6 * gyoto.core.GYOTO_SUN_MASS)
        print("Creating DirectionalDisk...")

        directional_disk.copyIntensity(pintensity, pgridshape)
        directional_disk.copyGridFreq(pfreq, nnu)
        directional_disk.copyGridCosi(pcosi, ni)
        directional_disk.copyGridRadius(pradius, nr)
        directional_disk.metric(metric)
        directional_disk.rMax(100)
        directional_disk.innerRadius(6)
        directional_disk.outerRadius(50)

        print("Printing DirectionalDisk:")
        screen = gyoto.core.Screen()
        screen.resolution(64)
        screen.time(1000. * metric.unitLength() / gyoto.core.GYOTO_C)
        screen.distance(100. * metric.unitLength())
        screen.set("FieldOfView", 1.1)
        screen.inclination(110. / 180. * np.pi)
        screen.set("PALN", np.pi)
        screen.freqObs(freqobs)
        print("Attaching DirectionalDisk to scenery...")
        sc = gyoto.core.Scenery()
        sc.metric(metric)
        sc.screen(screen)
        sc.astrobj(directional_disk)
        print("Saving data to fits file...")
        # directional_disk.fitsWrite("check-directionaldisk.fits.gz")
        print(f"{GYOTO_ARTIFACTS_DIR + "check-directionaldisk.fits.gz"}")
        directional_disk.fitsWrite("!check-directionaldisk.fits.gz",
                                   GYOTO_ARTIFACTS_DIR)
        # directional_disk.fitsWrite("!check-directionaldisk.fits.gz")

        print("Saving scenery to XML file...")
        gyoto.core.Factory(sc).write(
            GYOTO_ARTIFACTS_DIR + "check-directionaldisk.xml")

        # Reading the scenery we just saved
        print("Reading back scenery...")
        sc2 = gyoto.core.Factory(
            GYOTO_ARTIFACTS_DIR + "check-directionaldisk.xml"
        ).scenery()
        # sc2 = gyoto.core.Factory("check-directionaldisk.xml").scenery()

        # GYOTO_ARTIFACTS_DIR + "check-directionaldisk.xml").scenery()
        # Check
        # Compare Sceneries
        # assert sc2.screen().dMax() == sc.screen().dMax(),
        # "dmax was not conserved when writing and reading XML"
        # assert sc2.tMin() == sc.tMin(),
        # "tmin was not conserved when writing and reading XML"
        self.assertEqual(sc2.screen().get("FieldOfView"),
                         sc.screen().get("FieldOfView"), "different fov")
        print("Removing temporary files...")
        # os.remove(GYOTO_ARTIFACTS_DIR + "check-directionaldisk.xml")
        # os.remove(GYOTO_ARTIFACTS_DIR + "check-directionaldisk.fits.gz")

        file_output = PdfPages(GYOTO_ARTIFACTS_DIR + 'directionalDisk.pdf')

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

        # COMPARING WITH YORICK RESULTS
        # TO BE REMOVED ONCE DEBUGGED
        print("Reading back YORICK scenery...")
        sc2 = gyoto.core.Factory(
            GYOTO_YORICK_DIR + "check-directionaldisk.xml"
        ).scenery()
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

        # END OF TO BE REMOVED ONCE DEBUGGED

        # Ray-trace
        frame = sc[:,:]["Intensity"]
        plt.figure()
        plt.imshow(frame, origin='lower')
        # plt.savefig(frame)
        file_output.savefig()
        plt.close()
        file_output.close()


if __name__ == '__main__':
    unittest.main()
