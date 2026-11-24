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
import sys
import unittest
import gyoto.core
import numpy as np
import gyoto.metric
import gyoto.astrobj
import gyoto.spectrum
import gyoto.spectrometer
# import inspect
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

gyoto.core.requirePlugin('stdplug')

class TestPatternDisk(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts_dir = "tests/artifacts/"
        os.makedirs(cls.artifacts_dir, exist_ok=True)

        dirname, filename = os.path.split(os.path.abspath(__file__))

        # Create PatternDisk
        # Create opacity and intensity grids as numpy arrays.
        # Get pointers in a format that Gyoto understands.
        # Warning: here we assume that size_t is the same as uint64.
        gridshape = np.asarray((1, 3, 11), np.uint64)
        pgridshape = gyoto.core.array_size_t.fromnumpy1(gridshape)
        opacity = np.zeros(gridshape)
        popacity = gyoto.core.array_double.fromnumpy3(opacity)
        opacity[:, 0::2, 0::2] = 100.
        opacity[:, 1::2, 1::2] = 100.
        intensity = opacity * 0. + 1.
        pintensity = gyoto.core.array_double.fromnumpy3(intensity)
        # Create PatternDisk, attach grids, set some parameters
        metric = gyoto.std.KerrBL(Mass = (4e6, "sunmass"))
        cls.pd = gyoto.std.PatternDisk(
            Metric = metric,
            RMax = 50,
            InnerRadius = 3,
            OuterRadius = 28,
        )
        cls.pd.copyIntensity(pintensity, pgridshape)
        cls.pd.copyOpacity(popacity, pgridshape)
        cls.pd.repeatPhi(8)
        cls.sc = gyoto.core.Scenery(
            Metric = metric,
            Astrobj = cls.pd,
            Screen = gyoto.core.Screen(
                Metric = metric,
                Resolution = 64,
                Time = (1000., "geometrical_time"),
                Distance = (100., "geometrical"),
                FieldOfView = 30. / 100.,
                Inclination = (110., "degree"),
                PALN = (180., "degree")
            )
        )
        # needed for tests functions below
        cls.gridshape = gridshape
        cls.pgridshape = pgridshape
        cls.opacity = opacity
        cls.popacity = popacity
        cls.intensity = intensity
        cls.pintensity = pintensity

    def test_PatternDisk_create_scenery(self):
        self.assertIsInstance(self.sc, gyoto.core.Scenery)

    def test_fio_PatternDisk_simple(self):
        self.pd.fitsWrite("!test_patterndisk.fits.gz")
        self.sc.write("test_patterndisk.xml")
        gyoto.core.Scenery.read("test_patterndisk.xml")
        os.unlink("test_patterndisk.fits.gz")
        os.unlink("test_patterndisk.xml")

    def test_fio_PatternDisk_prefix(self):
        self.pd.fitsWrite("!test_patterndisk2.fits.gz",
                          self.artifacts_dir)
        self.sc.write(self.artifacts_dir + "test_patterndisk2.xml")
        sc2 = gyoto.core.Scenery.read(self.artifacts_dir
                                      + "test_patterndisk2.xml")
        self.assertEqual(sc2.screen().dMax(), self.sc.screen().dMax(),
                         "dmax was not conserved when RW XML file")
        self.assertEqual(sc2.tMin(), self.sc.tMin(),
                         "tmin was not conserved when RW XML file")

    def test_fio_PatternDisk3_different_path(self):
        os.makedirs(self.artifacts_dir + "fits/", exist_ok=True)
        os.makedirs(self.artifacts_dir + "xml/", exist_ok=True)
        self.pd.fitsWrite("!../fits/test_patterndisk3.fits.gz",
                          self.artifacts_dir + "xml/")

        TestPatternDisk.sc.write(self.artifacts_dir + "xml/test_patterndisk3.xml")

        gyoto.core.Scenery.read(self.artifacts_dir + "xml/test_patterndisk3.xml")


if __name__ == '__main__':
    unittest.main()
