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
import copy
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
        pd = gyoto.std.PatternDisk()
        pd.copyIntensity(pintensity, pgridshape)
        pd.copyOpacity(popacity, pgridshape)
        pd.innerRadius(3)
        pd.outerRadius(28)
        pd.repeatPhi(8)
        metric = cls.create_metric()
        # see if we have the good one
        pd.metric(metric)
        pd.rMax(50)
        # needed for tests functions below
        cls.gridshape = gridshape
        cls.pgridshape = pgridshape
        cls.opacity = opacity
        cls.popacity = popacity
        cls.intensity = intensity
        cls.pintensity = pintensity
        cls.pd = pd

    @staticmethod
    def create_metric():
        # Create a metric
        metric = gyoto.std.KerrBL()
        metric.mass(4e6, "sunmass")
        return metric

    def create_screen(self):
        # Create screen
        screen = gyoto.core.Screen()
        metric = self.create_metric()
        screen.metric(metric)
        screen.resolution(64)
        screen.time(1000., "geometrical_time")
        screen.distance(100., "geometrical")
        screen.fieldOfView(30. / 100.)
        screen.inclination(110., "degree")
        screen.PALN(180., "degree")
        return screen

    def create_scenery(self):
        # Create Scenery
        sc = gyoto.core.Scenery()
        metric = TestPatternDisk.create_metric(self)
        sc.metric(metric)
        sc.screen(TestPatternDisk.create_screen(self))
        sc.astrobj(TestPatternDisk.pd)

    def test_PatternDisk_create_scenery(self):
        # Create Scenery
        sc = gyoto.core.Scenery()
        metric = self.create_metric()
        sc.metric(metric)
        sc.screen(self.create_screen())
        sc.astrobj(self.pd)
        self.assertIsInstance(sc, gyoto.core.Scenery)
        TestPatternDisk.sc = copy.copy(sc)

    def test_fio_PatternDisk_simple(self):
        TestPatternDisk.pd.fitsWrite("!test_patterndisk.fits.gz")
        gyoto.core.Factory(TestPatternDisk.sc).write("test_patterndisk.xml")
        gyoto.core.Factory("test_patterndisk.xml").scenery()
        os.unlink("test_patterndisk.fits.gz")
        os.unlink("test_patterndisk.xml")

    def test_fio_PatternDisk_prefix(self):
        TestPatternDisk.pd.fitsWrite("!test_patterndisk2.fits.gz",
                                     self.artifacts_dir)
        gyoto.core.Factory(TestPatternDisk.sc).write(
            self.artifacts_dir + "test_patterndisk2.xml")
        sc2 = gyoto.core.Factory(self.artifacts_dir
                                 + "test_patterndisk2.xml").scenery()
        self.assertEqual(sc2.screen().dMax(), self.sc.screen().dMax(),
                         "dmax was not conserved when RW XML file")
        self.assertEqual(sc2.tMin(), self.sc.tMin(),
                         "tmin was not conserved when RW XML file")

    def test_fio_PatternDisk3_different_path(self):
        os.makedirs(self.artifacts_dir + "fits/", exist_ok=True)
        os.makedirs(self.artifacts_dir + "xml/", exist_ok=True)
        TestPatternDisk.pd.fitsWrite("!../fits/test_patterndisk3.fits.gz",
                                     self.artifacts_dir + "xml/")

        gyoto.core.Factory(TestPatternDisk.sc).write(
            self.artifacts_dir + "xml/test_patterndisk3.xml")

        gyoto.core.Factory(self.artifacts_dir
                           + "xml/test_patterndisk3.xml").scenery()


if __name__ == '__main__':
    unittest.main()
