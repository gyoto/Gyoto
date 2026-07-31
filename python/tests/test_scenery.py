#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Copyright 2011, 2014 Thibaut Paumard
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
import os
import unittest
import ctypes
import gyoto.core
import numpy as np
import gyoto.metric
import gyoto.astrobj
import gyoto.spectrum
import gyoto.spectrometer
# import inspect
import matplotlib.pyplot as plt
# from matplotlib.backends.backend_pdf import PdfPages

gyoto.core.requirePlugin('stdplug')

GYOTO_EXAMPLES_DIR = os.environ.get("GYOTO_EXAMPLES_DIR", "../doc/examples/")
GYOTO_MAX_THREADS =  int(os.environ.get("GYOTO_MAX_THREADS", os.cpu_count()))

class TestScenery(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts_dir = "tests/artifacts/"
        os.makedirs(cls.artifacts_dir, exist_ok=True)

    def test_create_obj_Scenery(self):
        sc = gyoto.core.Scenery()
        self.assertTrue(sc.getPointer(), 'failed to instanciate scenery')

        # test various ways of getting maxiter
        max1 = sc.maxiter()
        max2 = sc.MaxIter
        max3 = sc.get("MaxIter")
        self.assertEqual(max1, max2,
                         "sc.maxiter() and sc.MaxIter differ")
        self.assertEqual(max1, max3,
                         "sc.maxiter() and sc.get(\"MaxIter\") differ")

        # test various ways of setting it
        sc.maxiter(25)
        self.assertEqual(sc.MaxIter, 25, "sc.maxiter(...) failed")

        sc.MaxIter = 40
        self.assertEqual(sc.MaxIter, 40, "sc.Maxiter = ... failed")

        sc.set("MaxIter", 60)
        self.assertEqual(sc.MaxIter, 60, "sc.set('Maxiter',...) failed")

        # Set back the original value.
        sc.maxiter(max1)

        gg = gyoto.core.Metric("KerrBL")
        sc.metric(gg)

        ao = gyoto.std.Star()
        ao.metric(gg)
        ao.radius(0.5)
        ao.setInitCoord((0, 6, np.pi / 2., 0), (0, 1e-3, 0))
        sc.astrobj(ao)

        ao2 = sc.astrobj()

        self.assertEqual(ao.getPointer(), ao2.getPointer(),
                         "retrieved astrobj pointer is not same as set")

        sc.screen().time(10000)
        time2 = sc.screen().time()
        self.assertEqual(time2, 10000,
                         "retrieved time is not same as set")

        sc.screen().fieldOfView(np.pi / 4.)
        self.assertEqual(sc.Screen.FieldOfView, np.pi/4.)

        sc.screen().resolution(16)
        self.assertEqual(sc.screen().resolution(), 16)

        sc.screen().inclination(np.pi / 3.)
        self.assertEqual(sc.screen().inclination(), np.pi/3.)

        # just try writing a scenery
        gyoto.core.Factory(sc).write(self.artifacts_dir + "test_scenery.xml")

        # just try reading a scenery
        sc3 = gyoto.core.Factory(GYOTO_EXAMPLES_DIR +
                                 "example-moving-star.xml").scenery()

        sc3.Screen.Resolution = 32
        sc3.Astrobj.Radius = 2

        # Ray-tracing on 1 thread (sc())
        sc3.nThreads(1)
        im1 = sc3.rayTrace(quantities="Intensity")["Intensity"]

        # Ray-tracing on all threads (sc())
        sc3.nThreads(GYOTO_MAX_THREADS)
        im2 = sc3.rayTrace(quantities="Intensity")["Intensity"]

        # The two images are not strictly equal because Delta is
        # always reinitialised in the multithreading scenario
        self.assertLess(abs(im1-im2).max()/im1.max(), 1e-15,
                         "multithreading doesn't yield the same image")

        # Setting a mask
        pim2 = gyoto.core.array_double.fromnumpy2(im2)
        sc3.Screen.mask(pim2)

        # Checking the mask
        pmask = sc3.Screen.mask()
        c_ptr = ctypes.cast(int(pmask), ctypes.POINTER(ctypes.c_double))
        res = sc.Screen.Resolution
        mask = np.ctypeslib.as_array(c_ptr, (32, 32))

        self.assertEqual(abs(mask-im2).max()/im1.max(), 0,
                         "mask is not preserved")

        # Ray-tracing again with a mask
        im3 = sc3.rayTrace(quantities="Intensity")["Intensity"]

        self.assertLess(abs(im3-im2).max()/im1.max(), 1e-15,
                         "masking doesn't yield the same image")

        # initializing a photon
        coord = np.empty(8)
        sc3.Screen.getRayTriad(6, 19, coord)
        ph = gyoto.core.Photon(sc3.Metric, sc3.Astrobj,
                               gyoto.array_double.fromnumpy1(coord))
        self.assertTrue(ph.hit(), "this photon should hit")

        sc = gyoto.utils.readScenery(GYOTO_EXAMPLES_DIR +
                                     "example-complex-astrobj.xml")

        sc.NThreads = GYOTO_MAX_THREADS
        sc.NProcesses = 0

        data = sc[:, :]["Spectrum"]

        r1 = slice(7, 25, 4)
        r2 = slice(1, -2, 3)
        data2 = sc[r2, r1]["Spectrum"]
        self.assertEqual(abs(data[:, r2, r1] - data2).max(), 0,
                         "mask is not preserved")

if __name__ == '__main__':
    unittest.main()
