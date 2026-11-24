#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Copyright 2013-2014 Thibaut Paumard
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
import gyoto.core
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

gyoto.core.requirePlugin('stdplug')

GYOTO_EXAMPLES_DIR = os.environ.get("GYOTO_EXAMPLES_DIR", "../doc/examples/")
GYOTO_MAX_THREADS =  int(os.environ.get("GYOTO_MAX_THREADS", os.cpu_count()))

class TestStarTrace(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts_dir = "tests/artifacts/"
        os.makedirs(cls.artifacts_dir, exist_ok=True)

    def test_StarTrace(self):
        sc = gyoto.utils.readScenery(GYOTO_EXAMPLES_DIR
                                     + "example-moving-star.xml")
        sc.screen().mask(None)  # make sure no mask is set yet
        st = sc.astrobj()

        # Change resolution
        sc.Screen.Resolution = 32
        sc.Astrobj.Radius = 2

        # Instanciating StarTrace from Star
        st = gyoto.std.Star(st)
        stt = gyoto.std.StarTrace(st, 600, 800)

        stt.Adaptive = 0
        stt.Delta = 1.
        stt.OpticallyThin = False

        sc.Astrobj = stt
        sc.NThreads = GYOTO_MAX_THREADS

        # Ray-tracing StarTrace
        sc.Quantities = "Intensity"
        mask = sc[:, :]["Intensity"]
        sc.Astrobj = st

        # Ray-tracing Star without mask
        im1 = sc[:, :]["Intensity"]

        # Ray-tracing Star with mask
        pmask=gyoto.core.array_double.fromnumpy2(mask)
        sc.Screen.mask(pmask)
        im2 = sc[:, :]["Intensity"]

        # compare
        self.assertLess(abs((im1 - im2)/im1.max()).max(), 1e-15)

        # Create PDF file
        self.file_output = PdfPages(self.artifacts_dir + 'test_startrace.pdf')

        fig = plt.figure()
        plt.imshow(mask, origin='lower')
        self.file_output.savefig()

        plt.imshow(im1, origin='lower')
        self.file_output.savefig()

        plt.imshow(im2, origin='lower')
        self.file_output.savefig()

        plt.close()
        self.file_output.close()

if __name__ == '__main__':
    unittest.main()
