#!/usr/bin/env python
# -*- coding: utf-8 -*-

#
#    Copyright 2012-2014, 2019 Thibaut Paumard
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
#    along with Gyoto.  If not, see <http:# www.gnu.org/licenses/>.
#
import os
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

GYOTO_EXAMPLES_DIR = os.environ.get("GYOTO_EXAMPLES_DIR", "../doc/examples/")
PD_XML = GYOTO_EXAMPLES_DIR + "example-polish-doughnut.xml"

class TestPolishDoughnut(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts_dir = "tests/artifacts/"
        os.makedirs(cls.artifacts_dir, exist_ok=True)

    def test_create_obj_PolishDoughnut(self):
        pd = gyoto.std.PolishDoughnut()
        assert pd

        gg = gyoto.std.KerrBL()
        pd.metric(gg)

        self.assertIsNotNone(pd.metric())
        self.assertEqual(pd.metric().getPointer(), gg.getPointer())

    def test_create_PolishDoughnut_from_file(self):
        sc = gyoto.core.Factory(PD_XML).scenery()
        pd = gyoto.core.Factory(PD_XML).astrobj()
        assert pd

        spect = gyoto.spectrometer.Uniform()
        spect.kind("freqlog")
        spect.nSamples(10)
        spect.band([10, 14], "Hz")
        sc.screen().spectrometer(spect)
        sc.astrobj().opticallyThin(True)

        sc.requestedQuantitiesString("Spectrum[J.m-2.s-1.sr-1.Hz-1]")
        s1 = sc[14, 9]["Spectrum"]
        midpoints = sc.screen().spectrometer().midpoints("Hz")
        widths = sc.screen().spectrometer().widths("Hz")

        # Integrating one bin spectrum with radiative transfer...
        sc.requestedQuantitiesString("BinSpectrum[J.m-2.s-1.sr-1]")
        s2 = sc[14, 9]["BinSpectrum"]
        s22 = np.transpose((s2 ,s2)).flatten()
        channels = sc.Screen.Spectrometer.channelBoundaries()
        chanind =  sc.Screen.Spectrometer.channelIndices()
        chan2 = np.asarray(channels)[np.asarray(chanind)]

        # Create PDF file
        self.file_output = PdfPages(self.artifacts_dir + 'test_polish_doughnut.pdf')
        plt.figure()
        plt.loglog(midpoints, s1 * widths)
        plt.loglog(chan2, s22)
        plt.xlabel("Frequency [Hz]")
        self.file_output.savefig()
        plt.close()
        self.file_output.close()

if __name__ == '__main__':
    unittest.main()
