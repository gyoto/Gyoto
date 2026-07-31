#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#    Copyright 2011, 2014 Thibaut Paumard
#    Copyright 2026 Julien Brulé, Thibaut Paumard
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
import gyoto
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

gyoto.core.requirePlugin('stdplug')

class TestStar(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.artifacts_dir = "tests/artifacts/"
        os.makedirs(cls.artifacts_dir, exist_ok=True)

    def test_create_obj_Star(self):
        st = gyoto.std.Star()
        st2 = st.clone()
        self.assertTrue(st.getPointer())
        self.assertTrue(st2.getPointer())
        self.assertTrue(st.getPointer() != st2.getPointer())

        st.radius(0.5)
        st2.radius(2.)
        self.assertEqual(st.Radius, 0.5)
        self.assertEqual(st2.Radius, 2.)

        metric = gyoto.std.KerrBL()
        st.metric(metric)

        metric.spin(0.)
        self.assertEqual(metric.Spin, 0.)

        metric.spin(0.7)
        self.assertEqual(st.Metric.Spin, 0.7)

        metric.spin(0.995)
        self.assertEqual(metric.Spin, 0.995)

        metric.Mass = 4e6, "sunmass"

        st.setInitCoord((0., 10.791, 1.5708, 0.), (0., 0., 0.0166637))
        self.assertEqual(abs(np.asarray(st.InitCoord)
                             - (0.0, 10.791, 1.5708, 0.0, 1.1264104273542368,
                                0.0, 0.0, 0.018770165438302795)).max(),
                         0.,
                         "3-vel initcoord not correctly transformed to 4-vel"
                         )

        # integrating
        st.integrator("runge_kutta_fehlberg78")
        st.deltaMaxOverR(0.1)
        st.xFill(800.)

        # retrieving projected coordinates
        n = st.get_nelements()
        screen = gyoto.core.Screen()
        screen.metric(st.metric())
        x = np.ndarray(n)
        y = np.ndarray(n)
        z = np.ndarray(n)
        st.getSkyPos(screen, x, y, z)

        # Create PDF file
        self.file_output = PdfPages(self.artifacts_dir + 'test_star.pdf')
        fig = plt.figure()
        plt.plot(x, y)
        ax=fig.axes[0]
        ax.set_xlabel("rigth ascension offset [m]")
        ax.set_ylabel("declination offset [m]")
        ax.set_aspect('equal')
        ax.invert_xaxis()
        self.file_output.savefig()
        plt.close()
        self.file_output.close()

if __name__ == '__main__':
    unittest.main()
