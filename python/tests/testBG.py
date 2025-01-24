#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import gyoto.core
import gyoto.std
import unittest

class TestBGMetric(unittest.TestCase):
    def setUp(self):
        """Set up a BG metric instance before each test"""
        self.metric = gyoto.core.Metric("BG")
        # Set standard parameter values
        self.metric.set("V0", 0.000733333)
        self.metric.set("R", 100.0)
        self.metric.set("r0", 1.0)
        
    def test_parameters(self):
        """Test parameter getting/setting"""
        self.assertAlmostEqual(self.metric.get("V0"), 0.000733333)
        self.assertAlmostEqual(self.metric.get("R"), 100.0)
        self.assertAlmostEqual(self.metric.get("r0"), 1.0)
        
    def test_metric_components(self):
        """Test metric tensor components at a specific point"""
        # Test at r=50kpc, theta=pi/2 (equatorial plane)
        pos = np.array([0., 50., np.pi/2, 0.])
        g = np.zeros((4,4))
        self.metric.gmunu(g, pos)
        
        # Test known values for g_tt, g_tphi, g_rr, g_theta_theta, g_phi_phi
        self.assertAlmostEqual(g[0][0], -1.0)  # g_tt should be -1
        self.assertAlmostEqual(g[1][1], 1.0)   # g_rr should be 1
        
    def test_christoffel(self):
        """Test Christoffel symbols at a specific point"""
        pos = np.array([0., 50., np.pi/2, 0.])
        christ = np.zeros((4,4,4))
        self.metric.christoffel(christ, pos)
        
        # Test some known Christoffel symbols
        # Add specific tests based on known values
        
    def test_coordinate_singularities(self):
        """Test behavior near coordinate singularities"""
        # Test metric components as r ? 0
        pos_near_center = np.array([0., 1e-6, np.pi/2, 0.])
        g = np.zeros((4,4))
        self.metric.gmunu(g, pos_near_center)
        # Add assertions for expected behavior

if __name__ == '__main__':
    unittest.main()