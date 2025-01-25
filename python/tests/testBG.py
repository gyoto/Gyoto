#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import gyoto.core
import gyoto.std

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
        self.assertAlmostEqual(self.metric.get("V0"), 0.000733333, places=6)
        self.assertAlmostEqual(self.metric.get("R"), 100.0, places=6)
        self.assertAlmostEqual(self.metric.get("r0"), 1.0, places=6)
        
    def test_metric_components(self):
        """Test metric tensor components at a specific point"""
        pos = [0., 50., np.pi/2, 0.]
        g = [[0. for i in range(4)] for j in range(4)]
        self.metric.gmunu(g, pos)
        
        # Test all non-zero components with 6 decimal places
        self.assertAlmostEqual(g[0][0], -1.0, places=6)
        self.assertAlmostEqual(g[0][3], 0.027285, places=5)
        self.assertAlmostEqual(g[1][1], 1.0, places=6)
        self.assertAlmostEqual(g[2][2], 2500.0, places=4)
        self.assertAlmostEqual(g[3][0], 0.027285, places=5)
        self.assertAlmostEqual(g[3][3], 2500.0, places=4)
        
    def test_christoffel(self):
        """Test Christoffel symbols at a specific point"""
        pos = [0., 50., np.pi/2, 0.]
        christ = [[[0. for k in range(4)] for j in range(4)] for i in range(4)]
        self.metric.christoffel(christ, pos)
        
        # Test non-zero components with appropriate precision
        self.assertAlmostEqual(christ[0][0][1], 2.21e-9, places=11)
        self.assertAlmostEqual(christ[0][1][3], 0.000343, places=6)
        self.assertAlmostEqual(christ[1][0][3], -0.000203, places=6)
        self.assertAlmostEqual(christ[1][2][2], -50.0, places=4)
        self.assertAlmostEqual(christ[2][1][2], 0.02, places=3)  # 1/50.0
        self.assertAlmostEqual(christ[3][1][3], 0.02, places=3)
        # ... other Christoffel symbols with appropriate precision

if __name__ == '__main__':
    unittest.main()
