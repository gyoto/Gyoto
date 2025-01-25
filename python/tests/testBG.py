#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest
import numpy as np
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
        pos = gyoto.core.array_double(4)
        pos[0] = 0.
        pos[1] = 50.
        pos[2] = np.pi/2
        pos[3] = 0.
        
        g = gyoto.core.array_double(16)  # 4x4 matrix flattened
        self.metric.gmunu(g, pos)
        
        # Test all non-zero components with 6 decimal places
        self.assertAlmostEqual(g[0], -1.0, places=6)          # g[0,0]
        self.assertAlmostEqual(g[3], 0.027285, places=5)      # g[0,3]
        self.assertAlmostEqual(g[5], 1.0, places=6)           # g[1,1]
        self.assertAlmostEqual(g[10], 2500.0, places=4)       # g[2,2]
        self.assertAlmostEqual(g[12], 0.027285, places=5)     # g[3,0]
        self.assertAlmostEqual(g[15], 2500.0, places=4)       # g[3,3]

if __name__ == '__main__':
    unittest.main()
