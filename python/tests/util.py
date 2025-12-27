#!/usr/bin/env python
# -*- coding: utf-8 -*-
import unittest
import gyoto
from gyoto.util import convert
from gyoto.core import GYOTO_C, GYOTO_G_OVER_C_SQUARE, GYOTO_SUN_MASS, GYOTO_SUN_RADIUS, GYOTO_KPC

gyoto.core.requirePlugin('stdplug')

class TestConvert(unittest.TestCase):
    def test_convert(self):
        metric=gyoto.std.KerrBL()
        self.assertEqual(convert(1., "geometrical", "geometrical_time"), 1.)
        self.assertEqual(convert(1., "geometrical_time", "geometrical"), 1.)
        self.assertEqual(convert(1., "Hz", "Hz"), 1.)
        self.assertEqual(convert(1., "s", "m"), GYOTO_C)
        self.assertEqual(convert(1., "m", "s"), 1./GYOTO_C)
        self.assertEqual(convert(1., "m", "geometrical", metric), 1./GYOTO_G_OVER_C_SQUARE)
        self.assertEqual(convert(1., "geometrical", "m", metric), GYOTO_G_OVER_C_SQUARE)
        self.assertEqual(convert(1., "s", "geometrical", metric), GYOTO_C/GYOTO_G_OVER_C_SQUARE)
        self.assertEqual(convert(1., "geometrical", "s", metric), GYOTO_G_OVER_C_SQUARE/GYOTO_C)
        self.assertEqual(convert(1., "m", "geometrical_time", metric), 1./GYOTO_G_OVER_C_SQUARE)
        self.assertEqual(convert(1., "geometrical_time", "m", metric), GYOTO_G_OVER_C_SQUARE)
        self.assertEqual(convert(1., "s", "geometrical_time", metric), GYOTO_C/GYOTO_G_OVER_C_SQUARE)
        self.assertEqual(convert(1., "geometrical_time", "s", metric), GYOTO_G_OVER_C_SQUARE/GYOTO_C)
        self.assertEqual(convert(1., "sunmass", "kg", metric), GYOTO_SUN_MASS)
        self.assertEqual(convert(1., "sunradius", "m", metric), GYOTO_SUN_RADIUS)
        self.assertAlmostEqual(convert(1./GYOTO_KPC, "kpc", "m", metric), 1., 5)
