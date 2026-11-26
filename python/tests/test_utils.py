#!/usr/bin/env python
# -*- coding: utf-8 -*-
import unittest
import numpy

import gyoto
from gyoto.utils import convert
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

class TestMatrix4Invert(unittest.TestCase):
    """Regression tests for gyoto.core.matrix4Invert (Gyoto::matrix4Invert)."""
    # This test case was prepared by Mistral AI after fixing a major
    # bug in matrix4Invert

    SEED = 20261103  # fixed seed => reproducible random tests

    def setUp(self):
        # Local generator: reproducible, and does not touch numpy's global state
        # (default_rng/PCG64 stream is stable across numpy versions)
        self.rng = numpy.random.default_rng(self.SEED)

    def checkInverse(self, A, rtol=1e-12, atol=1e-12):
        # SWIG IN_ARRAY2 typemap expects a C-contiguous 4x4 float64 array
        A = numpy.ascontiguousarray(A, dtype=numpy.float64)
        A_copy = A.copy()
        Am1 = gyoto.core.matrix4Invert(A)
        self.assertEqual(Am1.shape, (4, 4))
        # C++ signature takes A as const: the input must be left unchanged
        numpy.testing.assert_array_equal(A, A_copy)
        # Reference: numpy's LAPACK-based inverse
        numpy.testing.assert_allclose(Am1, numpy.linalg.inv(A), rtol=rtol, atol=atol)
        # Defining property of the inverse
        I4 = numpy.identity(4)
        numpy.testing.assert_allclose(A @ Am1, I4, rtol=rtol, atol=atol)
        numpy.testing.assert_allclose(Am1 @ A, I4, rtol=rtol, atol=atol)

    def test_identity(self):
        self.checkInverse(numpy.identity(4))

    def test_diagonal(self):
        self.checkInverse(numpy.diag([2., -3., 4., 0.5]))

    def test_permutation(self):
        # Zero diagonal: forces pivoting (row/column swaps)
        P = numpy.zeros((4, 4))
        P[numpy.arange(4), [3, 2, 1, 0]] = 1.
        self.checkInverse(P)

    def test_random_seeded(self):
        # Reproducible: identical matrices at every run
        for _ in range(10):
            self.checkInverse(self.rng.random((4, 4)))

    def test_known_regression(self):
        # This matrix was inverted WRONGLY by the previous implementation
        # (error 0.3 on A*Am1): non-trivial pivoting at several steps
        A = numpy.array([[0., 4., 0., 2.],
                         [-2., -2., 2., 4.],
                         [0., 4., 2., 2.],
                         [1., 3., -3., -4.]])
        self.checkInverse(A)

    def test_ill_conditioned(self):
        # Hilbert 4x4, cond ~ 1.5e4: looser tolerance
        H = numpy.array([[1. / (i + j + 1) for j in range(4)] for i in range(4)])
        self.checkInverse(H, rtol=1e-10, atol=1e-10)

    def test_singular_returns_nan(self):
        # Fixed implementation fills Am1 with NaN for a singular matrix.
        # Binary-exact entries => the zero pivot is exact at any precision.
        # NOTE: like the random tests, this FAILS on the old buggy code.
        A = numpy.array([[1., 1., 1., 1.],
                         [2., 2., 2., 2.],
                         [0., 1., 0., 1.],
                         [1., 0., 1., 0.]])
        Am1 = gyoto.core.matrix4Invert(numpy.ascontiguousarray(A, dtype=numpy.float64))
        self.assertTrue(numpy.isnan(Am1).all())
