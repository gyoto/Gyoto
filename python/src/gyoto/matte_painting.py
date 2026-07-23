"""Matte Painting Module for Gyoto.

This module provides functionality to visualize gravitational lensing
effects on painted backgrounds. It includes various painter classes
that can apply geometric distortion effects to images based on
ray-traced photon directions.

The main function `matte_paint` takes a Gyoto Scenery or impact
coordinates and applies a painter function to render the final image
with lensing effects.

Author: Thibaut Paumard

"""

__all__ = ["matte_paint", "rotation", "Painter", "PMode", "Picture", "Panorama"]

import sys
from typing import Optional, Union, Callable, Tuple

import numpy as np
from numpy import pi, cos, sin, atan2

from gyoto import core
from gyoto.core import GYOTO_COORDKIND_SPHERICAL

# Maximum double precision floating point value
DBL_MAX = sys.float_info.max

def matte_paint(
    set: Union[core.Scenery, np.ndarray],
    paint: Callable[[np.ndarray, np.ndarray], np.ndarray],
    coordkind: Optional[int] = None,
    yaw: Optional[float] = None,
    pitch: Optional[float] = None,
    roll: Optional[float] = None,
    origin: str = "upper",
) -> np.ndarray:
    """Visualize lensing effects on a painted background.

    Only geometrical distortion effects are rendered. The function
    ray-traces the Scenery `set`, computes the origin direction of
    each Photon that originates from infinity, and passes these
    directions to the `paint` callable.

    Args:
        set: A Gyoto Scenery object or impact coordinates array as
            computed with sc(,,impactcoords=). If an array, it should
            have shape (NX, NY, 8) or (NX, NY, 16) where the last
            dimension contains coordinate data.
        paint: A callable that performs the actual painting. It will
            be called as: bg = paint(theta, phi) where theta and phi
            are the direction angles (in radians) from which each
            Photon originates, as seen from the center of the
            coordinate system.  The mask parameter is 0 for Photons
            that do not exist (e.g., photons that would fall into a
            black hole).
        coordkind: The coordinate kind (GYOTO_COORDKIND_SPHERICAL or
            GYOTO_COORDKIND_SPHERICAL).  Required if `set` is an array
            rather than a Scenery.

        yaw: Rotation angle about the Z-axis (in radians).
        pitch: Rotation angle about the Y-axis (in radians).
        roll: Rotation angle about the X-axis (in radians).

        origin: Origin convention for coordinate system, either
            "lower" or "upper". Defaults to "upper".

    Returns:
        np.ndarray: The rendered image as a numpy array.

    Raises:
        AssertionError: If origin is not "lower" or "upper".
        ValueError: If coordkind is None and set is not a Scenery.

    """
    assert origin in ["lower", "upper"], 'origin must be either "lower" or "upper"'

    # Adjust rotation angles based on origin convention
    if origin == "upper":
        pitch = -pitch
    else:
        roll = -roll

    # Extract coordinate data from Scenery or use provided array
    if isinstance(set, core.Scenery):
        res = set.Screen.Resolution
        ii = core.Range(1, res, 1)
        jj = core.Range(1, res, 1)
        grid = core.Grid(ii, jj, "\rj = ")
        ipct = np.zeros((res, res, 16), dtype=float)
        aop = core.AstrobjProperties()
        aop.impactcoords = core.array_double.fromnumpy3(ipct)
        aop.offset = res * res
        set.rayTrace(grid, aop)
        coordkind = set.Metric.coordKind()
        data = ipct[:, :, 8:]
    else:
        if coordkind is None:
            raise ValueError(
                'coordkind must be provided if `set` is not a Scenery'
            )
        if set.shape[2] == 8:
            data = set.copy()
        else:
            data = set[:, :, 8:]

    # Replace DBL_MAX (missing values) with NaN
    data[data == DBL_MAX] = np.nan

    # Prepare rotation matrix to transform from compact object
    # coordinates to painter coordinates.
    # Yaw, pitch, and roll are the angles defining the orientation
    # of the painter relative to the metric coordinate system.
    # The final matrix is: YAW @ PITCH @ ROLL.
    R = np.diag(np.ones(3))

    # Apply rotations in sequence: yaw, pitch, roll
    if yaw is not None:
        R = rotation(2, yaw)
    if pitch is not None:
        R @= rotation(1, pitch)
    if roll is not None:
        R @= rotation(0, roll)

    if coordkind == GYOTO_COORDKIND_SPHERICAL:
        # Unpack spherical coordinates
        t, r, theta, phi, tdot, rdot, thetadot, phidot = (
            data[..., i] for i in range(8)
        )

        # Handle stationary photons (tdot == 0)
        ind = (tdot == 0.)
        tdot[ind] = np.nan
        t[ind] = np.nan

        # Compute proper motion components
        taup = 1. / tdot
        rp = rdot * taup
        php = phidot * taup
        thp = thetadot * taup

        # Trigonometric values
        sth = np.sin(theta)
        cth = np.cos(theta)
        sph = np.sin(phi)
        cph = np.cos(phi)

        # Basis vectors in spherical coordinates
        ur = np.array([sth * cph, sth * sph, cth])
        uth = np.array([cth * cph, cth * sph, -sth])
        uph = np.array([-sph, cph, np.zeros_like(t)])

        # Partial derivatives
        durdph = sth * uph
        durdth = uth

        # Time derivative of direction vector
        durdt = php * durdph + thp * durdth

        # Velocity vector in Cartesian coordinates
        v = rp * ur + r * durdt

    else:
        # Unpack Cartesian coordinates
        t, x, y, z, tdot, xdot, ydot, zdot = (
            data[..., i] for i in range(8)
        )

        # Handle stationary photons (tdot == 0)
        ind = (tdot == 0.)
        tdot[ind] = np.nan
        t[ind] = np.nan

        # Compute proper motion components
        taup = 1. / tdot
        xp = xdot * taup
        yp = ydot * taup
        zp = zdot * taup

        # Velocity vector in Cartesian coordinates
        v = np.array([xp, yp, zp])

    # Express velocity in the painter basis
    # Equivalent to:
    # V = numpy.empty_like(v)
    # for i in range(v.shape[1]):
    #     for j in range(v.shape[2]):
    #         V[:, i, j] = v[:, i, j] @ R
    V = np.einsum('kij,kl->lij', v, R)

    # Transform to spherical coordinates in painter basis
    vproj2 = V[0, :, :]**2 + V[1, :, :]**2
    vr = np.sqrt(vproj2 + V[2, :, :]**2)
    vph = atan2(V[1, :, :], V[0, :, :])
    vth = atan2(np.sqrt(vproj2), V[2, :, :])

    # Convert to painter coordinate system
    theta = pi - vth
    phi = vph - pi

    # Ensure phi is within [-pi, pi]
    phi[phi < -pi] += 2. * pi

    # Apply painter function to get background
    bg = paint(theta, phi)

    # Replace NaN values with 0
    bg[np.isnan(bg)] = 0

    return bg

class Painter:
    """Base class for all painter implementations.

    This is an abstract base class that defines the interface for
    painter objects. Subclasses must implement the __call__ method to
    perform the actual painting operation.

    The painter takes theta and phi angles (in radians) and returns
    a background image array.

    """

    pass

class PMode(Painter):
    """Simple pattern painter for testing.

    Generates a sinusoidal pattern based on theta and phi angles.
    Useful for testing and visualization of coordinate
    transformations.

    Attributes:
        ntheta: Frequency multiplier for theta dimension.
        nphi: Frequency multiplier for phi dimension.

    """

    def __init__(self, ntheta: int = 1, nphi: int = 1) -> None:
        """Initialize PMode painter.

        Args:
            ntheta: Frequency multiplier for theta dimension. Defaults
                to 1.
            nphi: Frequency multiplier for phi dimension. Defaults to 1.

        """
        self.ntheta = ntheta
        self.nphi = nphi

    def __call__(self, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
        """Generate pattern based on theta and phi.

        Args:
            theta: Array of theta angles in radians.
            phi: Array of phi angles in radians.

        Returns:
            np.ndarray: Pattern image with values in range [0, 2].
        """
        bg = np.sin(self.ntheta * theta) * np.sin(self.nphi * phi) + 1.
        return bg

class Picture(Painter):
    """Painter that maps an input image using gnomonic projection.

    This painter takes a source image and projects it onto the
    celestial sphere using gnomonic (tangent plane) projection. The
    projection is centered at (phi1, lambda0) with a given field of
    view.

    Attributes:
        img: Source image to be projected.
        fov: Field of view in radians.
        phi1: Central latitude in radians.
        lambda0: Central longitude in radians.
        shape: Shape of the source image.
        ny: Number of rows in the source image.
        nx: Number of columns in the source image.

    """

    def __init__(
        self,
        img: np.ndarray,
        fov: float = 2. * np.arctan2(36., 100.),
        phi1: float = 0.,
        lambda0: float = 0.,
    ) -> None:
        """Initialize Picture painter.

        Args:
            img: Source image as a numpy array. Cannot be None.
            fov: Field of view in radians. Defaults to approximately
                38.94 degrees (2*arctan(36/100)).
            phi1: Central latitude in radians. Defaults to 0.
            lambda0: Central longitude in radians. Defaults to 0.

        Raises:
            AssertionError: If img is None.
        """
        assert img is not None, "img can't be None"
        self.img = img
        self.shape = img.shape
        self.ny = self.shape[0]
        self.nx = self.shape[1]
        self.fov = fov
        self.phi1 = phi1
        self.lambda0 = lambda0

    def __call__(
        self, theta_spherical: np.ndarray, phi_spherical: np.ndarray
    ) -> np.ndarray:
        """Project source image using gnomonic projection.

        Args:
            theta_spherical: Array of theta angles (polar angle from
                positive z-axis) in radians.
            phi_spherical: Array of phi angles (azimuthal angle in
                x-y plane) in radians.

        Returns:
            np.ndarray: Projected image with same shape as input angles
                and same number of channels as source image.
        """
        # Gnomonic projection:
        # http://mathworld.wolfram.com/GnomonicProjection.html
        #
        # Note: We use 'lamda' and 'lamda0' instead of 'lambda' as
        # lambda is a reserved word in Python.

        dmatte = self.shape
        ndims = len(dmatte)
        nx = self.nx
        ny = self.ny
        phi1 = self.phi1
        lamda0 = self.lambda0
        fov = self.fov
        img = self.img

        # Convert spherical coordinates to latitude/longitude
        phi = pi / 2. - theta_spherical  # latitude
        lamda = -phi_spherical  # longitude

        # Determine field of view in each dimension
        if nx > ny:
            xfov = fov
            yfov = ny * fov / nx
        else:
            yfov = fov
            xfov = nx * fov / ny

        # Gnomonic projection formulas
        cos_c = (
            sin(phi1) * sin(phi) + cos(phi1) * cos(phi) * cos(lamda - lamda0)
        )
        x = cos(phi) * sin(lamda - lamda0) / cos_c
        y = (
            (cos(phi1) * sin(phi) - sin(phi1) * cos(phi) * cos(lamda - lamda0))
            / cos_c
        )

        # Compute corner coordinates for scaling
        phi_corners = phi1 + yfov * 0.5 * np.array([-1., 1.])
        lamda_corners = lamda0 + xfov * 0.5 * np.array([-1., 1.])
        cos_c_corners = (
            sin(phi1) * sin(phi_corners) +
            cos(phi1) * cos(phi_corners) * cos(lamda_corners - lamda0)
        )
        x_corners = (
            cos(phi_corners) * sin(lamda_corners - lamda0) / cos_c_corners
        )
        y_corners = (
            (cos(phi1) * sin(phi_corners) -
             sin(phi1) * cos(phi_corners) * cos(lamda_corners - lamda0))
            / cos_c_corners
        )

        # Center of the image
        i0 = nx * 0.5
        j0 = ny * 0.5

        # Scale factor
        scale = nx / (x_corners[1] - x_corners[0])

        # Convert to pixel coordinates
        xp = i0 + x * scale
        yp = j0 + y * scale

        # Clamp coordinates to image bounds
        xp[xp < 0] = 0
        yp[yp < 0] = 0
        xp[xp >= nx - 1] = nx - 1
        yp[yp >= ny - 1] = ny - 1

        dd = phi.shape

        # Initialize output array
        if ndims == 3:
            bg = np.zeros((dd[0], dd[1], dmatte[2]), dtype=img.dtype)
        elif ndims == 2:
            bg = np.zeros((dd[0], dd[1]), like=img, dtype=img.dtype)
        else:
            raise ValueError('phi_spherical should have 2 or 3 dimensions')

        # Sample from source image
        for i in range(dd[0]):
            for j in range(dd[1]):
                if not np.isnan(xp[i, j]) and not np.isnan(yp[i, j]):
                    bg[i, j, ...] = img[round(yp[i, j]), round(xp[i, j])]

        return bg

class Panorama(Painter):
    """Paints a panoramic image using spherical projection.

    This painter takes a source image representing a full or partial
    spherical panorama and samples from it based on theta and phi
    angles.  The source image is assumed to be in equirectangular
    projection.

    Attributes:
        img: Source panoramic image.
        phi_fov: Field of view in phi (longitude) direction in radians.
        theta_fov: Field of view in theta (latitude) direction in
            radians.

    """

    def __init__(
        self,
        img: np.ndarray,
        phi_fov: float = 2. * pi,
        theta_fov: float = pi,
    ) -> None:
        """Initialize Panorama painter.

        Args:
            img: Source panoramic image as a numpy array. Cannot be
                None.
            phi_fov: Field of view in phi (longitude) direction in
                radians.  Defaults to 2*pi (full circle).
            theta_fov: Field of view in theta (latitude) direction in
                radians. Defaults to pi (180 degrees).

        Raises:
            AssertionError: If img is None.

        """
        assert img is not None, "img can't be None"
        self.img = img
        self.phi_fov = phi_fov
        self.theta_fov = theta_fov

    def __call__(self, theta: np.ndarray, phi: np.ndarray) -> np.ndarray:
        """Sample from panoramic image based on theta and phi angles.

        Args:
            theta: Array of theta angles (polar angle from positive
                z-axis) in radians.
            phi: Array of phi angles (azimuthal angle in x-y plane) in
                radians.

        Returns:
            np.ndarray: Sampled image with same shape as input angles
                and same number of channels as source image.

        """
        img = self.img
        phi_fov = self.phi_fov
        theta_fov = self.theta_fov

        dmatte = img.shape
        ndims = len(dmatte)
        nphi = dmatte[1]
        ntheta = dmatte[0]

        # Scaling factors (pixels per radian)
        # Negative signs because theta increases from top to bottom
        # and phi from right to left, as seen from inside the sphere.
        phi_scale = -nphi / phi_fov
        theta_scale = -ntheta / theta_fov

        # Center of the image (+0.5 for proper rounding)
        i0 = nphi * 0.5
        j0 = ntheta * 0.5

        # Convert angles to pixel coordinates
        x = np.round(i0 + phi * phi_scale)
        y = np.round(j0 + (theta - 0.5 * pi) * theta_scale)

        # Clamp coordinates to image bounds
        x[x < 0] = 0
        y[y < 0] = 0
        x[x >= nphi - 1] = nphi - 1
        y[y >= ntheta - 1] = ntheta - 1

        dd = phi.shape

        # Initialize output array
        if ndims == 3:
            bg = np.zeros((dd[0], dd[1], dmatte[2]), dtype=img.dtype)
        elif ndims == 2:
            bg = np.zeros((dd[0], dd[1]), like=img, dtype=img.dtype)
        else:
            raise ValueError('phi should have 2 or 3 dimensions')

        # Sample from source image
        for i in range(dd[0]):
            for j in range(dd[1]):
                if not np.isnan(x[i, j]) and not np.isnan(y[i, j]):
                    bg[i, j, ...] = img[round(y[i, j]), round(x[i, j])]

        return bg;

def rotation(axis: int, angle: float) -> np.ndarray:
    """Generate a 3D rotation matrix.

    Creates a rotation matrix about the specified axis (Ox=0, Oy=1,
    Oz=2) by the given angle.

    If vector V holds the coordinates of a point in the original
    coordinate system, then:
        - v @ R gives the coordinates of the same point in the rotated
          coordinate system.
        - R @ v gives the coordinates of the rotated vector in the
          original coordinate system.

    If S and T are two rotations, the combined rotation obtained by
    first applying S, then T, is given by R = T @ S.

    Args:
        axis: The axis of rotation (0 for x-axis, 1 for y-axis, 2 for
            z-axis).
        angle: The rotation angle in radians.

    Returns:
        np.ndarray: A 3x3 rotation matrix.

    """
    R = np.diag(np.ones(3))

    ca = np.cos(angle)
    sa = np.sin(angle)

    # 2D rotation matrix
    R2 = np.array([[ca, -sa], [sa, ca]])

    # Determine the plane of rotation (the two axes perpendicular to
    # the rotation axis)
    plane = (axis + 1 + np.array([0, 1])) % 3

    # Apply the 2D rotation to the plane
    R[np.ix_(plane, plane)] = R2

    return R
