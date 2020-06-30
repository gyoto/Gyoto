/**
 * \file GyotoUtils.h
 * \brief GYOTO utilities
 *
 *  Various utilities
 */

/*
    Copyright 2011, 2016 Thibaut Paumard

    This file is part of Gyoto.

    Gyoto is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Gyoto is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef __GyotoUtils_H_ 
#define __GyotoUtils_H_ 

#include "GyotoDefs.h"

#include <string>
#include <vector>

namespace Gyoto {
  /// Set debug mode
  /**
   * \param mode 1 to turn on debug mode, 0 to turn it off.
   */
  void debug(int mode);
  
  /// Get debug mode
  /**
   * \return >=1 if debug mode is on, else 0.
   */
  int debug();

  /// Set verbosity level
  /**
   * See standard verbosity levels defined in GyotoDefs.h:
   * 
   * - GYOTO_DEFAULT_DEBUG_MODE
   * - GYOTO_QUIET_VERBOSITY
   * - GYOTO_SEVERE_VERBOSITY
   * - GYOTO_WARNING_VERBOSITY
   * - GYOTO_DEFAULT_VERBOSITY
   * - GYOTO_INFO_VERBOSITY
   * - GYOTO_DEBUG_VERBOSITY
   */
  void verbose(int mode);
  
  /// Get verbosity level
  /**
   * See verbose(int mode).
   */
  int verbose();
  
  /// Convert lengths (deprecated)
  /**
   * \deprecated Will be removed once it is not used anymore in Gyoto
   * per se. Prefer Gyoto::Units framework.
   *
   * \param[in,out] x Lengths to convert, in geometrical units on
   * input, in specified unit on output.
   * \param[in] nelem Size of x array.
   * \param[in] mass_sun Black-hole mass in Solar masses.
   * \param[in] distance_kpc Distance from observer in kiloparsecs.
   * \param[in] unit One of "geometrical", "m", "km", "sun radius",
   * "rad", "degree", "arcmin", "arcsec", "mas", "uas".
   */
  void convert(double * const x, const std::size_t nelem,
	       const double mass_sun, const double distance_kpc,
	       const std::string unit);
  
  /// Interpret C string as double
  /**
   * Wrapper around std::atof() that also interprets DBL_MIN, DBL_MAX,
   * -DBL_MIN and -DBL_MAX.
   *
   * If str starts with "(-)DBL_M" and is not one of the four special
   * values, then an error is thrown.
   *
   * \param[in] str C string to interpret
   * \return  double valu represented by str.
   */
  double atof(const char * str);

  /// Print help on class
  /**
   * \param[in] class_name e.g. "Gyoto::Screen", "Gyoto::Astrobj::Torus".
   */
  void help(std::string class_name);

  /// Split string
  std::vector<std::string> split(std::string const &src, std::string const &delim);

  /// Bessel function computation
  /*
    boost Bessel function as implemented in
    #include <boost/math/special_functions/bessel.hpp> 
    are 50% longer than the following, while the following
    give results accurate at ~1e-6 which is more than enough.
   */
  double bessi0(double xx);///< Modified Bessel function I<SUB>0</SUB>
  double bessi1(double xx);///< Modified Bessel function I<SUB>1</SUB>
  double bessk0(double xx);///< Modified Bessel function K<SUB>0</SUB>
  double bessk1(double xx);///< Modified Bessel function K<SUB>1</SUB>
  double bessk(int nn, double xx);///< Modified Bessel function

  double hypergeom (double kappaIndex, double thetae); ///< Gauss hypergeometric 2F1 term for kappa-distribution synchrotron

  /// Tranform from Cartesian 3-position to spherical 3-position
  void cartesianToSpherical(double const cpos[3], double spos[3]);
  /// Tranform from spherical 3-position to Cartesian 3-position
  void sphericalToCartesian(double const spos[3], double cpos[3]);

  /// Invert 4x4 matrix
  /**
   * \param[in] IN_ARRAY2 the 4×4 matrix to invert
   * \param[out] ARGOUT_ARRAY2 the invert matrix of IN_ARRAY2
   */
  // Keep argument names for swig!
  void matrix4Invert(double ARGOUT_ARRAY2[4][4], double const IN_ARRAY2[4][4]);

  /// Invert 4x4 circular spacetime metric 
  /**
   * A circular spacetime metric (in the right coordinate system)
   *   - is symmetrical (like all metric matrices);
   *   - has only 6 non-zero element: the diagonal and the corners.
   *
   * \param[in] IN_ARRAY2 the 4×4 matrix to invert
   * \param[out] ARGOUT_ARRAY2 the invert matrix of IN_ARRAY2
   */
  // Keep argument names for swig!
  void matrix4CircularInvert(double ARGOUT_ARRAY2[4][4], double const IN_ARRAY2[4][4]);
}

#endif
