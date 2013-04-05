/**
 * \file GyotoUtils.h
 * \brief GYOTO utilities
 *
 *  Various utilities
 */

/*
    Copyright 2011 Thibaut Paumard

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
}

#endif
