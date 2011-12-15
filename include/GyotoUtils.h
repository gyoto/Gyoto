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

#include <string>
#include <GyotoDefs.h>

namespace Gyoto {
  class Photon;
  void debug(int mode); ///< mode=1 for debug output, 0 for no output
  int debug(); ///< return >=1 if in debug mode, else 0
  void verbose(int mode); ///< mode=1 for debug output, 0 for no output
  int verbose(); ///< return >=1 if in debug mode, else 0
  void convert(double * const x, const std::size_t nelem,
	       const double mass_sun, const double distance_kpc,
	       const std::string unit); /// Convert lengths
}

#define GYOTO_DEBUG(more) if (debug()) cerr << "DEBUG: "	   \
					    << __PRETTY_FUNCTION__ \
					    << ": " more	   \
					    << endl

#endif
