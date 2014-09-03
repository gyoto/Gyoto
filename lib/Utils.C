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

#include "GyotoUtils.h"
#include "GyotoError.h"
#include "GyotoPhoton.h"
#include <cmath>

using namespace Gyoto;
using namespace std;

static int gyoto_debug=GYOTO_DEFAULT_DEBUG_MODE;
#if GYOTO_DEFAULT_DEBUG_MODE
static int gyoto_verbosity=GYOTO_DEBUG_VERBOSITY;
static int gyoto_prev_verbosity=GYOTO_DEBUG_VERBOSITY;
#else
static int gyoto_verbosity=GYOTO_DEFAULT_VERBOSITY;
static int gyoto_prev_verbosity=GYOTO_DEBUG_VERBOSITY;
#endif


void Gyoto::debug(int mode) {
  if (mode != gyoto_debug) {
    if (mode) {
      gyoto_prev_verbosity=verbose();
      verbose(GYOTO_DEBUG_VERBOSITY);
    } else {
      verbose(gyoto_prev_verbosity);
    }
    gyoto_debug=mode;
  } 
}
int Gyoto::debug() { return gyoto_debug; }

void Gyoto::verbose(int mode) { gyoto_verbosity=mode; }
int Gyoto::verbose() { return gyoto_verbosity; }

void Gyoto::convert(double * const x, const size_t nelem, const double mass_sun, const double distance_kpc, const string unit) {
  /// Convert lengths
  
  double distance = distance_kpc*GYOTO_KPC;  // m
  double fact   = mass_sun * GYOTO_SUN_MASS * GYOTO_G_OVER_C_SQUARE; // m
  size_t i =0;


  if (!unit.compare("geometrical"))     return ;
  else if (!unit.compare("m"))          ;
  else if (!unit.compare("km"))         fact *=  1e-3 ;
  else if (!unit.compare("sun radius")) fact *=  1.      / GYOTO_SUN_RADIUS;
  else if (!unit.compare("rad"))        fact *=  1.      / (distance);
  else if (!unit.compare("degree"))     fact *=  180.    / (distance*M_PI);
  else if (!unit.compare("arcmin"))     fact *=  1.08e4  / (distance*M_PI);
  else if (!unit.compare("arcsec"))     fact *=  6.48e5  / (distance*M_PI);
  else if (!unit.compare("mas"))        fact *=  6.48e8  / (distance*M_PI);
  else if (!unit.compare("uas"))        fact *=  6.48e11 / (distance*M_PI);
  else throwError("Unknown unit.");

  for (i=0; i<nelem; ++i) x[i] *= fact ;

}

double Gyoto::atof(const char * str)
{
  GYOTO_DEBUG << "Gyoto::atof(\"" << str << "\")";
  ptrdiff_t offset=0;
  while (isspace(str[offset])) ++offset;
  bool positive=true;
  double retval=0.;
  if (str[offset] == '-') {
    positive=false;
    ++offset;
  }
  if (str[offset++]=='D' && str[offset++]=='B' && str[offset++]=='L' &&
      str[offset++]=='_' && str[offset++]=='M') {
    if (str[offset]=='A' && str[offset+1]=='X') {
      if (positive) retval = DBL_MAX;
      else retval = -DBL_MAX;
    } else if (str[offset]=='I' && str[offset+1]=='N') {
      if (positive) retval = DBL_MIN;
      else retval = -DBL_MIN;
    } else throwError("unrecognize double representation");
  } else retval = std::atof(str);

  GYOTO_DEBUG << "==" << retval << endl;

  return retval;
}
