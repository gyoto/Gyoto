/*
    Copyright 2013, 2016, 2018 Frederic Vincent, Thibaut Paumard

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
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); GYOTO_ERROR(ermsg); }

#include "GyotoPhoton.h"
#include "GyotoDynamicalDiskBolometric.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"


#include <fitsio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>
#include <limits>
#include <sstream>
#include <dirent.h>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(DynamicalDiskBolometric, "DynamicalDisk with bolometric emission")
GYOTO_PROPERTY_END(DynamicalDiskBolometric, DynamicalDisk::properties)


DynamicalDiskBolometric::DynamicalDiskBolometric() :
  DynamicalDisk()
{
  GYOTO_DEBUG << "DynamicalDiskBolometric Construction" << endl;
}

DynamicalDiskBolometric::DynamicalDiskBolometric(const DynamicalDiskBolometric& o) :
  DynamicalDisk(o)
{
  GYOTO_DEBUG << "DynamicalDiskBolometric Copy" << endl;
}
DynamicalDiskBolometric* DynamicalDiskBolometric::clone() const
{ return new DynamicalDiskBolometric(*this); }

DynamicalDiskBolometric::~DynamicalDiskBolometric() {
  GYOTO_DEBUG << "DynamicalDiskBolometric Destruction" << endl;
}

double DynamicalDiskBolometric::emission(double nu_em, double dsem,
					 state_t const &,
					 double const coord_obj[8]) const{
  GYOTO_ERROR("In DynamicalDiskBolometric::emission: "
	     "not implemented");
  return 0.;
}

double DynamicalDiskBolometric::bolometricEmission(double dsem,
						   state_t const &cph,
						   double const coord_obj[8]) const{
  //  cout << "in bolometricem, will return: " << DynamicalDisk::emission(0.,dsem,coord_obj,coord_obj) << endl;
  return DynamicalDisk::emission(0.,dsem,cph,coord_obj);  
  //never mind about 1st and 3rd elements
}

void DynamicalDiskBolometric::processHitQuantities(Photon* ph, 
						   state_t const & coord_ph_hit,
						   double const * coord_obj_hit, double dt,
						   Properties* data) const {
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
  /*
      NB: freqObs is the observer's frequency chosen in
      Screen::getRayCoord for the actual computation of the geodesic ;
      the physical value of nuobs will be used in spectrum
      computations by resorting to the xml specifications of the user
      (see below) ; this freqObs is used to transform the null
      worldline parameter dlambda (see below)
  */
  double freqObs=ph->freqObs(); // this is a useless quantity, always 1
  double dlambda = dt/coord_ph_hit[4]; //dlambda = dt/tdot
  double ggredm1 = -gg_->ScalarProd(&coord_ph_hit[0],&coord_obj_hit[4],
				    &coord_ph_hit[4]);// / 1.; 
                                       //this is nu_em/nu_obs
  double ggred = 1./ggredm1;           //this is nu_obs/nu_em
  double dsem = dlambda*ggredm1; // *1.
  double inc =0.;
  if (data) {
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "data requested. " 
	      << ", ggredm1=" << ggredm1
	      << ", ggred=" << ggred
	      << endl;
#endif

    if (data->redshift) {
      *data->redshift=ggred;
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->redshift);
#endif
    }
    if (data->time) {
      *data->time=coord_ph_hit[0];
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->time);
#endif
    }
    if (data->impactcoords) {
      if (coord_ph_hit.size()>8) GYOTO_ERROR("ImpactCoord incompatible with parallel transport");
      memcpy(data->impactcoords, coord_obj_hit, 8 * sizeof(double));
      memcpy(data->impactcoords+8, &coord_ph_hit[0], 8 * sizeof(double));
    }
#if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "dlambda = (dt="<< dt << ")/(tdot="<< coord_ph_hit[4]
		<< ") = " << dlambda << ", dsem=" << dsem << endl;
#endif
    if (data->intensity) GYOTO_ERROR("In DynamicalDiskBolometric::process: "
				    "unimplemented");
    else if (data->user4) {
      inc = (bolometricEmission(dsem, coord_ph_hit, coord_obj_hit))
	* (ph -> getTransmission(size_t(-1)))
	* ggred*ggred*ggred*ggred; // I/nu^4 invariant
      *data->user4 += inc;
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->user4);
#endif

    }
    if (data->binspectrum) GYOTO_ERROR("In DynamicalDiskBolometric::process: "
				      "unimplemented");
    if (data->spectrum)  GYOTO_ERROR("In DynamicalDiskBolometric::process: "
				    "unimplemented");
    /* update photon's transmission */
    ph -> transmit(size_t(-1),
		   transmission(freqObs*ggredm1, dsem,coord_ph_hit, coord_obj_hit));
  } else {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "NO data requested!" << endl;
#   endif
  }
}
