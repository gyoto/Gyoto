/*
    Copyright 2014 Frederic Vincent, Thibaut Paumard

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

#include <cstring>

#include <Gyoto.h>
#include "../ygyoto.h"
#include "yapi.h"
#include <GyotoFactory.h>

using namespace Gyoto;
using namespace Gyoto::Astrobj;

#include <iostream>
using namespace std;

// on_eval worker
void ygyoto_DirectionalDisk_eval(SmartPointer<Astrobj::Generic> *OBJ_, int argc) {

  static char const * knames[]={
    "unit",
#ifdef GYOTO_USE_CFITSIO
    "fitsread", 
#endif
    "copyintensity", "copygridfreq", "copygridcosi", "copygridradius",
#ifdef GYOTO_USE_CFITSIO
    "fitswrite",
#endif
    YGYOTO_THINDISK_GENERIC_KW,
    0
  };

#ifdef GYOTO_USE_CFITSIO
#define NKW 8
#else
#define NKW 6
#endif

  YGYOTO_WORKER_INIT(Astrobj, DirectionalDisk, knames,
		     YGYOTO_THINDISK_GENERIC_KW_N+NKW);

  YGYOTO_WORKER_SET_UNIT;
#ifdef GYOTO_USE_CFITSIO
  YGYOTO_WORKER_RUN( fitsRead(ygets_q(iarg)) );
#endif
  //YGYOTO_WORKER_GETSET_DOUBLE2(patternVelocity);
  //YGYOTO_WORKER_GETSET_LONG2(repeatPhi);
  //YGYOTO_WORKER_GETSET_DOUBLE2(nu0);
  //YGYOTO_WORKER_GETSET_DOUBLE2(dnu);
  //YGYOTO_WORKER_GETSET_DOUBLE2(phimin);
  //YGYOTO_WORKER_GETSET_DOUBLE2(phimax);

  /* INTENSITY */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "copyintensity=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      size_t ddims[3];
      (*OBJ) -> getIntensityNaxes(ddims);
      long dims[] = {3, long(ddims[0]), long(ddims[1]), long(ddims[2])};
      double * out = ypush_d(dims);
      memcpy(out, (*OBJ)->getIntensity(),
	     dims[1]*dims[2]*dims[3]*sizeof(double));
    } else {
      long ntot;
      long dims[Y_DIMSIZE];
      double const * const in = ygeta_d(iarg, &ntot, dims);
      if (dims[0]==0 && ntot && *in==0) (*OBJ) -> copyIntensity(NULL, 0);
      else if (dims[0]==3) {
	size_t ddims[] = {size_t(dims[1]), size_t(dims[2]), size_t(dims[3])};
	(*OBJ)->copyIntensity(in, ddims);
      } else
	y_error("COPYINTENSITY must be nil, 0, or array(double, nnu, nphi, nr");
    }
  }

  /* GRIDFREQ */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "copygridfreq=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      double const * const freq = (*OBJ)->getGridFreq();
      if (freq) {
	size_t ddims[3];
	(*OBJ) -> getIntensityNaxes(ddims);
	long dims[] = {1, long(ddims[0])};
	double * out = ypush_d(dims);
	memcpy(out, freq, ddims[0]*sizeof(double));
      } else ypush_long(0);
    } else {
      long ntot;
      long dims[Y_DIMSIZE];
      double const * const in = ygeta_d(iarg, &ntot, dims);
      if (dims[0]==0 && ntot && *in==0)  (*OBJ) -> copyGridFreq(NULL, 0);
      else if (dims[0]==1) (*OBJ) -> copyGridFreq(in, ntot);
      else y_error("COPYGRIDFREQ must be nil, 0, or array(double, nr");
    }
  }

  /* GRIDCOSI */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "copygridcosi=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      double const * const cosi = (*OBJ)->getGridCosi();
      if (cosi) {
	size_t ddims[3];
	(*OBJ) -> getIntensityNaxes(ddims);
	long dims[] = {1, long(ddims[1])};
	double * out = ypush_d(dims);
	memcpy(out, cosi, ddims[1]*sizeof(double));
      } else ypush_long(0);
    } else {
      long ntot;
      long dims[Y_DIMSIZE];
      double const * const in = ygeta_d(iarg, &ntot, dims);
      if (dims[0]==0 && ntot && *in==0)  (*OBJ) -> copyGridCosi(NULL, 0);
      else if (dims[0]==1) (*OBJ) -> copyGridCosi(in, ntot);
      else y_error("COPYGRIDCOSI must be nil, 0, or array(double, nr");
    }
  }

  /* GRIDRADIUS */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "copygridradius=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      double const * const radius = (*OBJ)->getGridRadius();
      if (radius) {
	size_t ddims[3];
	(*OBJ) -> getIntensityNaxes(ddims);
	long dims[] = {1, long(ddims[2])};
	double * out = ypush_d(dims);
	memcpy(out, radius, ddims[2]*sizeof(double));
      } else ypush_long(0);
    } else {
      long ntot;
      long dims[Y_DIMSIZE];
      double const * const in = ygeta_d(iarg, &ntot, dims);
      if (dims[0]==0 && ntot && *in==0)  (*OBJ) -> copyGridRadius(NULL, 0);
      else if (dims[0]==1) (*OBJ) -> copyGridRadius(in, ntot);
      else y_error("COPYGRIDRADIUS must be nil, 0, or array(double, nr");
    }
  }

#ifdef GYOTO_USE_CFITSIO
  YGYOTO_WORKER_RUN( fitsWrite(ygets_q(iarg)) );
#endif
  YGYOTO_WORKER_CALL_GENERIC(ThinDisk);
}

extern "C" {
  void Y__gyoto_DirectionalDisk_register_as_Astrobj(){
    ygyoto_Astrobj_register("DirectionalDisk",&ygyoto_DirectionalDisk_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_DirectionalDisk(int argc)
  {
    YGYOTO_CONSTRUCTOR_INIT2(Astrobj, Astrobj::Generic, DirectionalDisk, astrobj);
    if ((*OBJ)->kind().compare("DirectionalDisk"))
      y_error("Expecting Astrobj of kind DirectionalDisk");
    ygyoto_DirectionalDisk_eval(OBJ, argc);
  }

}
