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

#include <cstring>

#include <Gyoto.h>
#include "ygyoto.h"
#include "yapi.h"

using namespace Gyoto;
using namespace Gyoto::Astrobj;

#include <iostream>
using namespace std;

// on_eval worker
void ygyoto_PatternDisk_eval(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>
				*ao_, int argc) {
  int rvset[1]={0}, paUsed[1]={0};
  if (!ao_) { // Constructor mode
    GYOTO_DEBUG << "constructing object\n";
    ao_ = ypush_Astrobj();
    *ao_ = new PatternDisk();
  } else *ypush_Astrobj()=*ao_;

  SmartPointer<PatternDisk> *ao = (SmartPointer<PatternDisk> *)ao_;


  GYOTO_DEBUG << "processing keywords\n";
  static char const * knames[]={
    "fitsread", "patternvelocity", "repeatphi", "nu0", "dnu",
    "phimin", "phimax",
    "copyintensity", "copyopacity", "copyvelocity", "copygridradius",
    "fitswrite",
    YGYOTO_THINDISK_GENERIC_KW,
    0
  };
  static long kglobs[YGYOTO_THINDISK_GENERIC_KW_N+13];
  int kiargs[YGYOTO_THINDISK_GENERIC_KW_N+12];
  int piargs[]={-1,-1,-1,-1};
  
  yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
  
  int iarg=argc, parg=0;
  while (iarg>=1) {
    iarg = yarg_kw(iarg, kglobs, kiargs);
    if (iarg>=1) {
      if (parg<4) piargs[parg++]=iarg--;
      else y_error("gyoto_Astrobj takes at most 4 positional arguments");
    }
  }

  int k=-1;
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";
  // Call generic ThinDisk worker

  /* FITSREAD */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "fitsread=\n";
    iarg+=*rvset;
    (*ao)->fitsRead(ygets_q(iarg));
  }

  /* PATTERNVELOCITY */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getPatternVelocity());
    } else
      (*ao)->setPatternVelocity(ygets_d(iarg)) ;
  }

  /* REPEATPHI */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "repeatphi=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_long((*ao)->repeatPhi());
    } else
      (*ao)->repeatPhi(ygets_l(iarg)) ;
  }

  /* NU0 */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "nu0=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->nu0());
    } else
      (*ao)->nu0(ygets_d(iarg)) ;
  }

  /* DNU */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "dnu=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->dnu());
    } else
      (*ao)->dnu(ygets_d(iarg)) ;
  }

  /* PHIMIN */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "phimin=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->phimin());
    } else
      (*ao)->phimin(ygets_d(iarg)) ;
  }

  /* PHIMAX */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "phimax=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->phimax());
    } else
      (*ao)->phimax(ygets_d(iarg)) ;
  }

  /* INTENSITY */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "copyintensity=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      size_t ddims[3];
      (*ao) -> getIntensityNaxes(ddims);
      long dims[] = {3, ddims[0], ddims[1], ddims[2]};
      double * out = ypush_d(dims);
      memcpy(out, (*ao)->getIntensity(),
	     dims[1]*dims[2]*dims[3]*sizeof(double));
    } else {
      long ntot;
      long dims[Y_DIMSIZE];
      double const * const in = ygeta_d(iarg, &ntot, dims);
      if (dims[0]==0 && ntot && *in==0) (*ao) -> copyIntensity(NULL, 0);
      else if (dims[0]==3) {
	size_t ddims[] = {dims[1], dims[2], dims[3]};
	(*ao)->copyIntensity(in, ddims);
      } else
	y_error("COPYINTENSITY must be nil, 0, or array(double, nnu, nphi, nr");
    }
  }

  /* OPACITY */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "copyopacity=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      size_t ddims[3];
      (*ao) -> getIntensityNaxes(ddims);
      long dims[] = {3, ddims[0], ddims[1], ddims[2]};
      double * out = ypush_d(dims);
      memcpy(out, (*ao)->getOpacity(),
	     dims[1]*dims[2]*dims[3]*sizeof(double));
    } else {
      long ntot;
      long dims[Y_DIMSIZE];
      double const * const in = ygeta_d(iarg, &ntot, dims);
      if (dims[0]==0 && ntot && *in==0) (*ao) -> copyOpacity(NULL, 0);
      else if (dims[0]==3) {
	size_t ddims[] = {dims[1], dims[2], dims[3]};
	(*ao)->copyOpacity(in, ddims);
      } else
	y_error("COPYOPACITY must be nil, 0, or array(double, nnu, nphi, nr");
    }
  }

  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "copyvelocity=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      size_t ddims[3];
      (*ao) -> getIntensityNaxes(ddims);
      long dims[] = {3, 2, ddims[1], ddims[2]};
      double * out = ypush_d(dims);
      memcpy(out, (*ao)->getVelocity(),
	     2*dims[2]*dims[3]*sizeof(double));
    } else {
      long ntot;
      long dims[Y_DIMSIZE];
      double const * const in = ygeta_d(iarg, &ntot, dims);
      if (dims[0]==0 && ntot && *in==0) (*ao) -> copyVelocity(NULL, 0);
      else if (dims[0]==3 && dims[1]==2) {
	size_t ddims[] = {dims[2], dims[3]};
	(*ao)->copyVelocity(in, ddims);
      } else
	y_error("COPYVELOCITY must be nil, 0, or array(double, 2, nphi, nr");
    }
  }

  /* GRIDRADIUS */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "copygridradius=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      double const * const radius = (*ao)->getGridRadius();
      if (radius) {
	size_t ddims[3];
	(*ao) -> getIntensityNaxes(ddims);
	long dims[] = {1, ddims[2]};
	double * out = ypush_d(dims);
	memcpy(out, radius, ddims[2]*sizeof(double));
      } else ypush_long(0);
    } else {
      long ntot;
      long dims[Y_DIMSIZE];
      double const * const in = ygeta_d(iarg, &ntot, dims);
      if (dims[0]==0 && ntot && *in==0)  (*ao) -> copyGridRadius(NULL, 0);
      else if (dims[0]==1) (*ao) -> copyGridRadius(in, ntot);
      else y_error("COPYGRIDRADIUS must be nil, 0, or array(double, nr");
    }
  }

  /* FITSWRITE */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "fitswrite=\n";
    iarg+=*rvset;
    (*ao)->fitsWrite(ygets_q(iarg));
  }

  GYOTO_DEBUG << "calling ygyoto_ThinDisk_generic_eval\n";
  ygyoto_ThinDisk_generic_eval(ao_, kiargs+k+1, piargs, rvset, paUsed);
  GYOTO_DEBUG << "done\n";
}

extern "C" {
  void Y__gyoto_PatternDisk_register_as_Astrobj(){
    ygyoto_Astrobj_register("PatternDisk",&ygyoto_PatternDisk_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_PatternDisk(int argc)
  {
    SmartPointer<Astrobj::Generic> *ao = NULL;
    if (yarg_Astrobj(argc-1)) {
      ao = yget_Astrobj(--argc);
      if ((*ao)->getKind().compare("PatternDisk"))
	y_error("Expecting Astrobj of kind PatternDisk");
    }
    ygyoto_PatternDisk_eval(ao, argc);
  }

}
