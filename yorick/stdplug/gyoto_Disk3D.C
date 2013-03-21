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
#include "../ygyoto.h"
#include "yapi.h"

using namespace Gyoto;
using namespace Gyoto::Astrobj;

#include <iostream>
using namespace std;

#define OBJ ao

// on_eval worker
void ygyoto_Disk3D_eval(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>
				*ao_, int argc) {
  int rvset[1]={0}, paUsed[1]={0};
  if (!ao_) { // Constructor mode
    GYOTO_DEBUG << "constructing object\n";
    ao_ = ypush_Astrobj();
    *ao_ = new Disk3D();
  } else *ypush_Astrobj()=*ao_;

  SmartPointer<Disk3D> *ao = (SmartPointer<Disk3D> *)ao_;


  GYOTO_DEBUG << "processing keywords\n";
  static char const * knames[]={
    "unit",
    "fitsread", "repeatphi", "nu0", "dnu",
    "rin", "rout", "zmin", "zmax",
    "phimin", "phimax",
    "copyemissquant", "copyvelocity",
    "fitswrite",
    YGYOTO_ASTROBJ_GENERIC_KW,
    0
  };
  static long kglobs[YGYOTO_ASTROBJ_GENERIC_KW_N+14];
  int kiargs[YGYOTO_ASTROBJ_GENERIC_KW_N+13];
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
  char * unit=NULL;

  YGYOTO_WORKER_SET_UNIT;
  YGYOTO_WORKER_RUN(fitsRead, ygets_q(iarg));
  YGYOTO_WORKER_GETSET_LONG2(repeatPhi);
  YGYOTO_WORKER_GETSET_DOUBLE2(nu0);
  YGYOTO_WORKER_GETSET_DOUBLE2(dnu);
  YGYOTO_WORKER_GETSET_DOUBLE2(rin);
  YGYOTO_WORKER_GETSET_DOUBLE2(rout);
  YGYOTO_WORKER_GETSET_DOUBLE2(zmin);
  YGYOTO_WORKER_GETSET_DOUBLE2(zmax);
  YGYOTO_WORKER_GETSET_DOUBLE2(phimin);
  YGYOTO_WORKER_GETSET_DOUBLE2(phimax);

  /* EMISSQUANT */
  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "copyemissquant=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      size_t ddims[4];
      (*ao) -> getEmissquantNaxes(ddims);
      long dims[] = {4, ddims[0], ddims[1], ddims[2], ddims[3]};
      double * out = ypush_d(dims);
      memcpy(out, (*ao)->getEmissquant(),
	     dims[1]*dims[2]*dims[3]*dims[4]*sizeof(double));
    } else {
      long ntot;
      long dims[Y_DIMSIZE];
      double const * const in = ygeta_d(iarg, &ntot, dims);
      if (dims[0]==0 && ntot && *in==0) (*ao) -> copyEmissquant(NULL, 0);
      else if (dims[0]==4) {
	size_t ddims[] = {dims[1], dims[2], dims[3], dims[4]};
	(*ao)->copyEmissquant(in, ddims);
      } else
	y_error("COPYEMISSQUANT must be nil, 0, or array(double, nnu, nphi, nz, nr");
    }
  }

  /* VELOCITY */

  if ((iarg=kiargs[++k])>=0) {
    GYOTO_DEBUG << "copyvelocity=\n";
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      size_t ddims[4];
      (*ao) -> getEmissquantNaxes(ddims);
      long dims[] = {4, 3, ddims[1], ddims[2], ddims[3]};
      double * out = ypush_d(dims);
      memcpy(out, (*ao)->getVelocity(),
	     3*dims[2]*dims[3]*dims[4]*sizeof(double));
    } else {
      long ntot;
      long dims[Y_DIMSIZE];
      double const * const in = ygeta_d(iarg, &ntot, dims);
      if (dims[0]==0 && ntot && *in==0) (*ao) -> copyVelocity(NULL, 0);
      else if (dims[0]==4 && dims[1]==3) {
	size_t ddims[] = {dims[2], dims[3], dims[4]};
	(*ao)->copyVelocity(in, ddims);
      } else
	y_error("COPYVELOCITY must be nil, 0, or array(double, 3, nphi, nz, nr");
    }
  }

  /* FITSWRITE */
  YGYOTO_WORKER_RUN(fitsWrite, ygets_q(iarg));

  GYOTO_DEBUG << "calling ygyoto_Astrobj_generic_eval\n";
  ygyoto_Astrobj_generic_eval(ao_, kiargs+k+1, piargs, rvset, paUsed, unit);
  GYOTO_DEBUG << "done\n";
}

extern "C" {
  void Y__gyoto_Disk3D_register_as_Astrobj(){
    ygyoto_Astrobj_register("Disk3D",&ygyoto_Disk3D_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_Disk3D(int argc)
  {
    SmartPointer<Astrobj::Generic> *ao = NULL;
    if (yarg_Astrobj(argc-1)) {
      ao = yget_Astrobj(--argc);
      if ((*ao)->getKind().compare("Disk3D"))
	y_error("Expecting Astrobj of kind Disk3D");
    }
    ygyoto_Disk3D_eval(ao, argc);
  }

}
