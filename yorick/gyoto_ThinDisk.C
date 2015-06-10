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

#include <Gyoto.h>
#include "ygyoto.h"
#include "ygyoto_private.h"
#include "yapi.h"

using namespace Gyoto;
using namespace Gyoto::Astrobj;

#include <iostream>
using namespace std;

// on_eval worker
void ygyoto_ThinDisk_eval(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>* ao_, int argc) {
  int rvset[1]={0}, paUsed[1]={0};
  if (!ao_) { // Constructor mode
    ao_ = ypush_Astrobj();
    *ao_ = new ThinDisk();
  } else *ypush_Astrobj()=*ao_;

  static char const * knames[]={
    "unit",
    YGYOTO_THINDISK_GENERIC_KW,
    0
  };
  static long kglobs[YGYOTO_THINDISK_GENERIC_KW_N+6];
  int kiargs[YGYOTO_THINDISK_GENERIC_KW_N+5];
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

  char * unit = NULL;

  /* UNIT */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    GYOTO_DEBUG << "get unit" << endl;
    unit = ygets_q(iarg);
  }

  // Call generic ThinDisk worker
  ygyoto_ThinDisk_generic_eval(ao_, kiargs+k+1, piargs, rvset, paUsed, unit);
}


void
ygyoto_ThinDisk_generic_eval(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>*ao_,
			    int *kiargs, int *piargs,
			     int *rvset, int *paUsed, char * unit) {
  SmartPointer<ThinDisk> *ao = (SmartPointer<ThinDisk> *)ao_;
  int k=-1, iarg;
  char const * rmsg="Cannot set return value more than once";

  if (debug())
    for (int i=0; i<YGYOTO_THINDISK_GENERIC_KW_N; ++i)
      cerr << "DEBUG: Astrobj_generic_eval: kiargs[" << i << "]="
	   << kiargs[i] << endl;

  /* INNERRADIUS */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->innerRadius(unit?unit:""));
    } else
      (*ao)->innerRadius(ygets_d(iarg), unit?unit:"") ;
  }

  /* OUTERRADIUS */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->outerRadius(unit?unit:""));
    } else
      (*ao)->outerRadius(ygets_d(iarg), unit?unit:"") ;
  }

  /* THICKNESS */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->thickness(unit?unit:""));
    } else
      (*ao)->thickness(ygets_d(iarg), unit?unit:"") ;
  }

  /* DIR */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_long((*ao)->dir());
    } else
      (*ao)->dir(ygets_l(iarg)) ;
  }

  GYOTO_DEBUG << "calling ygyoto_Astrobj_generic_eval\n";
  ygyoto_Astrobj_generic_eval(ao_, kiargs+k+1, piargs, rvset, paUsed, unit);
  if (debug()) cerr << "DEBUG: out of ThinDisk_generic_eval"<< endl;
}

extern "C" {
  void Y__gyoto_ThinDisk_register_as_Astrobj(){
    ygyoto_Astrobj_register("ThinDisk",&ygyoto_ThinDisk_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_ThinDisk(int argc)
  {
    SmartPointer<Astrobj::Generic> *ao = NULL;
    if (yarg_Astrobj(argc-1)) {
      ao = yget_Astrobj(--argc);
      if ((*ao)->kind().compare("ThinDisk"))
	y_error("Expecting Astrobj of kind Star");
    }
    ygyoto_ThinDisk_eval(ao, argc);
  }

}
