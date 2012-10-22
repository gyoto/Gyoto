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
#include "../ygyoto.h"
#include "yapi.h"

using namespace Gyoto;
using namespace Gyoto::Astrobj;

#include <iostream>
using namespace std;

// on_eval worker
void ygyoto_PageThorneDisk_eval(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>
				*ao_, int argc) {
  int rvset[1]={0}, paUsed[1]={0};
  if (!ao_) { // Constructor mode
    ao_ = ypush_Astrobj();
    *ao_ = new PageThorneDisk();
  } else *ypush_Astrobj()=*ao_;

  SmartPointer<PageThorneDisk> *ao = (SmartPointer<PageThorneDisk> *)ao_;

  static char const * knames[]={
    "unit", "updatespin",
    YGYOTO_THINDISK_GENERIC_KW,
    0
  };
  static long kglobs[YGYOTO_THINDISK_GENERIC_KW_N+3];
  int kiargs[YGYOTO_THINDISK_GENERIC_KW_N+2];
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

  // Call generic ThinDisk worker

  /* UNIT */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    GYOTO_DEBUG << "get unit" << endl;
    unit = ygets_q(iarg);
  }

  /* UPDATESPIN */
  if ((iarg=kiargs[++k])>=0) (*ao) -> updateSpin();

  ygyoto_ThinDisk_generic_eval(ao_, kiargs+k+1, piargs, rvset, paUsed, unit);
}

extern "C" {
  void Y__gyoto_PageThorneDisk_register_as_Astrobj(){
    ygyoto_Astrobj_register("PageThorneDisk",&ygyoto_PageThorneDisk_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_PageThorneDisk(int argc)
  {
    SmartPointer<Astrobj::Generic> *ao = NULL;
    if (yarg_Astrobj(argc-1)) {
      ao = yget_Astrobj(--argc);
      if ((*ao)->getKind().compare("PageThorneDisk"))
	y_error("Expecting Astrobj of kind PageThorneDisk");
    }
    ygyoto_PageThorneDisk_eval(ao, argc);
  }

}
