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

#include <GyotoTorus.h>
#include "../ygyoto.h"
#include "yapi.h"

using namespace Gyoto;
using namespace Gyoto::Astrobj;

#include <iostream>
using namespace std;

#define OBJ ao

// on_eval worker
void ygyoto_Torus_eval(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>* ao_, int argc) {
  int rvset[1]={0}, paUsed[1]={0};
  if (!ao_) { // Constructor mode
    ao_ = ypush_Astrobj();
    *ao_ = new Torus();
  } else *ypush_Astrobj()=*ao_;

  SmartPointer<Torus> *ao = (SmartPointer<Torus> *)ao_;
  static char const * knames[]={
    "unit", "largeradius", "smallradius", "spectrum", "opacity",
    YGYOTO_ASTROBJ_GENERIC_KW,
    0
  };
  static long kglobs[YGYOTO_ASTROBJ_GENERIC_KW_N+6];
  int kiargs[YGYOTO_ASTROBJ_GENERIC_KW_N+5];
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
  YGYOTO_WORKER_GETSET_DOUBLE_UNIT(LargeRadius);
  YGYOTO_WORKER_GETSET_DOUBLE_UNIT(SmallRadius);
  YGYOTO_WORKER_GETSET_OBJECT(Spectrum);
# define ypush_Opacity ypush_Spectrum
# define yget_Opacity yget_Spectrum
  YGYOTO_WORKER_GETSET_OBJECT(Opacity);
# undef ypush_Opacity
# undef yget_Opacity
  // Call generic Astrobj worker
  ygyoto_Astrobj_generic_eval(ao_, kiargs+k+1, piargs, rvset, paUsed, unit);
}


extern "C" {
  void Y__gyoto_Torus_register_as_Astrobj(){
    ygyoto_Astrobj_register("Torus",&ygyoto_Torus_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_Torus(int argc)
  {
    SmartPointer<Astrobj::Generic> *ao = NULL;
    if (yarg_Astrobj(argc-1)) {
      ao = yget_Astrobj(--argc);
      if ((*ao)->getKind().compare("Torus"))
	y_error("Expecting Astrobj of kind Star");
    }
    ygyoto_Torus_eval(ao, argc);
  }

}
