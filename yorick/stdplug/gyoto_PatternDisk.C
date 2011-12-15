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
    ao_ = ypush_Astrobj();
    *ao_ = new PatternDisk();
  } else *ypush_Astrobj()=*ao_;

  SmartPointer<PatternDisk> *ao = (SmartPointer<PatternDisk> *)ao_;

  static char const * knames[]={
    "patternvelocity", "readfile",
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
  // Call generic ThinDisk worker

  /* PATTERNVELOCITY */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getPatternVelocity());
    } else
      (*ao)->setPatternVelocity(ygets_d(iarg)) ;
  }

  /* READFILE */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    (*ao)->readFile(ygets_q(iarg));
  }

  ygyoto_ThinDisk_generic_eval(ao_, kiargs+k+1, piargs, rvset, paUsed);
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
