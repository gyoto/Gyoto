/*
    Copyright 2013 Thibaut Paumard

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

#include <GyotoComplexSpectrometer.h>
#include "ygyoto.h"
#include "yapi.h"

#include <iostream>
using namespace std;

using namespace Gyoto;
using namespace Gyoto::Spectrometer;

namespace YGyoto {
  ygyoto_Spectrometer_eval_worker_t SpCplxEval;
}

// on_eval worker
void YGyoto::SpCplxEval(Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> *sp_, int argc) {
  GYOTO_DEBUG << endl;
  int rvset[1]={0}, paUsed[1]={0};
  if (!sp_) { // Constructor mode
    sp_ = ypush_Spectrometer();
    *sp_ = new Complex();
  } else *ypush_Spectrometer()=*sp_;

  SmartPointer<Complex> *sp = (SmartPointer<Complex> *)sp_;
  static char const * knames[]={
    "unit",
    "append", "remove",	"cardinal",		\
    YGYOTO_SPECTROMETER_GENERIC_KW,
    0
  };
  static long kglobs[YGYOTO_SPECTROMETER_GENERIC_KW_N+5];
  int kiargs[YGYOTO_SPECTROMETER_GENERIC_KW_N+4];
  int piargs[]={-1,-1,-1,-1};
  
  yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
  
  int iarg=argc, parg=0;
  while (iarg>=1) {
    iarg = yarg_kw(iarg, kglobs, kiargs);
    if (iarg>=1) {
      if (parg<4) piargs[parg++]=iarg--;
      else y_error("gyoto_Spectrometer takes at most 4 positional arguments");
    }
  }
  
  // Process specific keywords
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";
  char * unit=NULL;
  int k=-1;

  /* UNIT */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    GYOTO_DEBUG << "get unit" << endl;
    unit = ygets_q(iarg);
  }

  /* append */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    GYOTO_DEBUG << "append subspectro" << endl;
    (*sp)->append(*yget_Spectrometer(iarg));
  }

  /* remove */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    GYOTO_DEBUG << "reomve subspectro" << endl;
    (*sp)->remove(ygets_l(iarg));
  }

  /* cardinal */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if ((*rvset)++) y_error(rmsg);
    GYOTO_DEBUG << "get cardinal" << endl;
    ypush_long((*sp)->getCardinal());
  }
  
  ygyoto_Spectrometer_generic_eval(sp_, kiargs+k+1, piargs, rvset, paUsed, unit);
  
  if (!(*paUsed) && (iarg=piargs[0])>=0 && !yarg_nil(iarg)) {
    iarg+=*rvset;
    if ((*rvset)++) y_error(rmsg);
    ++(*paUsed);
    size_t cardinal=(*sp) -> getCardinal();
    long item = ygets_l(iarg);
    if (item > cardinal) y_error("index overreach array bounds");
    if (item <= 0) {
      item += cardinal;
      if (item <= 0) y_error("index overreach array bounds");
    }
    *ypush_Spectrometer() = (*(*sp))[item-1];
  }

  GYOTO_DEBUG << "done\n";
}


extern "C" {
  void Y__gyoto_SpCplx_register_as_Spectrometer(int argc){
    ygyoto_Spectrometer_register(Complex::Kind,&YGyoto::SpCplxEval);
  }

  // KERR CLASS
  // Constructor

  void
  Y_gyoto_SpectroComplex(int argc)
  {
      SmartPointer<Spectrometer::Generic> *sp = NULL;
    try {
      //    char *obj_type=(char*)yget_obj(argc-1,0);
      //    if (obj_type && //!strcmp(obj_type, "gyoto_Spectrometer")) {
      //    if (yget_obj(argc-1,0) && yarg_typeid(argc-1)==Y_OPAQUE) {
      if (yarg_Spectrometer(argc-1)) {
	sp = yget_Spectrometer(--argc);
	if ((*sp)->getKind() != Complex::Kind)
	  y_error("Expecting Spectrometer of kind Complex");
      }
    } YGYOTO_STD_CATCH;
    YGyoto::SpCplxEval(sp, argc);
  }

}
