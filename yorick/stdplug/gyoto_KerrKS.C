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

#include <GyotoKerrKS.h>
#include "ygyoto.h"
#include "yapi.h"

#include <iostream>
using namespace std;

using namespace Gyoto;
using namespace Gyoto::Metric;

// on_eval worker
void ygyoto_KerrKS_eval(Gyoto::SmartPointer<Gyoto::Metric::Generic> *gg_, int argc) {
  int rvset[1]={0}, paUsed[1]={0};
  if (!gg_) { // Constructor mode
    gg_ = ypush_Metric();
    *gg_ = new KerrKS();
  } else  *ypush_Metric()=*gg_;

  SmartPointer<KerrKS> *gg = (SmartPointer<KerrKS> *)gg_;
  static char const * knames[]={
    "unit", "spin",
    YGYOTO_METRIC_GENERIC_KW,
    0
  };
  static long kglobs[YGYOTO_METRIC_GENERIC_KW_N+4];
  int kiargs[YGYOTO_METRIC_GENERIC_KW_N+3];
  int piargs[]={-1,-1,-1,-1};
  
  yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
  
  int iarg=argc, parg=0;
  while (iarg>=1) {
    iarg = yarg_kw(iarg, kglobs, kiargs);
    if (iarg>=1) {
      if (parg<4) piargs[parg++]=iarg--;
      else y_error("gyoto_Metric takes at most 4 positional arguments");
    }
  }
  
  // Process specific GET keywords
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

  // SPIN
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      yarg_drop(1);
      ypush_double((*gg)->getSpin());
    } else
      (*gg)->setSpin(ygets_d(iarg)) ;
  }

  ygyoto_Metric_generic_eval(gg_, kiargs+k+1, piargs, rvset, paUsed, unit);
  
}


extern "C" {
  void Y__gyoto_KerrKS_register_as_Metric(){
    ygyoto_Metric_register("KerrKS",&ygyoto_KerrKS_eval);
  }

  // KERR CLASS
  // Constructor

  void
  Y_gyoto_KerrKS(int argc)
  {
    SmartPointer<Metric::Generic> *gg = NULL;
    try {
      if (yarg_Metric(argc-1)) {
	gg = yget_Metric(--argc);
	if ((*gg)->getKind() != "KerrKS")
	  y_error("Expecting Metric of kind KerrKS");
      }
    } YGYOTO_STD_CATCH;
    ygyoto_KerrKS_eval(gg, argc);
  }

}
