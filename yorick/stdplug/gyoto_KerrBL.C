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

#include <GyotoKerrBL.h>
#include "ygyoto.h"
#include "yapi.h"

#include <iostream>
using namespace std;

using namespace Gyoto;
using namespace Gyoto::Metric;

// on_eval worker
void ygyoto_KerrBL_eval(Gyoto::SmartPointer<Gyoto::Metric::Generic> *gg_, int argc) {
  if (debug()) cerr << "DEBUG: in ygyoto_KerrBL_eval()\n";
  int rvset[1]={0}, paUsed[1]={0};
  if (!gg_) { // Constructor mode
    gg_ = ypush_Metric();
    *gg_ = new KerrBL();
  } else *ypush_Metric()=*gg_;

  SmartPointer<KerrBL> *gg = (SmartPointer<KerrBL> *)gg_;
  static char const * knames[]={
    "spin", "makecoord",		\
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
  
  // Process specific keywords
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";
  int k=-1;

  // spin
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      yarg_drop(1);
      ypush_double((*gg)->getSpin());
    } else
      (*gg)->setSpin(ygets_d(iarg)) ;
  }

  if (kiargs[++k]>=0) { // makecoord
    if (debug()) cerr << "DEBUG: In ygyoto_KerrBL_eval(): get_coord" << endl;
    if ((*rvset)++) y_error(rmsg);
    if ((*paUsed)++) y_error(pmsg);
    long dims[]={1,8};
    long ntot=1;
    double * coord_=ygeta_d(kiargs[k],&ntot,0);
    if (ntot<7) y_error("YINIT should have >= 7 elements");
    double * cst_=ygeta_d(piargs[0],&ntot,0);
    if (ntot!=4) y_error("CST should have 4 elements");
    yarg_drop(1);
    double * coord=ypush_d(dims);
#   ifdef GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG_EXPR(gg);
#   endif
    (*gg)->MakeCoord(coord_,cst_,coord);
  }
  

  ygyoto_Metric_generic_eval(gg_, kiargs+k+1, piargs, rvset, paUsed);
  
  if (debug()) cerr << "DEBUG: ygyoto_KerrBL_eval() done\n";
}


extern "C" {
  void Y__gyoto_KerrBL_register_as_Metric(){
    ygyoto_Metric_register("KerrBL",&ygyoto_KerrBL_eval);
  }

  // KERR CLASS
  // Constructor

  void
  Y_gyoto_KerrBL(int argc)
  {
      SmartPointer<Metric::Generic> *gg = NULL;
    try {
      //    char *obj_type=(char*)yget_obj(argc-1,0);
      //    if (obj_type && //!strcmp(obj_type, "gyoto_Metric")) {
      //    if (yget_obj(argc-1,0) && yarg_typeid(argc-1)==Y_OPAQUE) {
      if (yarg_Metric(argc-1)) {
	gg = yget_Metric(--argc);
	if ((*gg)->getKind() != "KerrBL")
	  y_error("Expecting Metric of kind KerrBL");
      }
    } YGYOTO_STD_CATCH;
    ygyoto_KerrBL_eval(gg, argc);
  }

}
