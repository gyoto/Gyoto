/*
    Copyright (c) 2012 Thibaut Paumard

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

#include "../ygyoto.h"
#include "yapi.h"
#include <GyotoPolishDoughnut.h>
#ifdef GYOTO_USE_XERCES
#include <GyotoFactory.h>
#endif

#include <iostream>
using namespace std;

using namespace Gyoto;
using namespace Gyoto::Astrobj;

// on_eval worker
void ygyoto_PolishDoughnut_eval(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>* ao_, int argc) {
  int rvset[1]={0}, paUsed[1]={0}, builder=0;
  if (!ao_) { // Constructor mode
    ao_ = ypush_Astrobj();
    builder=1;
  } else *ypush_Astrobj()=*ao_;
  static char const * knames[]={
    "unit",
    "lambda", "tempratio", "centraldensity", "centraltempovervirial", "beta",
    "spectraloversampling",
    "l0", "Wsurface", "Wcentre", "rcusp", "rcentre",
    YGYOTO_ASTROBJ_GENERIC_KW,
    0
  };
  static long kglobs[YGYOTO_ASTROBJ_GENERIC_KW_N+13];
  int kiargs[YGYOTO_ASTROBJ_GENERIC_KW_N+12];
  int piargs[]={-1,-1,-1,-1};
  
  yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
  
  int iarg=argc, parg=0;
  while (iarg>=1) {
    iarg = yarg_kw(iarg, kglobs, kiargs);
    if (iarg>=1) {
      if (parg<4) piargs[parg++]=iarg--;
      else y_error("gyoto_PolishDoughnut takes at most 4 positional arguments");
    }
  }

  if (debug())
    for (int i=0; i<YGYOTO_ASTROBJ_GENERIC_KW_N+1; ++i)
      cerr << "DEBUG gyoto_PolishDoughnut_eval: kiargs[" << i << "]="
	   << kiargs[i] << endl;


  // if rvset==1, constructor mode:
  if (builder) {
    if (yarg_string(piargs[0])) {
#ifdef GYOTO_USE_XERCES
      try { *ao_ = Factory(ygets_q(piargs[0])).getAstrobj(); } 
      YGYOTO_STD_CATCH;
      *paUsed=1;
#else
      y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
    } else *ao_ = new PolishDoughnut();

  }
  SmartPointer<PolishDoughnut> *ao = (SmartPointer<PolishDoughnut> *)ao_;

  int k=-1;
  //// MEMBERS ////
  // "lambda", "tempratio", "centraldensity", "centraltempovervirial", "beta",
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";
  char * unit=NULL;

  /* UNIT */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    GYOTO_DEBUG << "get unit" << endl;
    unit = ygets_q(iarg);
  }

  if ((iarg=kiargs[++k])>=0) { // lambda
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getLambda());
    } else (*ao)->setLambda(ygets_d(iarg));
  }

  if ((iarg=kiargs[++k])>=0) { // tempratio
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getTemperatureRatio());
    } else (*ao)->setTemperatureRatio(ygets_d(iarg));
  }

  if ((iarg=kiargs[++k])>=0) { // centraldensity
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getCentralDensity(unit?unit:""));
    } else (*ao)->setCentralDensity(ygets_d(iarg), unit?unit:"");
  }

  if ((iarg=kiargs[++k])>=0) { // centraltempovervirial
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getCentralTempOverVirial());
    } else (*ao)->setCentralTempOverVirial(ygets_d(iarg));
  }

  if ((iarg=kiargs[++k])>=0) { // beta
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getBeta());
    } else (*ao)->setBeta(ygets_d(iarg));
  }

  if ((iarg=kiargs[++k])>=0) { // spectraloversampling
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_long((*ao)->getSpectralOversampling());
    } else (*ao)->setSpectralOversampling(ygets_l(iarg));
  }

  //    "l0", "Wsurface", "Wcentre", "rcusp", "rcentre",
  if ((iarg=kiargs[++k])>=0) { // l0
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getL0());
    } else y_error("l0 is read only");
  }

  if ((iarg=kiargs[++k])>=0) { // Wsurface
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getWsurface());
    } else y_error("Wsurface is read only");
  }

  if ((iarg=kiargs[++k])>=0) { // Wcentre
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getWcentre());
    } else y_error("Wcentre is read only");
  }

  if ((iarg=kiargs[++k])>=0) { // rcusp
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getRcusp());
    } else y_error("rcusp is read only");
  }

  if ((iarg=kiargs[++k])>=0) { // rcentre
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*ao)->getRcentre());
    } else y_error("rcentre is read only");
  }

  // Call generic Astrobj worker
  ygyoto_Astrobj_generic_eval(ao_, kiargs+k+1, piargs, rvset, paUsed, unit);
}


extern "C" {
  void Y__gyoto_PolishDoughnut_register_as_Astrobj(){
    ygyoto_Astrobj_register("PolishDoughnut",&ygyoto_PolishDoughnut_eval);
  }

  // Generic constructor/accessor
  void
  Y_gyoto_PolishDoughnut(int argc)
  {
    SmartPointer<Astrobj::Generic> *ao = NULL;
    //    char *obj_type=(char*)yget_obj(argc-1,0);
    //    if (obj_type && //!strcmp(obj_type, "gyoto_Astrobj")) {
    //    if (yget_obj(argc-1,0) && yarg_typeid(argc-1)==Y_OPAQUE) {
    if (yarg_Astrobj(argc-1)) {
      ao = yget_Astrobj(--argc);
      if ((*ao)->getKind()!="PolishDoughnut")
	y_error("Expecting Astrobj of kind PolishDoughnut");
      
    }
    ygyoto_PolishDoughnut_eval(ao, argc);
  }

}
