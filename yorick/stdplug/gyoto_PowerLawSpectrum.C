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

#include "ygyoto.h"
#include "yapi.h"
#include <iostream>
#include "GyotoFactory.h"
#include "GyotoPowerLawSpectrum.h"

namespace YGyoto {
  namespace Spectrum {
    ygyoto_Spectrum_eval_worker_t PowerLawYEval;
  }
}

using namespace Gyoto;
using namespace YGyoto;
using namespace Gyoto::Spectrum;
using namespace YGyoto::Spectrum;
using namespace std;

void YGyoto::Spectrum::PowerLawYEval(SmartPointer<Generic> * sp_, int argc) {
  if (debug()) cerr << "in PowerLawYEval()" << endl;
  int rvset[1]={0}, paUsed[1]={0}, constructor=0;

  // If needed, create the object.
  if (!sp_) { // Constructor mode
    constructor=1;
    sp_ = ypush_Spectrum();
  } else *ypush_Spectrum()=*sp_;

  // Parse arguments
  static char const * knames[]={
    "constant", "exponent",
    YGYOTO_SPECTRUM_GENERIC_KW,
    0
  };
#define nkw 2
  static long kglobs[YGYOTO_SPECTRUM_GENERIC_KW_N+nkw+1];
  int kiargs[YGYOTO_SPECTRUM_GENERIC_KW_N+nkw];
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

  // Constructor mode from XML file
  if (constructor) {
    if (yarg_string(piargs[0])) {
#ifdef GYOTO_USE_XERCES
      *sp_ = Factory(ygets_q(piargs[0])).getSpectrum();
      *paUsed=1;
#else
      y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
    } else *sp_ = new PowerLaw();
  }
  SmartPointer<PowerLaw> *sp = (SmartPointer<PowerLaw> *)sp_;

  // Process specific keywords
  int k=-1;
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";

  //// MEMBERS ////
  /* CONSTANT */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*sp)->getConstant());
    } else
      (*sp)->setConstant(ygets_d(iarg)) ;
  }
  
  /* EXPONENT */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) {
      if ((*rvset)++) y_error(rmsg);
      ypush_double((*sp)->getExponent());
    } else
      (*sp)->setExponent(ygets_d(iarg)) ;
  }
  
  
  // GENERIC WORKER
  ygyoto_Spectrum_generic_eval(sp_, kiargs+k+1, piargs, rvset, paUsed);

}

extern "C" {
  void
  Y_gyoto_PowerLawSpectrum(int argc)
  {
    if (debug()) cerr << "In Y_gyoto_PowerLawSpectrum" << endl;
    SmartPointer<Generic> *sp = NULL;
    if (yarg_Spectrum(argc-1)) {
      sp = yget_Spectrum(--argc);
      if ((*sp)->getKind().compare("PowerLaw"))
	y_error("Expecting Spectrum of kind PowerLaw");
    }
    PowerLawYEval(sp, argc);
  }

  void Y__gyoto_PowerLawSpectrum_register_as_Metric(){
    ygyoto_Spectrum_register("PowerLaw",&PowerLawYEval);
  }
}
