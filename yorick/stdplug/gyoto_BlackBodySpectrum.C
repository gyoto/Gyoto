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
#include "GyotoBlackBodySpectrum.h"

namespace YGyoto {
  namespace Spectrum {
    ygyoto_Spectrum_eval_worker_t BlackBodyYEval;
  }
}

using namespace Gyoto;
using namespace YGyoto;
using namespace Gyoto::Spectrum;
using namespace YGyoto::Spectrum;
using namespace std;

#define OBJ sp

void YGyoto::Spectrum::BlackBodyYEval(SmartPointer<Generic> * sp_, int argc) {
  if (debug()) cerr << "in BlackBodyYEval()" << endl;
  int rvset[1]={0}, paUsed[1]={0}, constructor=0;

  // If needed, create the object.
  if (!sp_) { // Constructor mode
    constructor=1;
    sp_ = ypush_Spectrum();
  } else *ypush_Spectrum()=*sp_;

  // Parse arguments
  static char const * knames[]={
    "temperature", "scaling",
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
    } else *sp_ = new BlackBody();
  }
  SmartPointer<BlackBody> *sp = (SmartPointer<BlackBody> *)sp_;

  // Process specific keywords
  int k=-1;
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";

  YGYOTO_WORKER_GETSET_DOUBLE(Temperature);
  YGYOTO_WORKER_GETSET_DOUBLE(Scaling);
  
  // GENERIC WORKER
  ygyoto_Spectrum_generic_eval(sp_, kiargs+k+1, piargs, rvset, paUsed);

}

extern "C" {
  void
  Y_gyoto_BlackBodySpectrum(int argc)
  {
    if (debug()) cerr << "In Y_gyoto_BlackBodySpectrum" << endl;
    SmartPointer<Generic> *sp = NULL;
    if (yarg_Spectrum(argc-1)) {
      sp = yget_Spectrum(--argc);
      if ((*sp)->getKind().compare("BlackBody"))
	y_error("Expecting Spectrum of kind BlackBody");
    }
    BlackBodyYEval(sp, argc);
  }

  void Y__gyoto_BlackBodySpectrum_register_as_Metric(){
    ygyoto_Spectrum_register("BlackBody",&BlackBodyYEval);
  }
}
