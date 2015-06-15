/*
    Copyright 2011, 2013 Thibaut Paumard

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
#include "ygyoto_private.h"
#include "yapi.h"
#include "pstdlib.h"

#ifdef GYOTO_USE_XERCES
#include <GyotoFactory.h>
#endif

#include <cstring>
using namespace std;
using namespace Gyoto;
using namespace Gyoto::Spectrum;

static char ygyoto_Spectrum_names[YGYOTO_TYPE_LEN][YGYOTO_MAX_REGISTERED]
={{0}};
static ygyoto_Spectrum_eval_worker_t *ygyoto_Spectrum_evals[YGYOTO_MAX_REGISTERED]
={0};
static int ygyoto_Spectrum_count=0;

YGYOTO_YUSEROBJ(Spectrum, Spectrum::Generic)
YGYOTO_BASE_CONSTRUCTOR1(Spectrum,spectrum)

extern "C" {
  void gyoto_Spectrum_eval(void *obj, int argc) {
    // If no parameters, return pointer
    if (argc==1 && yarg_nil(0)) {
      ypush_long((long)((gyoto_Spectrum*)obj)->smptr());
      return;
    }

    // Try calling kind-specific worker
    int n=0;
    SmartPointer<Spectrum::Generic> * OBJ_ = &(((gyoto_Spectrum*)obj)->smptr);
    const string kind = (*OBJ_)->kind();

    while (n<ygyoto_Spectrum_count && kind.compare(ygyoto_Spectrum_names[n]))
      ++n;

    if (n<ygyoto_Spectrum_count && ygyoto_Spectrum_evals[n]) {
      (*ygyoto_Spectrum_evals[n])(OBJ_, argc);
      return;
    }

    // Fall-back to default worker
    static char const * knames[]={
      "unit",
      YGYOTO_SPECTRUM_GENERIC_KW, 0
    };
    YGYOTO_WORKER_INIT(Spectrum, Generic,
		       knames, YGYOTO_SPECTRUM_GENERIC_KW_N+1);
    YGYOTO_WORKER_SET_UNIT;
    YGYOTO_WORKER_CALL_GENERIC(Spectrum);

  }
}

void ygyoto_Spectrum_register(char const*const name, ygyoto_Spectrum_eval_worker_t* on_eval){
  int n;
  if (ygyoto_Spectrum_count==YGYOTO_MAX_REGISTERED)
    y_error("Too many Spectra registered");
  for (n=0; n<ygyoto_Spectrum_count; ++n)
    if (!strcmp(ygyoto_Spectrum_names[n], name)) 
      return;

  strcpy(ygyoto_Spectrum_names[ygyoto_Spectrum_count], name);
  ygyoto_Spectrum_evals[ygyoto_Spectrum_count++]=on_eval;
}

void ygyoto_Spectrum_generic_eval(Gyoto::SmartPointer<Generic>*OBJ,
				int *kiargs, int *piargs,
				  int *rvset, int *paUsed, char* unit) {
  int k=-1, iarg;
  char const * rmsg="Cannot set return value more than once";

  if (debug())
    for (int i=0; i<YGYOTO_SPECTRUM_GENERIC_KW_N; ++i)
      cerr << "DEBUG: Spectrum_generic_eval: kiargs[" << i << "]="
	   << kiargs[i] << endl;

  YGYOTO_WORKER_XMLWRITE;

  // kind
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (!yarg_nil(iarg)) y_error("KIND is readonly");
    if (debug()) cerr << "kiargs=" << kiargs << endl;
    if ((*rvset)++) y_error(rmsg);
    char ** kind = ypush_q(0);
    *kind = p_strcpy((*OBJ)->kind().c_str());
  }

  YGYOTO_WORKER_SETPARAMETER;
  YGYOTO_WORKER_CLONE(Spectrum);
  YGYOTO_WORKER_HELP;

  /* INTEGRATE */
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    long ntot, dims[Y_DIMSIZE];
    double * freqs = ygeta_d(iarg, &ntot, dims);
    if (dims[0] !=1)
      y_error("gyoto_Spectrum(integrate=FREQS): FREQS must be 1-dimensional");
    if (ntot < 2)
      y_error("gyoto_Spectrum(integrate=FREQS): FREQS must gave >=2 elements");
    --dims[1]; --ntot;
    double * Inu1nu2 = ypush_d(dims);
    for (long i=0; i < ntot; ++i)
      Inu1nu2[i] = (*OBJ)->integrate(freqs[i], freqs[i+1]);
  }

  // GET SPECTRUM VALUE FOR WAVELENGTHS
  if (*rvset || *paUsed || (iarg=piargs[0]) < 0 || !yarg_number(iarg)) return;

  if (debug())
    cerr<<"DEBUG: gyoto_Spectrum_generic_eval evaluating Spectrum "
      "at frequency\n";
  long ntot, dims[Y_DIMSIZE];
  double * freqs = ygeta_d(iarg, &ntot, dims);
  double * Inu = ypush_d(dims);
  for (long i=0; i < ntot; ++i) Inu[i] = (**OBJ)(freqs[i]);

  if (debug()) cerr << "DEBUG: out of Spectrum_generic_eval"<< endl;

}

