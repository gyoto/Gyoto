/*
    Copyright 2011-2013 Thibaut Paumard

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

#include <cmath>
#include "ygyoto.h"
#include "ygyoto_private.h"
#include "yapi.h"
#include "pstdlib.h"
#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#endif
using namespace std;
using namespace Gyoto;
using namespace Gyoto::Spectrometer;

static char const * ygyoto_Spectrometer_names[YGYOTO_MAX_REGISTERED]
={0};
static ygyoto_Spectrometer_eval_worker_t *ygyoto_Spectrometer_evals[YGYOTO_MAX_REGISTERED]
={0};
static int ygyoto_Spectrometer_count=0;

YGYOTO_YUSEROBJ(Spectrometer, Spectrometer::Generic)
YGYOTO_BASE_CONSTRUCTOR1(Spectrometer,spectrometer)

extern "C" {
  void gyoto_Spectrometer_eval(void *obj, int argc) {
    SmartPointer<Spectrometer::Generic> *OBJ_ =
      &(((gyoto_Spectrometer*)obj)->smptr);
    // If no parameters, return pointer
    if (argc==1 && yarg_nil(0)) {
      ypush_long( (long) (*OBJ_)() );
      return;
    }

    // Try calling kind-specific worker
    int n=0;
    char const * const  kind = (*OBJ_)->kindid();

    while (n<ygyoto_Spectrometer_count &&
	   kind != ygyoto_Spectrometer_names[n]) ++n;

    if (n<ygyoto_Spectrometer_count && ygyoto_Spectrometer_evals[n]) {
      (*ygyoto_Spectrometer_evals[n])(OBJ_, argc);
      return;
    }

    // Fall-back to generic worker
    static char const * knames[]={
      "unit",
      YGYOTO_SPECTROMETER_GENERIC_KW, 0
    };

    YGYOTO_WORKER_INIT(Spectrometer,
		       Generic, knames, YGYOTO_METRIC_GENERIC_KW_N+1);

    YGYOTO_WORKER_SET_UNIT;

    YGYOTO_WORKER_CALL_GENERIC(Spectrometer);
  }

}


void ygyoto_Spectrometer_register(char const*const name, ygyoto_Spectrometer_eval_worker_t* on_eval){
  int n;
  if (ygyoto_Spectrometer_count==YGYOTO_MAX_REGISTERED)
    y_error("Too many Spectrometers registered");
  for (n=0; n<ygyoto_Spectrometer_count; ++n)
    if (ygyoto_Spectrometer_names[n]==name) return;

  ygyoto_Spectrometer_names[ygyoto_Spectrometer_count] = name;
  ygyoto_Spectrometer_evals[ygyoto_Spectrometer_count++]=on_eval;
  //if ((ygyoto_Spectrometer_count) < YGYOTO_METRIC_MAX_REGISTERED)
  //  strcpy(ygyoto_Spectrometer_names[ygyoto_Spectrometer_count], "");
}

void ygyoto_Spectrometer_generic_eval(SmartPointer<Spectrometer::Generic>*OBJ,
				int *kiargs, int *piargs,
				int *rvset, int *paUsed, char * unit) {
  int k=-1, iarg=-1;
  char const * rmsg="Cannot set return value more than once";

  /* METHODS */  

  // kind
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    if (!yarg_nil(iarg)) y_error("KIND is readonly");
    char ** kind = ypush_q(0);
    *kind = p_strcpy((*OBJ)->kindid());
  }

  // Process SET keywords
  // Save to file
  YGYOTO_WORKER_XMLWRITE;

  /* CLONE */
  YGYOTO_WORKER_CLONE(Spectrometer);

  /* HELP */
  YGYOTO_WORKER_HELP;

  // get members //
  /* NSAMPLES */
  if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) ypush_long((*OBJ)->nSamples());
      else y_error("NSAMPLES is readonly");
  }

  /* SETPARAMETER */
  YGYOTO_WORKER_SETPARAMETER;

  /* GET SPECTRO CHANNELS */
  if ((iarg=kiargs[++k])>=0) {
    if (!yarg_nil(iarg)) y_error("CHANNELS is readonly");
    if ((*rvset)++) y_error(rmsg);
    long nsamples = long((*OBJ) -> nSamples());
    if (nsamples) {
      long dims[] = {2, 2, nsamples};
      double converted[(*OBJ)->getNBoundaries()];
      (*OBJ) -> getChannelBoundaries(converted, unit?unit:"");
      GYOTO_DEBUG_ARRAY(converted, (*OBJ)->getNBoundaries());
      size_t const * const chanind = (*OBJ) -> getChannelIndices();
      GYOTO_DEBUG_ARRAY(chanind, 2*size_t(nsamples));
      double * ychannels = ypush_d(dims);
      for (long i=0; i<2*nsamples; ++i) {
	ychannels[i] = converted[chanind[i]];
	GYOTO_DEBUG << "ychannels["<< i << "]=" << ychannels[i] << endl;
      }

    } else ypush_nil();
  } 

  /* GET SPECTRO MIDPOINTS */
  if ((iarg=kiargs[++k])>=0) {
    if (!yarg_nil(iarg)) y_error("MIDPOINTS is readonly");
    if ((*rvset)++) y_error(rmsg);
    long nsamples = long((*OBJ) -> nSamples());
    if (nsamples) {
      long dims[] = {1, nsamples};
      double * ychannels = ypush_d(dims);
      (*OBJ)->getMidpoints(ychannels, unit?unit:"");
    } else ypush_nil();
  }

  /* GET SPECTRO WIDTHS */
  if ((iarg=kiargs[++k])>=0) {
    if (!yarg_nil(iarg)) y_error("WIDTHS is readonly");
    if ((*rvset)++) y_error(rmsg);
    long nsamples = long((*OBJ) -> nSamples());
    if (nsamples) {
      long dims[] = {1, nsamples};
      double * ywidths = ypush_d(dims);
      (*OBJ)->getWidths(ywidths, unit?unit:"");
    } else ypush_nil();
  } 

}

