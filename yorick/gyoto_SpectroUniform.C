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
#include "GyotoUniformSpectrometer.h"
#include "yapi.h"
#include "pstdlib.h"
#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#endif
using namespace std;
using namespace Gyoto;
using namespace Gyoto::Spectrometer;

namespace YGyoto {
  void SpectroUniformYEval(Gyoto::SmartPointer<Spectrometer::Generic>*sp, int argc);
}
using namespace YGyoto;

void
YGyoto::SpectroUniformYEval(Gyoto::SmartPointer<Spectrometer::Generic>*sp_,
			    int argc)
{
  int k=-1, rvset[1]={0}, paUsed[1]={0};
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";

  if (!sp_) {
    sp_ = ypush_Spectrometer();
    *sp_ = new Spectrometer::Uniform();
  } else {
    *ypush_Spectrometer() = *sp_;
  }

  SmartPointer<Uniform> *sp = (SmartPointer<Uniform> *)sp_;
  static char const * knames[]={
    "unit",
    "kind", "nsamples", "band", 
    "setparameter",
    "xmlwrite",
    "channels", "midpoints", "widths",
    0
  };
#define nkw 9
  static long kglobs[nkw+1];
  int kiargs[nkw];
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

  char *unit = NULL;
  /* UNIT */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    GYOTO_DEBUG << "get unit" << endl;
    unit = ygets_q(iarg);
  }

  /* SPECTRO_KIND */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) { // get spectro kind
      if ((*rvset)++) y_error(rmsg);
      *ypush_q(0) = p_strcpy((*sp) -> getKind() );
    } else { // set spectro kind
      (*sp) -> setKind(std::string(ygets_q(iarg)));
    }
  }

  /* NSAMPLES */
  if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) ypush_long((*sp)->getNSamples());
      else {
	(*sp) -> setNSamples( ygets_l(iarg) );
      }
  }

  /* BAND */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    double * boundaries = NULL;
    if (yarg_nil(iarg)) { // get spectro band
      if ((*rvset)++) y_error(rmsg);
      long dims[Y_DIMSIZE] = {1, 2};
      boundaries=ypush_d(dims);
      boundaries[0] = (*sp) -> getBand()[0];
      boundaries[1] = (*sp) -> getBand()[1];
    } else { // set spectro band
      long ntot;
      boundaries = ygeta_d(iarg, &ntot, NULL);
      if (ntot != 2)
	  y_error("BAND must have 2 elements");
      (*sp) -> setBand(boundaries, unit?unit:"");
    }
  }

  /* SETPARAMETER */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if ((*paUsed)++) y_error("pmsg");
    string name = ygets_q(iarg);
    string content = ygets_q(*piargs);
    (*sp)->setParameter(name, content,  unit?unit:"");
  }

  // Save to file
  if ((iarg=kiargs[++k])>=0) { // xmlwrite
    iarg+=*rvset;
#ifdef GYOTO_USE_XERCES
    char *filename=ygets_q(iarg);
    Factory(*sp_).write(filename);
#else
    y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
  }

    /* GET SPECTRO CHANNELS */
  if ((iarg=kiargs[++k])>=0) {
    if (!yarg_nil(iarg)) y_error("CHANNELS is readonly");
    if ((*rvset)++) y_error(rmsg);
    size_t nsamples = (*sp) -> getNSamples();
    if (nsamples) {
      long dims[] = {2, 2, nsamples};
      double converted[(*sp)->getNBoundaries()];
      (*sp) -> getChannelBoundaries(converted, unit?unit:"");
      GYOTO_DEBUG_ARRAY(converted, (*sp)->getNBoundaries());
      size_t const * const chanind = (*sp) -> getChannelIndices();
      GYOTO_DEBUG_ARRAY(chanind, 2*nsamples);
      double * ychannels = ypush_d(dims);
      for (size_t i=0; i<2*nsamples; ++i) {
	ychannels[i] = converted[chanind[i]];
	GYOTO_DEBUG << "ychannels["<< i << "]=" << ychannels[i] << endl;
      }

    } else ypush_nil();
  } 

  /* GET SPECTRO MIDPOINTS */
  if ((iarg=kiargs[++k])>=0) {
    if (!yarg_nil(iarg)) y_error("MIDPOINTS is readonly");
    if ((*rvset)++) y_error(rmsg);
    size_t nsamples = (*sp) -> getNSamples();
    if (nsamples) {
      long dims[] = {1, nsamples};
      double * ychannels = ypush_d(dims);
      (*sp)->getMidpoints(ychannels, unit?unit:"");
    } else ypush_nil();
  }

  /* GET SPECTRO WIDTHS */
  if ((iarg=kiargs[++k])>=0) {
    if (!yarg_nil(iarg)) y_error("WIDTHS is readonly");
    if ((*rvset)++) y_error(rmsg);
    size_t nsamples = (*sp) -> getNSamples();
    if (nsamples) {
      long dims[] = {1, nsamples};
      double * ywidths = ypush_d(dims);
      (*sp)->getWidths(ywidths, unit?unit:"");
    } else ypush_nil();
  } 

}


extern "C" {

  void Y__gyoto_SpectroUniform_register_as_Spectro(int argc){
    ygyoto_Spectrometer_register(Uniform::WaveKind,
				 &YGyoto::SpectroUniformYEval);
    ygyoto_Spectrometer_register(Uniform::WaveLogKind,
				 &YGyoto::SpectroUniformYEval);
    ygyoto_Spectrometer_register(Uniform::FreqKind,
				 &YGyoto::SpectroUniformYEval);
    ygyoto_Spectrometer_register(Uniform::FreqLogKind,
				 &YGyoto::SpectroUniformYEval);
  }

  void
  Y_gyoto_SpectroUniform(int argc)
  {
    SmartPointer<Spectrometer::Generic> *sp = NULL;
    try {
      //    char *obj_type=(char*)yget_obj(argc-1,0);
      //    if (obj_type && //!strcmp(obj_type, "gyoto_Metric")) {
      //    if (yget_obj(argc-1,0) && yarg_typeid(argc-1)==Y_OPAQUE) {
      if (yarg_Spectrometer(argc-1)) {
	sp = yget_Spectrometer(--argc);
	kind_t kind=(*sp)->getKind();
	if (kind != Uniform::WaveKind &&
	    kind != Uniform::WaveLogKind &&
	    kind != Uniform::FreqKind &&
	    kind != Uniform::FreqLogKind)
	  y_error("Expecting Spectrometer of kind Uniform");
      }
    } YGYOTO_STD_CATCH;
    YGyoto::SpectroUniformYEval(sp, argc);
  }

}
