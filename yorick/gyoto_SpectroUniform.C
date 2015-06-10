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

#define OBJ sp

namespace YGyoto {
  void SpectroUniformYEval(SmartPointer<Spectrometer::Generic>*sp, int argc);
}
using namespace YGyoto;

void
YGyoto::SpectroUniformYEval(Gyoto::SmartPointer<Spectrometer::Generic>*sp_,
			    int argc)
{

  static char const * knames[]={
    "unit",
    "kind", "nsamples", "band", 
    "setparameter",
    "xmlwrite", "clone",
    "channels", "midpoints", "widths",
    0
  };
  YGYOTO_WORKER_INIT(Spectrometer, Uniform,
		     knames, 10);

  YGYOTO_WORKER_SET_UNIT;

  /* SPECTRO_KIND */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) { // get spectro kind
      if ((*rvset)++) y_error(rmsg);
      *ypush_q(0) = p_strcpy((*sp) -> kindid() );
    } else { // set spectro kind
      (*sp) -> kind(std::string(ygets_q(iarg)));
    }
  }

  YGYOTO_WORKER_GETSET_LONG2(nSamples);

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
      (*sp) -> band(boundaries, unit?unit:"");
    }
  }

  YGYOTO_WORKER_SETPARAMETER;
#undef OBJ
#define OBJ sp_
  YGYOTO_WORKER_XMLWRITE;
#undef OBJ
#define OBJ sp
  YGYOTO_WORKER_CLONE(Spectrometer);

  /* GET SPECTRO CHANNELS */
  if ((iarg=kiargs[++k])>=0) {
    if (!yarg_nil(iarg)) y_error("CHANNELS is readonly");
    if ((*rvset)++) y_error(rmsg);
    long nsamples = long((*sp) -> nSamples());
    if (nsamples) {
      long dims[] = {2, 2, nsamples};
      double converted[(*sp)->getNBoundaries()];
      (*sp) -> getChannelBoundaries(converted, unit?unit:"");
      GYOTO_DEBUG_ARRAY(converted, (*sp)->getNBoundaries());
      size_t const * const chanind = (*sp) -> getChannelIndices();
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
    long nsamples = long((*sp) -> nSamples());
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
    long nsamples = long((*sp) -> nSamples());
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
    YGYOTO_CONSTRUCTOR_INIT2(Spectrometer, Spectrometer::Generic,
			     Uniform, spectrometer);
    kind_t kind=(*sp)->kindid();
    if (kind != Uniform::WaveKind &&
	kind != Uniform::WaveLogKind &&
	kind != Uniform::FreqKind &&
	kind != Uniform::FreqLogKind)
      y_error("Expecting Spectrometer of kind Uniform");
    YGyoto::SpectroUniformYEval(OBJ, argc);
  }

}
