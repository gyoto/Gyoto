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

#include <cmath>
#include "ygyoto.h"
#include "yapi.h"
#include "pstdlib.h"
#ifdef GYOTO_USE_XERCES
#include "GyotoFactory.h"
#endif
using namespace std;
using namespace Gyoto;

namespace YGyoto {
  void SpectroYEval(Gyoto::SmartPointer<Spectrometer>*sp, int argc);
}
using namespace YGyoto;

void YGyoto::SpectroYEval(Gyoto::SmartPointer<Spectrometer>*sp, int argc) {
  int k=-1, rvset[1]={0}, paUsed[1]={0};
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";

  if (!sp) {
    sp = ypush_Spectrometer();
    *sp = new Spectrometer();
  } else {
    *ypush_Spectrometer() = *sp;
  }

  static char const * knames[]={
    "kind", "nsamples", "band", 
    "xmlwrite",
    "channels", "midpoints", "widths",
    0
  };
#define nkw 7
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

  /* SPECTRO_KIND */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (yarg_nil(iarg)) { // get spectro kind
      if ((*rvset)++) y_error(rmsg);
      *ypush_q(0) = p_strcpy((*sp) -> getKindStr() . c_str());
    } else { // set spectro kind
      (*sp) -> setKind(ygets_q(iarg));
    }
  }

  /* NSAMPLES */
  if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) ypush_long((*sp)->getNSamples());
      else (*sp) -> setNSamples( ygets_l(iarg) );
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
      (*sp) -> setBand(boundaries);
    }
  }

  // Save to file
  if ((iarg=kiargs[++k])>=0) { // xmlwrite
    iarg+=*rvset;
#ifdef GYOTO_USE_XERCES
    char *filename=ygets_q(iarg);
    Factory(*sp).write(filename);
#else
    y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
  }

  /* GET SPECTRO CHANNELS */
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    size_t nsamples = (*sp) -> getNSamples();
    if (nsamples) {
      long dims[] = {1, nsamples+1};
      double const * const nuchannels = (*sp) -> getChannels();
      SpectroKind_t kind= (*sp) -> getKind();
      double * ychannels = ypush_d(dims);
      for (size_t i=0; i<=nsamples; ++i) {
	ychannels[i] = nuchannels[i];
	if (kind==GYOTO_SPECTRO_KIND_WAVE ||
	    kind==GYOTO_SPECTRO_KIND_WAVELOG)
	  ychannels[i]=GYOTO_C/ychannels[i];
	if (kind==GYOTO_SPECTRO_KIND_FREQLOG ||
	    kind==GYOTO_SPECTRO_KIND_WAVELOG)
	  ychannels[i]=log10(ychannels[i]);
      }

    } else ypush_nil();
  } 

  /* GET SPECTRO MIDPOINTS */
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    size_t nsamples = (*sp) -> getNSamples();
    if (nsamples) {
      long dims[] = {1, nsamples};
      double const * const nuchannels = (*sp) -> getMidpoints();
      SpectroKind_t kind= (*sp) -> getKind();
      double * ychannels = ypush_d(dims);
      for (size_t i=0; i<nsamples; ++i) {
	ychannels[i] = nuchannels[i];
	if (kind==GYOTO_SPECTRO_KIND_WAVE ||
	    kind==GYOTO_SPECTRO_KIND_WAVELOG)
	  ychannels[i]=GYOTO_C/ychannels[i];
	if (kind==GYOTO_SPECTRO_KIND_FREQLOG ||
	    kind==GYOTO_SPECTRO_KIND_WAVELOG)
	  ychannels[i]=log10(ychannels[i]);
      }

    } else ypush_nil();
  }

  /* GET SPECTRO WIDTHS */
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    size_t nsamples = (*sp) -> getNSamples();
    if (nsamples) {
      long dims[] = {1, nsamples};
      double const * const nuwidths = (*sp) -> getWidths();
      SpectroKind_t kind= (*sp) -> getKind();
      double * ywidths = ypush_d(dims);
      for (size_t i=0; i<nsamples; ++i) {
	ywidths[i] = nuwidths[i];
      }
    } else ypush_nil();
  } 

}


extern "C" {
  typedef struct gyoto_Spectro {
    SmartPointer<Spectrometer> spectro;
  } gyoto_Spectro;

  
  void gyoto_Spectro_free(void *obj) {
    if (((gyoto_Spectro*)obj)->spectro) {
      ((gyoto_Spectro*)obj)->spectro = NULL;
    } else printf("null pointer\n");
  }

  void gyoto_Spectro_print(void *obj) {
#ifdef GYOTO_USE_XERCES
    string rest="", sub="";
    size_t pos=0, len=0;
    rest = Factory(((gyoto_Spectro*)obj)->spectro).format();
    while (len=rest.length())  {
      sub=rest.substr(0, pos=rest.find_first_of("\n",0));
      rest=rest.substr(pos+1, len-1);
      y_print( sub.c_str(),1 );
    }
#else
    y_print("GYOTO Spectrum object of type ",0);
    y_print(((gyoto_Spectro*)obj)->spectro->getKindStr().c_str(),0);
#endif
  }
  void gyoto_Spectro_eval(void *obj, int argc) {
    // If no parameters, return pointer
    if (argc==1 && yarg_nil(0)) {
      ypush_long((long)((gyoto_Spectro*)obj)->spectro());
      return;
    }
    SmartPointer<Spectrometer> * sp = &(((gyoto_Spectro*)obj)->spectro);
    YGyoto::SpectroYEval(sp, argc);
  }
  static y_userobj_t gyoto_Spectro_obj =
    {const_cast<char*>("gyoto_Spectrometer"), &gyoto_Spectro_free, &gyoto_Spectro_print,
     &gyoto_Spectro_eval, 0, 0};

  void
  Y_gyoto_Spectrometer(int argc)
  {
    if (debug()) cerr << "In Y_gyoto_Spectrometer" << endl;
    SmartPointer<Spectrometer> *sp = NULL;
    if (yarg_Spectrometer(argc-1)) {
      sp = yget_Spectrometer(--argc);
    }
    SpectroYEval(sp, argc);
  }
}

int
yarg_Spectrometer(int iarg) {
  return yget_obj(iarg,0)==gyoto_Spectro_obj.type_name;
}

SmartPointer<Spectrometer>* yget_Spectrometer(int iarg) {
  return &(((gyoto_Spectro*)yget_obj(iarg, &gyoto_Spectro_obj))->spectro);
}

SmartPointer<Spectrometer>* ypush_Spectrometer() {
  return &(((gyoto_Spectro*)(ypush_obj(&gyoto_Spectro_obj,
					sizeof(gyoto_Spectro))))->spectro);
}

