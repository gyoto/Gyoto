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
using namespace Gyoto::Spectrometer;

namespace YGyoto {
  void SpectroUniformYEval(Gyoto::SmartPointer<Spectrometer::Generic>*sp, int argc);
}
using namespace YGyoto;

static char const * ygyoto_Spectrometer_names[YGYOTO_MAX_REGISTERED]
={0};
static ygyoto_Spectrometer_eval_worker_t *ygyoto_Spectrometer_evals[YGYOTO_MAX_REGISTERED]
={0};
static int ygyoto_Spectrometer_count=0;




void YGyoto::SpectroUniformYEval(Gyoto::SmartPointer<Spectrometer::Generic>*sp_, int argc) {
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
    "xmlwrite",
    "channels", "midpoints", "widths",
    0
  };
#define nkw 8
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
    if ((*rvset)++) y_error(rmsg);
    size_t nsamples = (*sp) -> getNSamples();
    if (nsamples) {
      long dims[] = {1, nsamples+1};
      double const * const nuchannels = (*sp) -> getChannelBoundaries();
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
    SmartPointer<Spectrometer::Generic> spectro;
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
    SmartPointer<Spectrometer::Generic> sp = (((gyoto_Spectro*)obj)->spectro);
    // If no parameters, return pointer
    if (argc==1 && yarg_nil(0)) {
      ypush_long( (long) sp() );
      return;
    }

    // Try calling kind-specific worker
    int n=0;
    char const * const  kind = sp->getKind();

    while (n<ygyoto_Spectrometer_count &&
	   kind != ygyoto_Spectrometer_names[n]) ++n;

    if (n<ygyoto_Spectrometer_count && ygyoto_Spectrometer_evals[n]) {
      (*ygyoto_Spectrometer_evals[n])(&sp, argc);
      return;
    }

    // Fall-back to generic worker
    static char const * knames[]={
      "unit",
      YGYOTO_SPECTROMETER_GENERIC_KW, 0
    };
    static long kglobs[YGYOTO_SPECTROMETER_GENERIC_KW_N+2];
    int kiargs[YGYOTO_SPECTROMETER_GENERIC_KW_N+1];
    int piargs[]={-1,-1,-1,-1};
    // push back spectrometer by default
    *ypush_Spectrometer()=sp;
    yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
    
    int iarg=argc, parg=0;
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<4) piargs[parg++]=iarg--;
	else y_error("gyoto_Spectrometer takes at most 4 positional arguments");
      }
    }

    int rvset[1]={0}, paUsed[1]={0};
    char * unit=NULL;
    int k=-1;

    /* UNIT */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      GYOTO_DEBUG << "get unit" << endl;
      unit = ygets_q(iarg);
    }
  
    ygyoto_Spectrometer_generic_eval(&sp, kiargs+k+1, piargs, rvset, paUsed, unit);

  }
  static y_userobj_t gyoto_Spectro_obj =
    {const_cast<char*>("gyoto_Spectrometer"), &gyoto_Spectro_free, &gyoto_Spectro_print,
     &gyoto_Spectro_eval, 0, 0};

  void
  Y_gyoto_Spectrometer(int argc)
  {
    if (debug()) cerr << "In Y_gyoto_Spectrometer" << endl;
    int rvset[1]={0}, paUsed[1]={0};
    SmartPointer<Spectrometer::Generic> *sp = NULL;
    int builder=0;

    if (yarg_Spectrometer(argc-1)) {
      sp = yget_Spectrometer(--argc);
      // Try calling kind-specific worker
      int n=0;
      char const * const kind = (*sp)->getKind();
      while (n<ygyoto_Spectrometer_count &&
	     kind != ygyoto_Spectrometer_names[n])
	++n;
      if (n<ygyoto_Spectrometer_count &&
	  ygyoto_Spectrometer_evals[n]) {
	(*ygyoto_Spectrometer_evals[n])(sp, argc);
	return;
      }
    
      // push back Spectrometer
      *ypush_Spectrometer()=*sp;
    } else { // Constructor mode
      sp = ypush_Spectrometer();
      *rvset=1;
    }

    static char const * knames[]={
      "unit",
      YGYOTO_SPECTROMETER_GENERIC_KW, 0
    };
    static long kglobs[YGYOTO_SPECTROMETER_GENERIC_KW_N+2];
    int kiargs[YGYOTO_SPECTROMETER_GENERIC_KW_N+1];
    int piargs[]={-1,-1,-1,-1};
    yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
    
    int iarg=argc, parg=0;
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<4) piargs[parg++]=iarg--;
	else y_error("gyoto_Spectrometer takes at most 4 positional arguments");
      }
    }

    // if builder==1, constructor mode:
    if (builder) {
      if (yarg_string(piargs[0])) {
#ifdef GYOTO_USE_XERCES
	char * fname = ygets_q(piargs[0]);
	Spectrometer::Subcontractor_t *sub = Spectrometer::getSubcontractor(fname, 1);
	paUsed[0]=1;
	if (sub) *sp=(*sub)(NULL);
	else y_error("not implemented");
#else
	y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
      } else y_error("Cannot allocate object of virtual class Spectrometer");
    }

    char * unit=NULL;
    int k=-1;

    /* UNIT */
    if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      GYOTO_DEBUG << "get unit" << endl;
      unit = ygets_q(iarg);
    }

    ygyoto_Spectrometer_generic_eval(sp, kiargs+k+1, piargs, rvset, paUsed, unit);
    //  SpectroYEval(sp, argc);
  }

  void Y__gyoto_SpectroUniform_register_as_Spectro(int argc){
    ygyoto_Spectrometer_register(GYOTO_SPECTRO_KIND_WAVE,
				 &YGyoto::SpectroUniformYEval);
    ygyoto_Spectrometer_register(GYOTO_SPECTRO_KIND_WAVELOG,
				 &YGyoto::SpectroUniformYEval);
    ygyoto_Spectrometer_register(GYOTO_SPECTRO_KIND_FREQ,
				 &YGyoto::SpectroUniformYEval);
    ygyoto_Spectrometer_register(GYOTO_SPECTRO_KIND_FREQLOG,
				 &YGyoto::SpectroUniformYEval);
  }

  // Uniform CLASS
  // Constructor

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
	SpectroKind_t kind=(*sp)->getKind();
	if (kind != GYOTO_SPECTRO_KIND_WAVE &&
	    kind != GYOTO_SPECTRO_KIND_WAVELOG &&
	    kind != GYOTO_SPECTRO_KIND_FREQ &&
	    kind != GYOTO_SPECTRO_KIND_FREQLOG)
	  y_error("Expecting Spectrometer of kind Uniform");
      }
    } YGYOTO_STD_CATCH;
    YGyoto::SpectroUniformYEval(sp, argc);
  }

}

int
yarg_Spectrometer(int iarg) {
  return yget_obj(iarg,0)==gyoto_Spectro_obj.type_name;
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

void ygyoto_Spectrometer_generic_eval(SmartPointer<Spectrometer::Generic>*sp,
				int *kiargs, int *piargs,
				int *rvset, int *paUsed, char * unit) {
  int k=-1, iarg=-1;
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";
  long  ntot, dims[Y_DIMSIZE];

  /* METHODS */  

  // kind
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    if (!yarg_nil(iarg)) y_error("KIND is readonly");
    char ** kind = ypush_q(0);
    *kind = p_strcpy((*sp)->getKind());
  }

  // Process SET keywords
  // Save to file
  if ((iarg=kiargs[++k])>=0) { // xmlwrite
    iarg+=*rvset;
    char *filename=ygets_q(iarg);
#ifdef GYOTO_USE_XERCES
    Factory(*sp).write(filename);
#else
    y_error("This GYOTO was compiled without xerces: no xml i/o");
#endif
  }

  /* CLONE */
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    *ypush_Spectrometer() = (*sp)->clone();
  }

  // get members //
  /* NSAMPLES */
  if ((iarg=kiargs[++k])>=0) {
      iarg+=*rvset;
      if (yarg_nil(iarg)) ypush_long((*sp)->getNSamples());
      else y_error("NSAMPLES is readonly");
  }

    /* GET SPECTRO CHANNELS */
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    size_t nsamples = (*sp) -> getNSamples();
    if (nsamples) {
      long dims[] = {2, nsamples, 2};
      double const * const nuchannels = (*sp) -> getChannelBoundaries();
      size_t const * const chanind = (*sp) -> getChannelIndices();
      double * ychannels = ypush_d(dims);
      for (size_t i=0; i<=2*nsamples; ++i) {
	ychannels[i] = nuchannels[chanind[i]];
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
      double * ychannels = ypush_d(dims);
      for (size_t i=0; i<nsamples; ++i) {
	ychannels[i] = nuchannels[i];
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
      double * ywidths = ypush_d(dims);
      for (size_t i=0; i<nsamples; ++i) {
	ywidths[i] = nuwidths[i];
      }
    } else ypush_nil();
  } 

}

SmartPointer<Spectrometer::Generic>* yget_Spectrometer(int iarg) {
  return &(((gyoto_Spectro*)yget_obj(iarg, &gyoto_Spectro_obj))->spectro);
}

SmartPointer<Spectrometer::Generic>* ypush_Spectrometer() {
  return &(((gyoto_Spectro*)(ypush_obj(&gyoto_Spectro_obj,
					sizeof(gyoto_Spectro))))->spectro);
}

