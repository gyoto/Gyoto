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
#include "pstdlib.h"

#ifdef GYOTO_USE_XERCES
#include <GyotoFactory.h>
#endif

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Spectrum;

static char ygyoto_Spectrum_names[YGYOTO_TYPE_LEN][YGYOTO_MAX_REGISTERED]
={{0}};
static ygyoto_Spectrum_eval_worker_t *ygyoto_Spectrum_evals[YGYOTO_MAX_REGISTERED]
={0};
static int ygyoto_Spectrum_count=0;

extern "C" {
  // SPECTRUM CLASS
  // Opaque Yorick object
  typedef struct gyoto_Spectrum {
    SmartPointer<Spectrum::Generic> spectrum;
    //char type[YGYOTO_TYPE_LEN];
  } gyoto_Spectrum;
  void gyoto_Spectrum_free(void *obj) {
    if (((gyoto_Spectrum*)obj)->spectrum) {
      ((gyoto_Spectrum*)obj)->spectrum = NULL;
    } else printf("null pointer\n");
  }
  void gyoto_Spectrum_print(void *obj) {
#ifdef GYOTO_USE_XERCES
    string rest="", sub="";
    size_t pos=0, len=0;
    rest = Factory(((gyoto_Spectrum*)obj)->spectrum).format();
    while (len=rest.length())  {
      sub=rest.substr(0, pos=rest.find_first_of("\n",0));
      rest=rest.substr(pos+1, len-1);
      y_print( sub.c_str(),1 );
    }
#else
    y_print("GYOTO Spectrum object of type ",0);
    y_print(((gyoto_Spectrum*)obj)->spectrum->getKind().c_str(),0);
#endif
  }
  void gyoto_Spectrum_eval(void *obj, int argc) {
    // If no parameters, return pointer
    if (argc==1 && yarg_nil(0)) {
      ypush_long((long)((gyoto_Spectrum*)obj)->spectrum());
      return;
    }

    // Try calling kind-specific worker
    int n=0;
    SmartPointer<Spectrum::Generic> * sp = &(((gyoto_Spectrum*)obj)->spectrum);
    const string kind = (*sp)->getKind();

    while (n<ygyoto_Spectrum_count && kind.compare(ygyoto_Spectrum_names[n]))
      ++n;

    if (n<ygyoto_Spectrum_count && ygyoto_Spectrum_evals[n]) {
      (*ygyoto_Spectrum_evals[n])(sp, argc);
      return;
    }

    // Fall-back to default worker
    static char * knames[]={
      YGYOTO_SPECTRUM_GENERIC_KW, 0
    };
    static long kglobs[YGYOTO_SPECTRUM_GENERIC_KW_N+1];
    int kiargs[YGYOTO_SPECTRUM_GENERIC_KW_N];
    int piargs[]={-1,-1,-1,-1};
    // push default return value
    *ypush_Spectrum()=*sp;
    yarg_kw_init(knames, kglobs, kiargs);
    
    int iarg=argc, parg=0;
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<4) piargs[parg++]=iarg--;
	else y_error("gyoto_Spectrum takes at most 4 positional arguments");
      }
    }

    int rvset[1]={0}, paUsed[1]={0};
    ygyoto_Spectrum_generic_eval(sp, kiargs, piargs, rvset, paUsed);


  }
  static y_userobj_t gyoto_Spectrum_obj =
    {"gyoto_Spectrum", &gyoto_Spectrum_free, &gyoto_Spectrum_print, &gyoto_Spectrum_eval, 0, 0};

}

SmartPointer<Spectrum::Generic>* yget_Spectrum(int iarg) {
  return &(((gyoto_Spectrum*)yget_obj(iarg, &gyoto_Spectrum_obj))->spectrum);
}

SmartPointer<Spectrum::Generic>* ypush_Spectrum() {
  return &(((gyoto_Spectrum*)(ypush_obj(&gyoto_Spectrum_obj,
					sizeof(gyoto_Spectrum))))->spectrum);
}

int yarg_Spectrum(int iarg) {
  return yget_obj(iarg,0)==gyoto_Spectrum_obj.type_name;
}


void ygyoto_Spectrum_register(char* name, ygyoto_Spectrum_eval_worker_t* on_eval){
  int n;
  if (ygyoto_Spectrum_count==YGYOTO_MAX_REGISTERED)
    y_error("Too many Spectra registered");
  for (n=0; n<ygyoto_Spectrum_count; ++n)
    if (!strcmp(ygyoto_Spectrum_names[n], name)) 
      return;

  strcpy(ygyoto_Spectrum_names[ygyoto_Spectrum_count], name);
  ygyoto_Spectrum_evals[ygyoto_Spectrum_count++]=on_eval;
}

void ygyoto_Spectrum_generic_eval(Gyoto::SmartPointer<Generic>*sp,
				int *kiargs, int *piargs,
				int *rvset, int *paUsed) {
  int k=-1, iarg;
  char * rmsg="Cannot set return value more than once";
  char * pmsg="Cannot use positional argument more than once";

  if (debug())
    for (int i=0; i<YGYOTO_SPECTRUM_GENERIC_KW_N; ++i)
      cerr << "DEBUG: Spectrum_generic_eval: kiargs[" << i << "]="
	   << kiargs[i] << endl;

  // Save to file
  if ((iarg=*(kiargs++))>=0) { // xmlwrite
    iarg+=*rvset;
#ifdef GYOTO_USE_XERCES
    char *filename=ygets_q(iarg);
    Factory(*sp).write(filename);
#else
    y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
  }
    
    // kind
  if ((iarg=*(kiargs++))>=0) {
    iarg+=*rvset;
    if (!yarg_nil(iarg)) y_error("KIND is readonly");
    if (debug()) cerr << "kiargs=" << kiargs << endl;
    if ((*rvset)++) y_error(rmsg);
    char ** kind = ypush_q(0);
    *kind = p_strcpy((*sp)->getKind().c_str());
  }

  /* CLONE */
  if ((iarg=kiargs[++k])>=0) {
    if ((*rvset)++) y_error(rmsg);
    *ypush_Spectrum() = (*sp)->clone();
  }

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
      Inu1nu2[i] = (*sp)->integrate(freqs[i], freqs[i+1]);
  }

  // GET SPECTRUM VALUE FOR WAVELENGTHS
  if (*rvset || *paUsed || (iarg=piargs[0]) < 0 || !yarg_number(iarg)) return;

  if (debug())
    cerr<<"DEBUG: gyoto_Spectrum_generic_eval evaluating Spectrum "
      "at frequency\n";
  long ntot, dims[Y_DIMSIZE];
  double * freqs = ygeta_d(iarg, &ntot, dims);
  double * Inu = ypush_d(dims);
  for (long i=0; i < ntot; ++i) Inu[i] = (**sp)(freqs[i]);

  if (debug()) cerr << "DEBUG: out of Spectrum_generic_eval"<< endl;

}

extern "C" {
  void Y_gyoto_Spectrum(int argc) 
  {
    int rvset[1]={0}, paUsed[1]={0};
    SmartPointer<Spectrum::Generic> *sp = NULL;

    if (yarg_Spectrum(argc-1)) {
      sp = yget_Spectrum(--argc);
      *ypush_Spectrum() = *sp; // push back spectrum
    } else { // Constructor mode
      sp = ypush_Spectrum();
      *rvset=1;
    }

    static char * knames[]={
      YGYOTO_SPECTRUM_GENERIC_KW,
      0
    };
    static long kglobs[YGYOTO_SPECTRUM_GENERIC_KW_N+12];
    int kiargs[YGYOTO_SPECTRUM_GENERIC_KW_N+11];
    int piargs[]={-1,-1,-1,-1};
  
    yarg_kw_init(knames, kglobs, kiargs);
  
    int iarg=argc, parg=0;
    while (iarg>=1) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      if (iarg>=1) {
	if (parg<4) piargs[parg++]=iarg--;
	else y_error("gyoto_Spectrum takes at most 4 positional arguments");
      }
    }

    // if rvset==1, constructor mode:
    if (rvset[1]) {
      if (yarg_string(piargs[0])) {
#ifdef GYOTO_USE_XERCES
	*sp = Factory(ygets_q(piargs[0])).getSpectrum(); 
	paUsed[1]=1;
#else
      y_error("This GYOTO was compiled without XERCES: no xml i/o");
#endif
      } else y_error("Cannot allocate object of virtual class Spectrum");
    }

    ygyoto_Spectrum_generic_eval(sp, kiargs, piargs, rvset, paUsed);

  }

  void
  Y_is_gyoto_Spectrum(int argc)
  {
    ypush_long(yarg_Spectrum(0));
  }

}
