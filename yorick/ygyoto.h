/*
    Copyright 2011-2015 Thibaut Paumard

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

#ifndef __YGYOTO_H
#define __YGYOTO_H

#include "GyotoConfig.h"
#include "GyotoMetric.h"
#include "GyotoAstrobj.h"
#include "GyotoThinDisk.h"
#include "GyotoSpectrum.h"
#include "GyotoScreen.h"
#include "GyotoPhoton.h"
#include "GyotoScenery.h"
#include "GyotoScreen.h"

#include <cstring>
#include <sstream>

#define YGYOTO_TYPE_LEN 21
#define YGYOTO_STD_CATCH catch(Gyoto::Error e) \
  { y_error(e.get_message().c_str()); }
#define YGYOTO_MAX_REGISTERED 20

typedef void
ygyoto_Metric_eval_worker_t(Gyoto::SmartPointer<Gyoto::Metric::Generic>*, int);
typedef void
ygyoto_Astrobj_eval_worker_t(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>*, int);
typedef void
ygyoto_Spectrum_eval_worker_t(Gyoto::SmartPointer<Gyoto::Spectrum::Generic>*,\
			      int);
typedef void
ygyoto_Spectrometer_eval_worker_t(Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>*, int);

#ifndef YGYOTO_LOCAL_SUPPLIER
/*
  To push or get a Gyoto Object or to check if iarg is a Gyoto object
  of a given type, use the following wrappers around the YAPI
  functions ypush_obj() and yget_obj(). yarg_OBJTYPE return 1 if iarg
  is an object of that type.
 */
Gyoto::SmartPointer<Gyoto::Metric::Generic> *yget_Metric(int iarg);
Gyoto::SmartPointer<Gyoto::Metric::Generic> *ypush_Metric();
int yarg_Metric(int iarg);
Gyoto::SmartPointer<Gyoto::Astrobj::Generic>* yget_Astrobj(int iarg);
Gyoto::SmartPointer<Gyoto::Astrobj::Generic>* ypush_Astrobj();
int yarg_Astrobj(int iarg);
Gyoto::SmartPointer<Gyoto::Screen>* yget_Screen(int iarg);
Gyoto::SmartPointer<Gyoto::Screen>* ypush_Screen();
int yarg_Screen(int iarg);
Gyoto::SmartPointer<Gyoto::Photon>* yget_Photon(int iarg);
Gyoto::SmartPointer<Gyoto::Photon>* ypush_Photon();
int yarg_Photon(int iarg) ;
Gyoto::SmartPointer<Gyoto::Scenery>* yget_Scenery(int iarg);
Gyoto::SmartPointer<Gyoto::Scenery>* ypush_Scenery();
int yarg_Scenery(int iarg);
Gyoto::SmartPointer<Gyoto::Screen>* yget_Screen(int iarg);
Gyoto::SmartPointer<Gyoto::Screen>* ypush_Screen();
int yarg_Screen(int iarg);
Gyoto::SmartPointer<Gyoto::Spectrum::Generic>* yget_Spectrum(int iarg);
Gyoto::SmartPointer<Gyoto::Spectrum::Generic>* ypush_Spectrum();
int yarg_Spectrum(int iarg);
Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>* yget_Spectrometer(int iarg);
Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>* ypush_Spectrometer();
int yarg_Spectrometer(int iarg);
void ypush_property(Gyoto::SmartPointer<Gyoto::SmartPointee>,
		    Gyoto::Property const&, std::string, std::string);
void yget_property(Gyoto::SmartPointer<Gyoto::SmartPointee>,
		   Gyoto::Property const&, int, std::string, std::string);

/*
  You can register your own on_eval worker. It will be called when a
  Metric object GG_OBJ of your specific kind is called as
  gg_obj(arguments). (Resp. Astrobj.) In your own routins, you can
  also call the generic eval (and you certainly should, it's the one
  that processes the generic keywords below). The generic Metric on
  eval must be called at the end of your specific on eval because it
  will try to return the metric coefficients if no other outpur value
  has been pushed.
 */ 
void ygyoto_Metric_register(char const * const kind, ygyoto_Metric_eval_worker_t* on_eval);
void ygyoto_Metric_generic_eval(Gyoto::SmartPointer<Gyoto::Metric::Generic>*,
				int *kiargs, int *piargs, int *rvset,
				int *paUsed, char * unit);


void ygyoto_Astrobj_register(char const * const kind, ygyoto_Astrobj_eval_worker_t* on_eval);
void ygyoto_Astrobj_generic_eval(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>*,
				 int *kiargs, int *piargs, int *rvset,
				 int *paUsed, char * unit);
void ygyoto_ThinDisk_generic_eval(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>*,
				 int *kiargs, int *piargs, int *rvset,
				  int *paUsed, char * unit);


void ygyoto_Spectrum_register(char const * const kind,
			      ygyoto_Spectrum_eval_worker_t*on_eval);
void ygyoto_Spectrum_generic_eval
(Gyoto::SmartPointer<Gyoto::Spectrum::Generic> *,
 int *kiargs, int *piargs, int *rvset, int *paUsed, char *unit);

void ygyoto_Spectrometer_register(char const * const kind,
			      ygyoto_Spectrometer_eval_worker_t*on_eval);
void ygyoto_Spectrometer_generic_eval
(Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> *,
 int *kiargs, int *piargs, int *rvset, int *paUsed, char * unit);


#endif

/*
  The generic on_eval will process these keywords. You need to have
  initiated the kiargs table. You should call the generic on_eval
  after processing your SET keywords and before processing your GET
  keywords. rvset and paUsed are input _and_ output. *rvset == 1 if a
  return value (ypush_*) is already defined, *rvset == 0 if
  not. Likewise, *paUsed means that (at least) one positional argument
  has already been used. Attempting to set the return value twice or
  to use a positional argument twice is to be considered an error
  condition.
*/
// Keywords processed by ygyoto_Metric_generic_eval
#define YGYOTO_METRIC_GENERIC_KW "prime2tdot",				\
    "nullifycoord", "kind", "coordkind", "setparameter", "scalarprod",	\
    "mass", "deltamin", "deltamax", "keplerian",					\
    "unitlength", "circularvelocity",					\
    "christoffel", "xmlwrite", "clone", "help"
// Number of those keywords
#define YGYOTO_METRIC_GENERIC_KW_N 16

// Keywords processed by ygyoto_Astrobj_generic_eval
#define YGYOTO_ASTROBJ_GENERIC_KW "metric", "rmax", "opticallythin",	\
    "xmlwrite", "kind", "setparameter", "clone", "help"
// number of those keywords
#define YGYOTO_ASTROBJ_GENERIC_KW_N 8

// Keywords processed by ygyoto_ThinDisk_generic_eval
#define YGYOTO_THINDISK_GENERIC_KW \
  "innerradius", "outerradius", "thickness", "dir",	\
    YGYOTO_ASTROBJ_GENERIC_KW
// number of those keywords
#define YGYOTO_THINDISK_GENERIC_KW_N 4+YGYOTO_ASTROBJ_GENERIC_KW_N

// maximum number of keywords accepted by a base Astrobj class
#define YGYOTO_ASTROBJ_BASE_MAX_KW_N YGYOTO_THINDISK_GENERIC_KW_N

// Keywords processed by ygyoto_Spectrum_generic_eval
#define YGYOTO_SPECTRUM_GENERIC_KW				\
  "xmlwrite", "kind", "setparameter", "clone", "help", "integrate"
// number of those keywords
#define YGYOTO_SPECTRUM_GENERIC_KW_N 6

// Keywords processed by ygyoto_Spectrometer_generic_eval
#define YGYOTO_SPECTROMETER_GENERIC_KW			\
  "kind", "xmlwrite", "clone", "help", "nsamples",		\
    "setparameter", "channels", "midpoints", "widths"
// number of those keywords
#define YGYOTO_SPECTROMETER_GENERIC_KW_N 9

/*

  The following are neat for writing workers. The worker needs to
  abide by my unwritten coding style. Assume MyKind is a subclass of
  MyBase. MyBase is one of Metric, Astrobj, Spectrum, Spectrometer.

  void ygyoto_MyKind_eval(SmartPointer<Astrobj::Generic>* OBJ_, int argc) {
    // Define keywords
    static char const * knames[]={
      "unit",
      "keyword1", "keyword2", ..., "keywordN",
      YGYOTO_MYBASE_GENERIC_KW,
      0
    };

    // Parse arguments
    YGYOTO_WORKER_INIT(MyBase, MyKind, knames, YGYOTO_MYBASE_GENERIC_KW_N+N+1);

    // Read unit keyword
    YGYOTO_WORKER_SET_UNIT;

    // process individual keywords
    YGYOTO_WORKER_GETSET_DOUBLE_UNIT(Keyword1);
    ...
    YGYOTO_WORKER_GETSET...(KeywordN);

    // Call generic worker
    YGYOTO_WORKER_CALL_GENERIC(MyBase);
  }
 */

#define YGYOTO_STR1(X) #X
#define YGYOTO_STR(X) YGYOTO_STR1(X)
#define YGYOTO_CAT1(X, Y) X##Y
#define YGYOTO_CAT(X, Y) YGYOTO_CAT1(X,Y)

#ifdef GYOTO_USE_XERCES
# define YGYOTO_CONSTRUCTOR_INIT2_XML(BASENAME, GETTER)			\
  char * fname = ygets_q(argc-1);					\
  OBJ = ypush_##BASENAME();						\
  GYOTO_DEBUG_EXPR(OBJ);						\
  * OBJ = Gyoto::Factory(fname).GETTER();				\
  GYOTO_DEBUG << "Swapping object for filename\n";	                \
  yarg_swap(0, argc);			\
  GYOTO_DEBUG << "Dropping filename from stack\n";			\
  yarg_drop(1);							\
  GYOTO_DEBUG << "Dropped filename from stack\n";			\
  --argc;
#else
# define YGYOTO_CONSTRUCTOR_INIT2_XML(BASENAME, GETTER)	\
  y_error("this gyoto was built without XML support");
#endif


#define YGYOTO_CONSTRUCTOR_INIT2(BASENAME, BASECLASS, DERIVEDCLASS, GETTER)\
  Gyoto::SmartPointer<BASECLASS> *OBJ = NULL;				\
  if (yarg_##BASENAME(argc-1)) {					\
    OBJ = yget_##BASENAME(--argc);					\
    GYOTO_DEBUG_EXPR(OBJ);						\
  } else if (yarg_string(argc-1)) {					\
    YGYOTO_CONSTRUCTOR_INIT2_XML(BASENAME, GETTER);			\
  } else {								\
    OBJ = ypush_##BASENAME();						\
    GYOTO_DEBUG_EXPR(OBJ);						\
    *OBJ = new DERIVEDCLASS();						\
    GYOTO_DEBUG << "object created" << endl;				\
    for (int arg=0; arg<argc; ++arg)					\
      yarg_swap(arg, arg+1);						\
  }									\
  if (argc==1 && yarg_nil(0)) {						\
    yarg_drop(1);							\
    --argc;								\
  }
#define YGYOTO_CONSTRUCTOR_INIT1(BASENAME, BASECLASS, DERIVEDCLASS)    \
  YGYOTO_CONSTRUCTOR_INIT2(BASENAME, BASECLASS, DERIVEDCLASS, get##BASENAME)
#define YGYOTO_CONSTRUCTOR_INIT(BASE, KIND) \
  YGYOTO_CONSTRUCTOR_INIT1(BASE, Gyoto::BASE::Generic, Gyoto::BASE::KIND)

#define YGYOTO_WORKER_INIT(BASE, KIND, KNAMES, NKW)			\
  YGYOTO_WORKER_INIT1(BASE, Gyoto::BASE::KIND, KNAMES, NKW)
#define YGYOTO_WORKER_INIT1(BASE, CLASS, KNAMES, NKW)			\
  int rvset[1]={0}, paUsed[1]={0};					\
  *ypush_##BASE()=*YGYOTO_CAT(OBJ, _);					\
  {									\
  long kidx;								\
  Gyoto::Property const * prop;						\
  std::string pname="", unit="";					\
  bool unit_found, prop_pushed=false;					\
  ++argc;								\
  while ( (  argc > 0                 ) &&				\
	  ( (kidx=yarg_key(argc-1)) >=0 ) &&				\
	  ( (prop=(*YGYOTO_CAT(OBJ, _))					\
	     ->property(pname=yfind_name(kidx)))) ) {			\
    if ( (kidx=yarg_key(argc-3)) >=0 &&					\
	 !strcmp(yfind_name(kidx),"unit") ) {				\
      unit=ygets_q(argc-4);						\
      unit_found=true;							\
    } else {								\
      unit_found=false;							\
      unit="";								\
    }									\
    if (yarg_nil(argc-2)) {						\
      if (prop_pushed) y_error("Can push only one return value");	\
      prop_pushed=true;							\
      yarg_drop(1);							\
      ypush_property(*YGYOTO_CAT(OBJ, _), *prop, pname, unit);	\
    }									\
    else yget_property(*YGYOTO_CAT(OBJ, _), *prop, argc-2, pname, unit); \
    argc -= unit_found?4:2;						\
  }									\
  if (prop_pushed) ++*rvset;						\
  --argc;								\
  }									\
  Gyoto::SmartPointer<CLASS> *OBJ =					\
    (Gyoto::SmartPointer<CLASS> *)YGYOTO_CAT(OBJ, _);			\
  static long kglobs[NKW+1];						\
  int kiargs[NKW];							\
  int piargs[]={-1,-1,-1,-1,-1};					\
  yarg_kw_init(const_cast<char**>(KNAMES), kglobs, kiargs);		\
  int iarg=argc, parg=0;						\
  while (iarg>=1) {							\
    iarg = yarg_kw(iarg, kglobs, kiargs);				\
    if (iarg>=1) {							\
      if (parg<5) piargs[parg++]=iarg--;				\
      else y_error( #CLASS " worker takes at most 5 positional arguments"); \
    }									\
  }									\
  GYOTO_DEBUG_ARRAY(piargs, 5);						\
  GYOTO_DEBUG_ARRAY(kiargs, NKW);					\
  int k=-1;								\
  char const * rmsg="Cannot set return value more than once";		\
  char const * pmsg="Cannot use positional argument more than once";	\
  char * unit=NULL;							\
  /* make sure we use the variables at least once */			\
  if (OBJ && rmsg && pmsg && unit && k && *paUsed) {}


#define YGYOTO_WORKER_GETSET_VECTOR(MEMBER, N)			  \
  if ((iarg=kiargs[++k])>=0) {					  \
    iarg+=*rvset;						  \
    if (yarg_nil(iarg)) {					  \
      if ((*rvset)++) y_error(rmsg);				  \
      long dims[] = {1,N};					  \
      double * coord=ypush_d(dims);				  \
      (*OBJ)-> get##MEMBER (coord);				  \
    } else {							  \
      long ntot;						  \
      double * pos = ygeta_d(iarg, &ntot, NULL);		  \
      if (ntot<N) y_error("POS must have at least " #N " elements");	  \
      (*OBJ) -> set##MEMBER (pos);				  \
    }								  \
  }

#define YGYOTO_WORKER_GETSET4(MEMBER)				  \
  YGYOTO_WORKER_GETSET_VECTOR(MEMBER, 4)

#define YGYOTO_WORKER_GETSET_DOUBLE_UNIT(MEMBER)			   \
  if ((iarg=kiargs[++k])>=0) {					   \
      iarg+=*rvset;						   \
      if (yarg_nil(iarg)) {					   \
	if ((*rvset)++) y_error("Only one return value possible"); \
	ypush_double((*OBJ) -> get##MEMBER (unit?unit:""));	   \
      } else							   \
	(*OBJ) -> set##MEMBER (ygets_d(iarg), unit?unit:"");	   \
    }

#define YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(MEMBER)			   \
  if ((iarg=kiargs[++k])>=0) {					   \
      iarg+=*rvset;						   \
      if (yarg_nil(iarg)) {					   \
	if ((*rvset)++) y_error("Only one return value possible"); \
	ypush_double((*OBJ) -> MEMBER (unit?unit:""));	   \
      } else							   \
	(*OBJ) -> MEMBER (ygets_d(iarg), unit?unit:"");	   \
    }

#define YGYOTO_WORKER_GETSET_DOUBLE(MEMBER)			   \
  if ((iarg=kiargs[++k])>=0) {					   \
      iarg+=*rvset;						   \
      if (yarg_nil(iarg)) {					   \
	if ((*rvset)++) y_error("Only one return value possible"); \
	ypush_double((*OBJ) -> get##MEMBER ());		   \
      } else							   \
	(*OBJ) -> set##MEMBER (ygets_d(iarg));		   \
    }

#define YGYOTO_WORKER_GET_DOUBLE(MEMBER)			   \
  if ((iarg=kiargs[++k])>=0) {					   \
      iarg+=*rvset;						   \
      if (yarg_nil(iarg)) {					   \
	if ((*rvset)++) y_error("Only one return value possible"); \
	ypush_double((*OBJ) -> MEMBER ());		   \
      } else							   \
	y_error( #MEMBER " is read only"); \
    }

#define YGYOTO_WORKER_GETSET_DOUBLE2(MEMBER)			   \
  if ((iarg=kiargs[++k])>=0) {					   \
      iarg+=*rvset;						   \
      if (yarg_nil(iarg)) {					   \
	if ((*rvset)++) y_error("Only one return value possible"); \
	ypush_double((*OBJ) -> MEMBER ());			   \
      } else							   \
	(*OBJ) -> MEMBER (ygets_d(iarg));			   \
    }

#define YGYOTO_WORKER_GETSET_STRING2(MEMBER)			   \
  if ((iarg=kiargs[++k])>=0) {					   \
      iarg+=*rvset;						   \
      if (yarg_nil(iarg)) {					   \
	if ((*rvset)++) y_error("Only one return value possible"); \
	*ypush_q(0) = p_strcpy((*OBJ) -> MEMBER () . c_str() );	   \
      } else							   \
	(*OBJ) -> MEMBER (std::string(ygets_q(iarg)));		   \
    }

#define YGYOTO_WORKER_GETSET_LONG(MEMBER)			   \
  if ((iarg=kiargs[++k])>=0) {					   \
      iarg+=*rvset;						   \
      if (yarg_nil(iarg)) {					   \
	if ((*rvset)++) y_error("Only one return value possible"); \
	ypush_long((*OBJ) -> get##MEMBER ());		   \
      } else							   \
	(*OBJ) -> set##MEMBER (ygets_l(iarg));		   \
    }

#define YGYOTO_WORKER_GETSET_LONG2(MEMBER)			   \
  if ((iarg=kiargs[++k])>=0) {					   \
      iarg+=*rvset;						   \
      if (yarg_nil(iarg)) {					   \
	if ((*rvset)++) y_error("Only one return value possible"); \
	ypush_long((*OBJ) -> MEMBER ());		   \
      } else							   \
	(*OBJ) -> MEMBER (ygets_l(iarg));		   \
    }

#define YGYOTO_WORKER_GETSET_OBJECT0(SETMEMBER,GETMEMBER,YMEMBER)   \
  if ((iarg=kiargs[++k])>=0) {					   \
    iarg+=*rvset;						   \
    if (yarg_nil(iarg)) {					   \
      GYOTO_DEBUG << "pushing " #YMEMBER << std::endl;		   \
      if ((*rvset)++) y_error("Only one return value possible");   \
      void * tmp = (*OBJ) -> GETMEMBER ();			   \
      if (tmp)							   \
	*ypush_##YMEMBER () = (*OBJ) -> GETMEMBER ();		   \
      else ypush_long(0);					   \
    } else {							   \
      GYOTO_DEBUG << "setting " #YMEMBER << std::endl;		   \
      if (yarg_number(iarg)) {					   \
        if (ygets_l(iarg) != 0)					   \
	  y_error("Argument should be a " #YMEMBER " or 0");	   \
	(*OBJ) -> SETMEMBER (0);				   \
      } else							   \
	(*OBJ) -> SETMEMBER (*yget_##YMEMBER (kiargs[k]));	   \
    }								   \
  }

#define YGYOTO_WORKER_GETSET_OBJECT(MEMBER)			   \
  YGYOTO_WORKER_GETSET_OBJECT0(set##MEMBER,get##MEMBER,MEMBER)

#define YGYOTO_WORKER_GETSET_OBJECT2(MEMBER,YMEMBER)		   \
  YGYOTO_WORKER_GETSET_OBJECT0(MEMBER,MEMBER,YMEMBER)

#define YGYOTO_WORKER_SET_UNIT		 \
  if ((iarg=kiargs[++k])>=0) {		 \
    iarg+=*rvset;			 \
    GYOTO_DEBUG << "set unit" << std::endl;	\
    unit = ygets_q(iarg);		 \
  }

#define YGYOTO_WORKER_SETPARAMETER			\
  if ((iarg=kiargs[++k])>=0) {				\
    iarg+=*rvset; 					\
    if ((*paUsed)++) y_error("pmsg");			\
    string name = ygets_q(iarg);			\
    string content = "";				\
    if (piargs[0] >= 0) content = ygets_q(*piargs);		\
    try	{							\
      (*OBJ)->setParameter(name, content,  unit?unit:"");	\
    } YGYOTO_STD_CATCH;						\
  }

#ifdef GYOTO_USE_XERCES
# define YGYOTO_WORKER_XMLWRITE			\
  if ((iarg=kiargs[++k])>=0) {			\
    iarg+=*rvset;				\
    char *filename=ygets_q(iarg);		\
    Gyoto::Factory(*OBJ).write(filename);	\
  }
#else
# define YGYOTO_WORKER_XMLWRITE						\
  if ((iarg=kiargs[++k])>=0) {						\
    y_error("This GYOTO was compiled without xerces: no xml i/o");	\
  }
#endif

#define YGYOTO_WORKER_CLONE(TYPE)		\
  if ((iarg=kiargs[++k])>=0) {			\
    if ((*rvset)++) y_error(rmsg);		\
    *ypush_##TYPE () = (*OBJ)->clone();		\
  }

#define YGYOTO_WORKER_RUN(METHOD)			\
  if ((iarg=kiargs[++k])>=0) {				\
    GYOTO_DEBUG << #METHOD << std::endl ;		\
    iarg+=*rvset;					\
    (*OBJ)-> METHOD ;			\
  }

#define YGYOTO_WORKER_HELP				\
  if ((iarg=kiargs[++k])>=0) {				\
    GYOTO_DEBUG << "help" << std::endl ;		\
    iarg+=*rvset;					\
    (*OBJ)-> help() ;					\
    ypush_nil();					\
  }

#define YGYOTO_WORKER_CALL_GENERIC(BASE) \
  ygyoto_##BASE##_generic_eval(YGYOTO_CAT(OBJ,_), \
			       kiargs+k+1, piargs, rvset, paUsed, unit);

/*

  The following are needed to export the ABI to other plug-ins.

  You usually don't need to read below this line.

 */

typedef Gyoto::SmartPointer<Gyoto::Metric::Generic> *ygyoto_yget_Metric_t(int);
typedef Gyoto::SmartPointer<Gyoto::Metric::Generic> *ygyoto_ypush_Metric_t();

typedef void ygyoto_ypush_property_t
(Gyoto::SmartPointer<Gyoto::SmartPointee>,
 Gyoto::Property const&, std::string, std::string);
typedef void ygyoto_yget_property_t
(Gyoto::SmartPointer<Gyoto::SmartPointee>,
 Gyoto::Property const&, int, std::string, std::string);

typedef int yarg_OBJTYPE_t(int);
typedef void ygyoto_Metric_register_t(char const * const, ygyoto_Metric_eval_worker_t*);
typedef void ygyoto_Metric_generic_eval_t(Gyoto::SmartPointer<Gyoto::Metric::Generic>*, \
					  int *, int *, int*, int*, char *);

typedef Gyoto::SmartPointer<Gyoto::Astrobj::Generic> *ygyoto_yget_Astrobj_t(int);
typedef Gyoto::SmartPointer<Gyoto::Astrobj::Generic> *ygyoto_ypush_Astrobj_t();
//typedef int yarg_Astrobj_t(int);
typedef void ygyoto_Astrobj_register_t(char const * const, ygyoto_Astrobj_eval_worker_t*);
typedef void ygyoto_Astrobj_generic_eval_t \
(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>*, int *, int *, int*, int*, char *);

typedef Gyoto::SmartPointer<Gyoto::Spectrum::Generic> *\
ygyoto_yget_Spectrum_t(int);
typedef Gyoto::SmartPointer<Gyoto::Spectrum::Generic> *\
ygyoto_ypush_Spectrum_t();
//typedef int yarg_Spectrum::Generic_t(int);
typedef void ygyoto_Spectrum_register_t\
(char const*const, ygyoto_Spectrum_eval_worker_t*);
typedef void ygyoto_Spectrum_generic_eval_t \
(Gyoto::SmartPointer<Gyoto::Spectrum::Generic>*, int *, int *, int*, int*, char*);

typedef Gyoto::SmartPointer<Gyoto::Screen> *ygyoto_yget_Screen_t(int);
typedef Gyoto::SmartPointer<Gyoto::Screen> *ygyoto_ypush_Screen_t();

typedef Gyoto::SmartPointer<Gyoto::Scenery> *ygyoto_yget_Scenery_t(int);
typedef Gyoto::SmartPointer<Gyoto::Scenery> *ygyoto_ypush_Scenery_t();

typedef Gyoto::SmartPointer<Gyoto::Photon> *ygyoto_yget_Photon_t(int);
typedef Gyoto::SmartPointer<Gyoto::Photon> *ygyoto_ypush_Photon_t();

typedef Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> *\
ygyoto_yget_Spectrometer_t(int);
typedef Gyoto::SmartPointer<Gyoto::Spectrometer::Generic> *\
ygyoto_ypush_Spectrometer_t();
//typedef int yarg_Spectrometer::Generic_t(int);
typedef void ygyoto_Spectrometer_register_t\
(char const*const, ygyoto_Spectrometer_eval_worker_t*);
typedef void ygyoto_Spectrometer_generic_eval_t \
(Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>*, int *, int *, int*, int*, char*);

typedef struct YGyotoSupplier {
  // Metric
  ygyoto_yget_Metric_t  *yget_Metric;
  ygyoto_ypush_Metric_t *ypush_Metric;
  yarg_OBJTYPE_t         *yarg_Metric;
  ygyoto_Metric_register_t *ygyoto_Metric_register;
  ygyoto_Metric_generic_eval_t* ygyoto_Metric_generic_eval;
  // Astrobj
  ygyoto_yget_Astrobj_t  *yget_Astrobj;
  ygyoto_ypush_Astrobj_t *ypush_Astrobj;
  yarg_OBJTYPE_t         *yarg_Astrobj;
  ygyoto_Astrobj_register_t *ygyoto_Astrobj_register;
  ygyoto_Astrobj_generic_eval_t* ygyoto_Astrobj_generic_eval;
  ygyoto_Astrobj_generic_eval_t* ygyoto_ThinDisk_generic_eval;
  // Spectrum
  ygyoto_yget_Spectrum_t  *yget_Spectrum;
  ygyoto_ypush_Spectrum_t *ypush_Spectrum;
  yarg_OBJTYPE_t         *yarg_Spectrum;
  ygyoto_Spectrum_register_t *ygyoto_Spectrum_register;
  ygyoto_Spectrum_generic_eval_t* ygyoto_Spectrum_generic_eval;
  // Screen
  ygyoto_yget_Screen_t  *yget_Screen;
  ygyoto_ypush_Screen_t *ypush_Screen;
  yarg_OBJTYPE_t        *yarg_Screen;
  // Scenery
  ygyoto_yget_Scenery_t  *yget_Scenery;
  ygyoto_ypush_Scenery_t *ypush_Scenery;
  yarg_OBJTYPE_t         *yarg_Scenery;
  // Photon
  ygyoto_yget_Photon_t  *yget_Photon;
  ygyoto_ypush_Photon_t *ypush_Photon;
  yarg_OBJTYPE_t        *yarg_Photon;
  // Spectrometer
  ygyoto_yget_Spectrometer_t  *yget_Spectrometer;
  ygyoto_ypush_Spectrometer_t *ypush_Spectrometer;
  yarg_OBJTYPE_t        *yarg_Spectrometer;
  ygyoto_Spectrometer_register_t *ygyoto_Spectrometer_register;
  ygyoto_Spectrometer_generic_eval_t* ygyoto_Spectrometer_generic_eval;
  // Properties
  ygyoto_ypush_property_t *ypush_property;
  ygyoto_yget_property_t *yget_property;
} YGyotoSupplier_t;


#ifdef YGYOTO_LOCAL_SUPPLIER
// The above ABI is exposed as macros in external plug-ins

extern YGyotoSupplier_t* YGYOTO_LOCAL_SUPPLIER;

#define yget_Metric(iarg) YGYOTO_LOCAL_SUPPLIER -> yget_Metric(iarg)
#define ypush_Metric()   YGYOTO_LOCAL_SUPPLIER  -> ypush_Metric()
#define yget_property(a, b, c, d, e) \
  YGYOTO_LOCAL_SUPPLIER -> yget_property(a, b, c, d, e)
#define ypush_property(a, b, c, d)\
  YGYOTO_LOCAL_SUPPLIER  -> ypush_property(a, b, c, d)
#define yarg_Metric(iarg) YGYOTO_LOCAL_SUPPLIER -> yarg_Metric(iarg)
#define ygyoto_Metric_register(kind, on_eval) \
                          YGYOTO_LOCAL_SUPPLIER -> \
		  ygyoto_Metric_register(kind, on_eval)
#define ygyoto_Metric_generic_eval(gg, kiargs, piargs, rvset, paUsed, unit) \
        YGYOTO_LOCAL_SUPPLIER -> \
        ygyoto_Metric_generic_eval(gg, kiargs, piargs, rvset, paUsed, unit)

#define yget_Astrobj(iarg) YGYOTO_LOCAL_SUPPLIER -> yget_Astrobj(iarg)
#define ypush_Astrobj()   YGYOTO_LOCAL_SUPPLIER  -> ypush_Astrobj()
#define yarg_Astrobj(iarg) YGYOTO_LOCAL_SUPPLIER -> yarg_Astrobj(iarg)
#define ygyoto_Astrobj_register(kind, on_eval) \
                          YGYOTO_LOCAL_SUPPLIER -> \
		  ygyoto_Astrobj_register(kind, on_eval)
#define ygyoto_Astrobj_generic_eval(gg, kiargs, piargs, rvset, paUsed, unit) \
                          YGYOTO_LOCAL_SUPPLIER -> \
			  ygyoto_Astrobj_generic_eval(gg, kiargs, piargs, rvset, paUsed, unit)
#define ygyoto_ThinDisk_generic_eval(gg, kiargs, piargs, rvset, paUsed, unit) \
                          YGYOTO_LOCAL_SUPPLIER -> \
			  ygyoto_ThinDisk_generic_eval(gg, kiargs, piargs, rvset, paUsed, unit)

#define yget_Spectrum(iarg) YGYOTO_LOCAL_SUPPLIER -> yget_Spectrum(iarg)
#define ypush_Spectrum()   YGYOTO_LOCAL_SUPPLIER  -> ypush_Spectrum()
#define yarg_Spectrum(iarg) YGYOTO_LOCAL_SUPPLIER -> yarg_Spectrum(iarg)
#define ygyoto_Spectrum_register(kind, on_eval) \
                          YGYOTO_LOCAL_SUPPLIER -> \
		  ygyoto_Spectrum_register(kind, on_eval)
#define ygyoto_Spectrum_generic_eval(gg, kiargs, piargs, rvset, paUsed, unit) \
                          YGYOTO_LOCAL_SUPPLIER -> \
			  ygyoto_Spectrum_generic_eval(gg, kiargs, piargs, rvset, paUsed, unit)

#define yget_Screen(iarg) YGYOTO_LOCAL_SUPPLIER -> yget_Screen(iarg)
#define ypush_Screen()   YGYOTO_LOCAL_SUPPLIER  -> ypush_Screen()
#define yarg_Screen(iarg) YGYOTO_LOCAL_SUPPLIER -> yarg_Screen(iarg)

#define yget_Scenery(iarg) YGYOTO_LOCAL_SUPPLIER -> yget_Scenery(iarg)
#define ypush_Scenery()   YGYOTO_LOCAL_SUPPLIER  -> ypush_Scenery()
#define yarg_Scenery(iarg) YGYOTO_LOCAL_SUPPLIER -> yarg_Scenery(iarg)

#define yget_Photon(iarg) YGYOTO_LOCAL_SUPPLIER -> yget_Photon(iarg)
#define ypush_Photon()   YGYOTO_LOCAL_SUPPLIER  -> ypush_Photon()
#define yarg_Photon(iarg) YGYOTO_LOCAL_SUPPLIER -> yarg_Photon(iarg)

#define yget_Spectrometer(iarg) YGYOTO_LOCAL_SUPPLIER -> yget_Spectrometer(iarg)
#define ypush_Spectrometer()   YGYOTO_LOCAL_SUPPLIER  -> ypush_Spectrometer()
#define yarg_Spectrometer(iarg) YGYOTO_LOCAL_SUPPLIER -> yarg_Spectrometer(iarg)
#define ygyoto_Spectrometer_register(kind, on_eval) \
                          YGYOTO_LOCAL_SUPPLIER -> \
		  ygyoto_Spectrometer_register(kind, on_eval)
#define ygyoto_Spectrometer_generic_eval(gg, kiargs, piargs, rvset, paUsed) \
                          YGYOTO_LOCAL_SUPPLIER -> \
		  ygyoto_Spectrometer_generic_eval(gg, kiargs, piargs, rvset,paUsed)

#endif


#endif
