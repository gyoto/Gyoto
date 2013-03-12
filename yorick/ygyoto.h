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
 int *kiargs, int *piargs, int *rvset, int *paUsed);

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
#define YGYOTO_METRIC_GENERIC_KW "prime2tdot", \
    "nullifycoord", "kind", "setparameter", 				\
    "mass", "unitlength", "circularvelocity", "xmlwrite", "clone"
// Number of those keywords
#define YGYOTO_METRIC_GENERIC_KW_N 9

// Keywords processed by ygyoto_Astrobj_generic_eval
#define YGYOTO_ASTROBJ_GENERIC_KW "metric", "rmax", "opticallythin",	\
    "xmlwrite", "kind", "setparameter", "clone"
// number of those keywords
#define YGYOTO_ASTROBJ_GENERIC_KW_N 7

// Keywords processed by ygyoto_ThinDisk_generic_eval
#define YGYOTO_THINDISK_GENERIC_KW \
  "innerradius", "outerradius", "thickness", "dir",	\
    YGYOTO_ASTROBJ_GENERIC_KW
// number of those keywords
#define YGYOTO_THINDISK_GENERIC_KW_N 4+YGYOTO_ASTROBJ_GENERIC_KW_N

// maximum number of keywords accepted by a base Astrobj class
#define YGYOTO_ASTROBJ_BASE_MAX_KW_N YGYOTO_THINDISK_GENERIC_KW_N

// Keywords processed by ygyoto_Spectrum_generic_eval
#define YGYOTO_SPECTRUM_GENERIC_KW "xmlwrite", "kind", "setparameter", "clone", "integrate"
// number of those keywords
#define YGYOTO_SPECTRUM_GENERIC_KW_N 5

// Keywords processed by ygyoto_Spectrometer_generic_eval
#define YGYOTO_SPECTROMETER_GENERIC_KW "kind", "xmlwrite", "clone", "nsamples",\
    "channels", "midpoints", "widths"
// number of those keywords
#define YGYOTO_SPECTROMETER_GENERIC_KW_N 7

/*

  The following are needed to export the ABI to other plug-ins.

 */

typedef Gyoto::SmartPointer<Gyoto::Metric::Generic> *ygyoto_yget_Metric_t(int);
typedef Gyoto::SmartPointer<Gyoto::Metric::Generic> *ygyoto_ypush_Metric_t();
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
(Gyoto::SmartPointer<Gyoto::Spectrum::Generic>*, int *, int *, int*, int*);

typedef Gyoto::SmartPointer<Gyoto::Screen> *ygyoto_yget_Screen_t(int);
typedef Gyoto::SmartPointer<Gyoto::Screen> *ygyoto_ypush_Screen_t();

typedef Gyoto::SmartPointer<Gyoto::Scenery> *ygyoto_yget_Scenery_t(int);
typedef Gyoto::SmartPointer<Gyoto::Scenery> *ygyoto_ypush_Scenery_t();

typedef Gyoto::SmartPointer<Gyoto::Photon> *ygyoto_yget_Photon_t(int);
typedef Gyoto::SmartPointer<Gyoto::Photon> *ygyoto_ypush_Photon_t();

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
} YGyotoSupplier_t;


#ifdef YGYOTO_LOCAL_SUPPLIER
// The above ABI is exposed as macros in external plug-ins

extern YGyotoSupplier_t* YGYOTO_LOCAL_SUPPLIER;

#define yget_Metric(iarg) YGYOTO_LOCAL_SUPPLIER -> yget_Metric(iarg)
#define ypush_Metric()   YGYOTO_LOCAL_SUPPLIER  -> ypush_Metric()
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
#define ygyoto_Spectrum_generic_eval(gg, kiargs, piargs, rvset, paUsed) \
                          YGYOTO_LOCAL_SUPPLIER -> \
		  ygyoto_Spectrum_generic_eval(gg, kiargs, piargs, rvset,paUsed)

#define yget_Screen(iarg) YGYOTO_LOCAL_SUPPLIER -> yget_Screen(iarg)
#define ypush_Screen()   YGYOTO_LOCAL_SUPPLIER  -> ypush_Screen()
#define yarg_Screen(iarg) YGYOTO_LOCAL_SUPPLIER -> yarg_Screen(iarg)

#define yget_Scenery(iarg) YGYOTO_LOCAL_SUPPLIER -> yget_Scenery(iarg)
#define ypush_Scenery()   YGYOTO_LOCAL_SUPPLIER  -> ypush_Scenery()
#define yarg_Scenery(iarg) YGYOTO_LOCAL_SUPPLIER -> yarg_Scenery(iarg)

#define yget_Photon(iarg) YGYOTO_LOCAL_SUPPLIER -> yget_Photon(iarg)
#define ypush_Photon()   YGYOTO_LOCAL_SUPPLIER  -> ypush_Photon()
#define yarg_Photon(iarg) YGYOTO_LOCAL_SUPPLIER -> yarg_Photon(iarg)

#endif


#endif
