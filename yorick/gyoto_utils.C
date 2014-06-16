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

#include <GyotoUtils.h>
#include <GyotoRegister.h>
#include <yapi.h>
#include <cstring>
#include <iostream>
#include <signal.h>

#include "ygyoto.h"


using namespace std;

static YGyotoSupplier_t *YGyotoGlobalSupplier = NULL;

void ygyotoErrorHandler (const Gyoto::Error e) { y_error(e); }

extern "C" {

  void
  Y_gyoto_dontcatchSIGFPE(int argc)
  {
    signal(SIGFPE, SIG_DFL);
  }


  void
  Y_gyoto_dontcatchSIGSEGV(int argc)
  {
    signal(SIGSEGV, SIG_DFL);
  }


  void
  Y_gyoto_debug(int argc)
  {
    ypush_long(Gyoto::debug());
    if (argc && !yarg_nil(argc)) Gyoto::debug(ygets_l(1));
  }

  void
  Y_gyoto_verbose(int argc)
  {
    ypush_long(Gyoto::verbose());
    if (argc && !yarg_nil(argc)) Gyoto::verbose(ygets_l(1));
  }


  void
  Y_gyoto_loadPlugin(int argc)
  {
    // Step 1: determine whether nofail is set (to true)
    int nofail=0;
    static char const * knames[2] = { "nofail", 0 };
    static long kglobs[2];
    int kiargs[1];
    yarg_kw_init(const_cast<char**>(knames), kglobs, kiargs);
    int iarg=argc-1;
    while (iarg>=0) {
      iarg = yarg_kw(iarg, kglobs, kiargs);
      iarg--;
    }
    if (kiargs[0]>=0) {// nofail= present
      nofail=yarg_true(kiargs[0]);
    }

    // Step 2: load plug-ins
    long ntot=0;
    long dims[Y_DIMSIZE];
    ystring_t * plugins = 0;
    for (int iarg=argc-1; iarg>=0; iarg--) {
      if (kiargs[0]<0 ||(iarg!=kiargs[0] && iarg!=kiargs[0]+1)) {
	plugins = ygeta_q(iarg, &ntot, dims);
	for (long i=0; i<ntot; ++i) Gyoto::loadPlugin(plugins[i], nofail);
      }
    }
    ypush_nil();
    //    Gyoto::Register::init();
  }

  void
  Y___gyoto_initRegister(int argc)
  {
#if defined GYOTO_USE_XERCES
    Gyoto::Register::init();
#endif
  }

  void
  Y_gyoto_listRegister(int argc)
  {
#if defined GYOTO_USE_XERCES
    Gyoto::Register::list();
#endif
  }

  void
  Y_gyoto_haveXerces(int)
  {
    ypush_long(
#if defined GYOTO_USE_XERCES
	       1
#else
	       0
#endif
	       );
  }

  void
  Y_gyoto_haveBoost(int)
  {
    ypush_long(
#if defined HAVE_BOOST
	       1
#else
	       0
#endif
	       );
  }

  void
  Y___gyoto_exportSupplier(int argc)
  {
    if (!YGyotoGlobalSupplier) {
      YGyotoGlobalSupplier = new YGyotoSupplier_t();

      //Plug Metric ABI
      YGyotoGlobalSupplier -> yget_Metric  = &yget_Metric;
      YGyotoGlobalSupplier -> ypush_Metric = &ypush_Metric;
      YGyotoGlobalSupplier -> yarg_Metric  = &yarg_Metric;
      YGyotoGlobalSupplier -> ygyoto_Metric_register = &ygyoto_Metric_register;
      YGyotoGlobalSupplier -> ygyoto_Metric_generic_eval
                                   = &ygyoto_Metric_generic_eval;
      // Plug Astrobj ABI
      YGyotoGlobalSupplier -> yget_Astrobj  = &yget_Astrobj;
      YGyotoGlobalSupplier -> ypush_Astrobj = &ypush_Astrobj;
      YGyotoGlobalSupplier -> yarg_Astrobj  = &yarg_Astrobj;
      YGyotoGlobalSupplier -> ygyoto_Astrobj_register
	                           = &ygyoto_Astrobj_register;
      YGyotoGlobalSupplier -> ygyoto_Astrobj_generic_eval
                                   = &ygyoto_Astrobj_generic_eval;
      YGyotoGlobalSupplier -> ygyoto_ThinDisk_generic_eval
                                   = &ygyoto_ThinDisk_generic_eval;
      // Plug Spectrum ABI
      YGyotoGlobalSupplier -> yget_Spectrum  = &yget_Spectrum;
      YGyotoGlobalSupplier -> ypush_Spectrum = &ypush_Spectrum;
      YGyotoGlobalSupplier -> yarg_Spectrum  = &yarg_Spectrum;
      YGyotoGlobalSupplier -> ygyoto_Spectrum_register
	                           = &ygyoto_Spectrum_register;
      YGyotoGlobalSupplier -> ygyoto_Spectrum_generic_eval
                                   = &ygyoto_Spectrum_generic_eval;
      // Plug Screen ABI
      YGyotoGlobalSupplier -> yget_Screen  = &yget_Screen;
      YGyotoGlobalSupplier -> ypush_Screen = &ypush_Screen;
      YGyotoGlobalSupplier -> yarg_Screen  = &yarg_Screen;

      // Plug Scenery ABI
      YGyotoGlobalSupplier -> yget_Scenery  = &yget_Scenery;
      YGyotoGlobalSupplier -> ypush_Scenery = &ypush_Scenery;
      YGyotoGlobalSupplier -> yarg_Scenery  = &yarg_Scenery;

      // Plug Spectrometer ABI
      YGyotoGlobalSupplier -> yget_Spectrometer  = &yget_Spectrometer;
      YGyotoGlobalSupplier -> ypush_Spectrometer = &ypush_Spectrometer;
      YGyotoGlobalSupplier -> yarg_Spectrometer  = &yarg_Spectrometer;
      YGyotoGlobalSupplier -> ygyoto_Spectrometer_register
	                           = &ygyoto_Spectrometer_register;
      YGyotoGlobalSupplier -> ygyoto_Spectrometer_generic_eval
                                   = &ygyoto_Spectrometer_generic_eval;
    }
    ypush_long(long(YGyotoGlobalSupplier));
  }

  void
  Y___gyoto_setErrorHandler(int argc)
  { Gyoto::Error::setHandler(&ygyotoErrorHandler); }

}
