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

#include "GyotoUtils.h"
#include "GyotoRegister.h"
#include "GyotoProperty.h"
#include "GyotoValue.h"
#include "GyotoObject.h"
#include "GyotoMetric.h"
#include "GyotoAstrobj.h"
#include "GyotoSpectrum.h"
#include "GyotoSpectrometer.h"
#include "GyotoScreen.h"
#include <yapi.h>
#include <pstdlib.h>
#include <cstring>
#include <iostream>
#include <sstream>
#include <signal.h>
#include <vector>

#if defined HAVE_FENV_H
# include <fenv.h>
#endif

#if defined HAVE_MPI
# include <mpi.h>
#endif

#include "ygyoto.h"
#include "ygyoto_private.h"


using namespace std;

static YGyotoSupplier_t *YGyotoGlobalSupplier = NULL;

void ygyotoErrorHandler (const Gyoto::Error e) { y_error(e); }

#if defined HAVE_MPI
void ygyotoMPIErrorHandlerFcn(MPI_Comm * comm, int * error_code, ...) {
  char error_string[MPI_MAX_ERROR_STRING];
  int error_string_length;
  MPI_Error_string(*error_code, error_string, &error_string_length);
  error_string[error_string_length] = '\0';
  y_error(error_string);
}
MPI_Errhandler ygyotoMPIErrorHandler;
#endif

extern "C" {

  void
  Y_gyoto_dontcatchSIGFPE(int)
  {
    signal(SIGFPE, SIG_DFL);
  }

  void Y_gyoto_FE(int argc) {
    std::string name=ygets_q(0);
#if defined HAVE_FENV_H
# if defined FE_DIVBYZERO
    if (name=="DIVBYZERO") {
      ypush_int(FE_DIVBYZERO);
      return;
    }
# endif
# if defined FE_INEXACT
    if (name=="INEXACT") {
      ypush_int(FE_INEXACT);
      return;
    }
# endif
# if defined FE_INVALID
    if (name=="INVALID") {
      ypush_int(FE_INVALID);
      return;
    }
# endif
# if defined FE_OVERFLOW
    if (name=="OVERFLOW") {
      ypush_int(FE_OVERFLOW);
      return;
    }
# endif
# if defined FE_UNDERFLOW
    if (name=="UNDERFLOW") {
      ypush_int(FE_UNDERFLOW);
      return;
    }
# endif
# if defined FE_ALL_EXCEPT
    if (name=="ALL_EXCEPT") {
      ypush_int(FE_ALL_EXCEPT);
      return;
    }
# endif
    y_errorq("No such exception: FE_%s", name.c_str());
#else
    GYOTO_WARNING << "no GNU fenv.h in this Gyoto\n";
#endif
  }

  void
  Y_gyoto_dontcatchSIGSEGV(int)
  {
    signal(SIGSEGV, SIG_DFL);
  }

  void
  Y_gyoto_fedisableexcept(int argc)
  {
#if defined HAVE_FENV_H
    int excepts=FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID;
    if (argc && !yarg_nil(0)) excepts=ygets_i(0);
    ypush_int(fedisableexcept(excepts));
#else
    GYOTO_WARNING << "no GNU fenv.h in this Gyoto\n";
#endif
  }

  void
  Y_gyoto_feenableexcept(int argc)
  {
#if defined HAVE_FENV_H
    int excepts=FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID;
    if (argc && !yarg_nil(0)) excepts=ygets_i(0);
    ypush_int(feenableexcept(excepts));
#else
    GYOTO_WARNING << "no GNU fenv.h in this Gyoto\n";
#endif
  }

  void
  Y_gyoto_debug(int argc)
  {
    ypush_long(Gyoto::debug());
    if (argc && !yarg_nil(argc)) Gyoto::debug(int(ygets_l(1)));
  }

  void
  Y_gyoto_verbose(int argc)
  {
    ypush_long(Gyoto::verbose());
    if (argc && !yarg_nil(argc)) Gyoto::verbose(int(ygets_l(1)));
  }


  void
  Y_gyoto_havePlugin(int argc)
  {
    ypush_long(Gyoto::havePlugin(ygets_q(0)));
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
    for (iarg=argc-1; iarg>=0; iarg--) {
      if (kiargs[0]<0 ||(iarg!=kiargs[0] && iarg!=kiargs[0]+1)) {
	plugins = ygeta_q(iarg, &ntot, dims);
	for (long i=0; i<ntot; ++i) Gyoto::loadPlugin(plugins[i], nofail);
      }
    }
    ypush_nil();
    //    Gyoto::Register::init();
  }

  void
  Y_gyoto_requirePlugin(int argc)
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
    for (iarg=argc-1; iarg>=0; iarg--) {
      if (kiargs[0]<0 ||(iarg!=kiargs[0] && iarg!=kiargs[0]+1)) {
	plugins = ygeta_q(iarg, &ntot, dims);
	for (long i=0; i<ntot; ++i) Gyoto::requirePlugin(plugins[i], nofail);
      }
    }
    ypush_nil();
    //    Gyoto::Register::init();
  }

  void
  Y___gyoto_initRegister(int argc)
  {
    const char * pluglist = NULL;
    if (argc && !yarg_nil(argc-1)) pluglist=ygets_q(argc-1);
    Gyoto::Register::init(pluglist);
  }

  void
  Y_gyoto_listRegister(int)
  {
    Gyoto::Register::list();
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
  Y_gyoto_haveUDUNITS(int)
  {
    ypush_long(
#if defined HAVE_UDUNITS
	       1
#else
	       0
#endif
	       );
  }

  void
  Y_gyoto_haveCFITSIO(int)
  {
    ypush_long(
#if defined GYOTO_USE_CFITSIO
	       1
#else
	       0
#endif
	       );
  }

  void
  Y_gyoto_havePTHREAD(int)
  {
    ypush_long(
#if defined HAVE_PTHREAD
	       1
#else
	       0
#endif
	       );
  }

  void
  Y_gyoto_haveMPI(int)
  {
    ypush_long(
#if defined HAVE_MPI
	       1
#else
	       0
#endif
	       );
  }

  void
  Y_gyoto_haveFENV(int)
  {
    ypush_long(
#if defined HAVE_FENV_H
	       1
#else
	       0
#endif
	       );
  }

  void
  Y_gyoto_MPI_Init(int argc)
  {
#if defined HAVE_MPI
    long int mpiargcl=0;
    char **mpiargv=NULL;
    long index=-1;
    if (argc>1) y_error("gyoto.MPI_Init() takes at most one argument");
    if (argc) {
      index=yget_ref(0);
      if (!yarg_nil(0)) mpiargv=ygeta_q(0, &mpiargcl, NULL);
    }
    int mpiargc=mpiargcl;
    ypush_long(MPI_Init(&mpiargc, &mpiargv));
    if (index>=0) {
      long dims[]={1, mpiargc};
      ystring_t * out=ypush_q(dims);
      for (long i=0; i<mpiargc; ++i) out[i] = p_strcpy(mpiargv[i]);
      yput_global(index, 0);
      yarg_drop(1);
    }
    MPI_Comm_create_errhandler(&ygyotoMPIErrorHandlerFcn,
			       &ygyotoMPIErrorHandler);
    MPI_Comm_set_errhandler(MPI_COMM_WORLD, ygyotoMPIErrorHandler);

#else
    ypush_long(1);
#endif
  }

  void
  Y_gyoto_MPI_Initialized(int argc)
  {
#if defined HAVE_MPI
    int flag=0;
    MPI_Initialized(&flag);
    ypush_long(flag);
#else
    ypush_long(0);
#endif
  }

  void
  Y_gyoto_MPI_Finalize(int)
  {
#if defined HAVE_MPI
    ypush_long(MPI_Finalize());
#else
    ypush_long(1);
#endif
  }

  void
  Y_gyoto_MPI_Finalized(int argc)
  {
#if defined HAVE_MPI
    int flag=0;
    MPI_Finalized(&flag);
    ypush_long(flag);
#else
    ypush_long(0);
#endif
  }

  void
  Y___gyoto_exportSupplier(int)
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

      // Property ABI
      YGyotoGlobalSupplier -> yget_property = &yget_property;
      YGyotoGlobalSupplier -> ypush_property = &ypush_property;
    }
    ypush_long(long(YGyotoGlobalSupplier));
  }

  void
  Y___gyoto_setErrorHandler(int)
  { Gyoto::Error::setHandler(&ygyotoErrorHandler); }

}


/* Don't overuse the two below, the are not exported with the rest of the ABI */
/* The point is to cache the variable names and global indices used by
   the closure on_eval operator */
char const * const __ygyoto_var_name(size_t id) {
  static std::vector<std::string> names;
  if (id >= names.size()) {
    size_t cursize=names.size();
    names.resize(id+1);
    for (size_t k=cursize; k<=id; ++k) {
      stringstream ss;
      ss << "__gyoto_var" << k;
      names[k]=ss.str();
    }
  } 
  return names[id].c_str();
}

long int __ygyoto_var_idx(size_t id) {
  static std::vector<long> ids;
  if (id >= ids.size()) {
    size_t cursize=ids.size();
    ids.resize(id+1);
    for (size_t k=cursize; k<=id; ++k)
      ids[k]=yget_global(__ygyoto_var_name(k), 0);
  } 
  return ids[id];
}

void ypush_property(Gyoto::SmartPointer<Gyoto::SmartPointee> ptr,
		    Gyoto::Property const& p,
		    std::string name, std::string unit) {
  Gyoto::Value val;

  Gyoto::SmartPointee * smptee = (Gyoto::SmartPointee*) ptr();
  Gyoto::Object * object = dynamic_cast<Gyoto::Object*> (smptee);
  Gyoto::Astrobj::Generic * ao=NULL;

  if (!smptee) GYOTO_ERROR("NULL SmartPointee*");

  // Some Astrobj (in particular Star) inherit twice from Object.
  if (!object && (ao=dynamic_cast<Gyoto::Astrobj::Generic*> (smptee)) )
    object = dynamic_cast<Gyoto::Object*> (ao);

  if (!object)
    GYOTO_ERROR("dynamic_cast from SmartPointee* to Object* failed");

  if (p.type == Gyoto::Property::double_t ||
      p.type == Gyoto::Property::vector_double_t)
    val = object -> get(p, unit);
  else
    val = object -> get(p);

  switch(p.type) {
  case Gyoto::Property::bool_t:
    ypush_long(name==p.name?bool(val):!val);
    break;
  case Gyoto::Property::long_t:
    ypush_long(long(val));
    break;
  case Gyoto::Property::unsigned_long_t:
    ypush_long(long((unsigned long)(val)));
    break;
  case Gyoto::Property::size_t_t:
    ypush_long(long(size_t(val)));
    break;
  case Gyoto::Property::double_t:
    ypush_double(val);
    break;
  case Gyoto::Property::string_t:
  case Gyoto::Property::filename_t:
    *ypush_q(0) = p_strcpy(string(val).c_str());
    break;
  case Gyoto::Property::vector_double_t:
    {
      std::vector<double> vval = val;
      size_t n = vval.size();
      long dims[]={1, long(n)};
      double * buf = ypush_d(dims);
      for (size_t i=0; i<n; ++i) buf[i]=vval[i];
    }
    break;
  case Gyoto::Property::vector_unsigned_long_t:
    {
      std::vector<unsigned long> vval = val;
      size_t n = vval.size();
      long dims[]={1, long(n)};
      long * buf = ypush_l(dims);
      for (size_t i=0; i<n; ++i) buf[i]=vval[i];
    }
    break;
  case Gyoto::Property::metric_t:
    *ypush_Metric() = Gyoto::SmartPointer<Gyoto::Metric::Generic>(val);
    break;
  case Gyoto::Property::astrobj_t:
    *ypush_Astrobj() = Gyoto::SmartPointer<Gyoto::Astrobj::Generic>(val);
    break;
  case Gyoto::Property::spectrum_t:
    *ypush_Spectrum() = Gyoto::SmartPointer<Gyoto::Spectrum::Generic>(val);
    break;
  case Gyoto::Property::spectrometer_t:
    *ypush_Spectrometer() = Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>(val);
    break;
  case Gyoto::Property::screen_t:
    *ypush_Screen() = Gyoto::SmartPointer<Gyoto::Screen>(val);
    break;
  default:
    y_error("Property type unimplemented in ypush_property()");
   }
}

void yget_property(Gyoto::SmartPointer<Gyoto::SmartPointee> ptr,
		   Gyoto::Property const& p, int iarg, std::string name,
		   std::string unit) {
  Gyoto::Value val;

  Gyoto::SmartPointee * smptee = (Gyoto::SmartPointee*) ptr();
  Gyoto::Object * object = dynamic_cast<Gyoto::Object*> (smptee);
  Gyoto::Astrobj::Generic * ao=NULL;

  if (!smptee) GYOTO_ERROR("NULL SmartPointee*");

  // Some Astrobj (in particular Star) inherit twice from Object.
  if (!object && (ao=dynamic_cast<Gyoto::Astrobj::Generic*> (smptee)) )
    object = dynamic_cast<Gyoto::Object*> (ao);

  if (!object)
    GYOTO_ERROR("dynamic_cast from SmartPointee* to Object* failed");

  switch(p.type) {
  case Gyoto::Property::bool_t:
    {
      val=bool(ygets_l(iarg));
      if (name != p.name) val = !val;
    }
    break;
  case Gyoto::Property::long_t:
    val = long(ygets_l(iarg));
    break;
  case Gyoto::Property::unsigned_long_t:
    val = (unsigned long)(ygets_l(iarg));
    break;
  case Gyoto::Property::size_t_t:
    val = size_t(ygets_l(iarg));
    break;
  case Gyoto::Property::double_t:
    object->set(p, ygets_d(iarg), unit);
    return;
  case Gyoto::Property::filename_t:
  case Gyoto::Property::string_t:
    val =  string(ygets_q(iarg));
    break;
  case Gyoto::Property::vector_double_t:
    {
      long n;
      double *buf = ygeta_d(iarg, &n, NULL);
      std::vector<double> vval(n, 0.);
      for (long i=0; i<n; ++i) vval[i]=buf[i];
      object->set(p, vval, unit);
    }
    return;
  case Gyoto::Property::vector_unsigned_long_t:
    {
      long n;
      long *buf = ygeta_l(iarg, &n, NULL);
      std::vector<unsigned long> vval=val;
      vval.resize(n);
      for (long i=0; i<n; ++i) vval[i]=buf[i];
    }
    break;
  case Gyoto::Property::metric_t:
    if (yarg_number(iarg) && ygets_l(iarg)==0)
      val=(Gyoto::SmartPointer<Gyoto::Metric::Generic>(NULL));
    else
      val = *yget_Metric(iarg);
    break;
  case Gyoto::Property::astrobj_t:
    if (yarg_number(iarg) && ygets_l(iarg)==0)
      val=(Gyoto::SmartPointer<Gyoto::Astrobj::Generic>(NULL));
    else
      val = *yget_Astrobj(iarg);
    break;
  case Gyoto::Property::spectrum_t:
    if (yarg_number(iarg) && ygets_l(iarg)==0)
      val=(Gyoto::SmartPointer<Gyoto::Spectrum::Generic>(NULL));
    else
      val = *yget_Spectrum(iarg);
    break;
  case Gyoto::Property::spectrometer_t:
    if (yarg_number(iarg) && ygets_l(iarg)==0)
      val=(Gyoto::SmartPointer<Gyoto::Spectrometer::Generic>(NULL));
    else
      val = *yget_Spectrometer(iarg);
    break;
  case Gyoto::Property::screen_t:
    if (yarg_number(iarg) && ygets_l(iarg)==0)
      val=(Gyoto::SmartPointer<Gyoto::Screen>(NULL));
    else
      val = *yget_Screen(iarg);
    break;
  default:
    y_error("Property type unimplemented in yget_property()");
   }
  object->set(p, val);
}
