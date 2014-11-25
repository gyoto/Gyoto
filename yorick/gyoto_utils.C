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
#include <pstdlib.h>
#include <cstring>
#include <iostream>
#include <sstream>
#include <signal.h>
#include <vector>

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


  void
  Y_gyoto_dontcatchSIGSEGV(int)
  {
    signal(SIGSEGV, SIG_DFL);
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
  Y___gyoto_initRegister(int)
  {
    Gyoto::Register::init();
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
char const * const __ygyoto_var_name(long id) {
  static std::vector<std::string> names;
  if (id >= names.size()) {
    long cursize=names.size();
    names.resize(id+1);
    for (long k=cursize; k<=id; ++k) {
      stringstream ss;
      ss << "__gyoto_var" << k;
      names[k]=ss.str();
    }
  } 
  return names[id].c_str();
}

long int __ygyoto_var_idx(long id) {
  static std::vector<long> ids;
  if (id >= ids.size()) {
    long cursize=ids.size();
    ids.resize(id+1);
    for (long k=cursize; k<=id; ++k)
      ids[k]=yget_global(__ygyoto_var_name(k), 0);
  } 
  return ids[id];
}

void ypush_property(Gyoto::SmartPointer<Gyoto::SmartPointee> ptr,
		    Gyoto::Property const& p, int iarg,
		    std::string name, std::string unit) {
  switch(p.type) {
  case Gyoto::Property::bool_t:
    {
      bool val;
      ptr->get(p, val);
      ypush_long(name==p.name?val:!val);
    }
    break;
  case Gyoto::Property::double_t:
    {
      double val;
      ptr->get(p, val, unit);
      ypush_double(val);
    }
    break;
  default:
    y_error("Property type unimplemented in ypush_property()");
   }
}

void yget_property(Gyoto::SmartPointer<Gyoto::SmartPointee> ptr,
		   Gyoto::Property const& p, int iarg, std::string name,
		   std::string unit) {
  cerr << "dummy: should get property" << endl;
  switch(p.type) {
  case Gyoto::Property::bool_t:
    {
      bool val=ygets_l(iarg);
      if (name != p.name) val = !val;
      ptr->set(p, val);
    }
    break;
  case Gyoto::Property::double_t:
    ptr->set(p, ygets_d(iarg), unit);
    break;
  default:
    y_error("Property type unimplemented in yget_property()");
   }
}
