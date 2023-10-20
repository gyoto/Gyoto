/*
    Copyright 2011, 2013-2015, 2018 Thibaut Paumard & Frédéric Vincent

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

#include "GyotoPhoton.h"
#include "GyotoFactory.h"
#include "GyotoDefs.h"
#include "ygyoto.h"
#include "ygyoto_private.h"
#include "yapi.h"
#include "pstdlib.h"

#define OBJ ph

using namespace Gyoto;
using namespace std;

#define YGYOTO_PHOTON_GENERIC_KW "unit", "metric", "initcoord", "astrobj", \
    "spectro", "tmin",	"delta", "adaptive", "maxiter", \
    "integrator", "deltamin", "deltamax", "deltamaxoverr", "abstol", "reltol", \
    "setparameter", "reset", "xfill", "save_txyz", "xmlwrite", "is_hit", \
    "get_txyz", "get_coord", "get_cartesian", "clone", "help"
#define YGYOTO_PHOTON_GENERIC_KW_N 26

void ygyoto_Photon_generic_eval(Gyoto::SmartPointer<Gyoto::Photon>* ph,
				 int *kiargs, int *piargs, int *rvset, int *paUsed) {
  if (debug()) cerr << "\nDEBUG: In ygyoto_Photon_generic_eval: ";
  int k=-1, iarg;
  char const * rmsg="Cannot set return value more than once";
  char const * pmsg="Cannot use positional argument more than once";
  char * unit = NULL;

  YGYOTO_WORKER_SET_UNIT;
  YGYOTO_WORKER_GETSET_OBJECT2(metric,Metric);

  /* INITCOORD */
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (debug()) cerr << "initcoord=";
    if (yarg_nil(iarg)) { // Getting initcoord
      if (debug()) cerr << "     get_initcoord=1" << endl;
      if ((*rvset)++) y_error(rmsg);
      state_t v(8);
      (*ph)->getInitialCoord(v);
      long dims[]= {1, 8};
      double *coord = ypush_d(dims);
      for (int i=0; i<8; ++i) coord[i]=v[i];
      //TODO: support polarization (v.size()==16)
    } else {          // Setting initcoord
      long ptot=1;
      double coord[8], Ephi[4], Etheta[4];
      bool compute_polar_basis=false;

      if (yarg_number(iarg)) { // initcoord=coord or initcoord=pos, vel
	double *pos=ygeta_d(iarg, &ptot, 0);
	if (ptot!=4 && ptot !=8)
	  y_error("POS should have either 4 or 8 elements.");
	for (int ii=0; ii<ptot; ++ii) coord[ii]=pos[ii];
	if (ptot==4) {
	  if ((*paUsed)++) y_error(pmsg);
	  long vtot=1;
	  double *vel=ygeta_d(piargs[0]+*rvset, &vtot, 0);
	  if (vtot==4) for (int ii=0; ii<4; ++ii) coord[ii+4]=vel[ii];
	  else if (vtot==3) {
	    for (int ii=0; ii<3; ++ii) coord[ii+5]=vel[ii];
	    if (!(*ph)->metric())
	      y_error("METRIC should have been set already");
	    (*ph)->metric()->nullifyCoord(coord);
	  } else y_error("VEL should have 3 or 4 elements");
	}
      } else {
	SmartPointer<Screen> sc = NULL;
	if (yarg_Scenery(iarg)) { // initcoord=scenery,i,j or senery,da,dd
	  SmartPointer<Scenery> *scenery = yget_Scenery(iarg);
	  (*ph) -> metric((*scenery)->metric());
	  (*ph) -> astrobj((*scenery)->astrobj());
	  sc = (*scenery)->screen();
	} else sc = *yget_Screen(iarg); //initcoord=screen,i,j or screen, da, dd
	if (yarg_number(piargs[0]+*rvset)==2) {
	  double da = ygets_d(piargs[0]+*rvset);
	  double dd = ygets_d(piargs[1]+*rvset);
	  if (debug()) cerr << "screen, dalpha="<<da <<", ddelta="<<dd <<endl;
	  sc -> getRayTriad(da, dd, 
			    coord,
			    compute_polar_basis,
			    Ephi, Etheta);
	} else {
	  size_t i = ygets_l(piargs[0]+*rvset);
	  size_t j = ygets_l(piargs[1]+*rvset);
	  if (debug()) cerr << "screen, i="<<i <<", j="<<j << endl;
	  sc -> getRayTriad(i, j, 
			    coord,
			    compute_polar_basis,
			    Ephi, Etheta);
	}
	(*ph) -> spectrometer(sc->spectrometer());
      }

      if (debug()) {
	cerr << "DEBUG:        initcoord=[" ;
	for (int ii=0; ii<7; ++ii) cerr << coord[ii] << ", ";
	cerr << coord[7] ;
      }
      (*ph)->setInitCoord(coord) ;
      if (debug()) cerr << "]" << endl;
    }
    if (debug()) cerr << "... " << endl;
  }

  YGYOTO_WORKER_GETSET_OBJECT2(astrobj,Astrobj);
  YGYOTO_WORKER_GETSET_OBJECT2(spectrometer,Spectrometer);
  YGYOTO_WORKER_GETSET_DOUBLE2(tMin);
  YGYOTO_WORKER_GETSET_DOUBLE2_UNIT(delta);
  YGYOTO_WORKER_GETSET_LONG2(adaptive);
  YGYOTO_WORKER_GETSET_LONG2( maxiter );
  YGYOTO_WORKER_GETSET_STRING2( integrator );
  YGYOTO_WORKER_GETSET_DOUBLE2( deltaMin );
  YGYOTO_WORKER_GETSET_DOUBLE2( deltaMax );
  YGYOTO_WORKER_GETSET_DOUBLE2( deltaMaxOverR );
  YGYOTO_WORKER_GETSET_DOUBLE2( absTol );
  YGYOTO_WORKER_GETSET_DOUBLE2( relTol );
  YGYOTO_WORKER_SETPARAMETER;
  YGYOTO_WORKER_RUN( reset() );
  YGYOTO_WORKER_RUN( xFill(ygets_d(iarg), false) );

  // save_txyz=filename, t1, mass_sun, distance_kpc, unit, screen
  if ((iarg=kiargs[++k])>=0) {
    iarg+=*rvset;
    if (debug()) cerr << "     save_txyz=";
    int myarg=0;
    if (debug()) cerr << "save_txyz=";
    char*  filename          = ygets_q(iarg);
    double t1                = 0.;//ygets_d(piargs[0]+*rvset);
    double mass_sun          = 1.;//ygets_d(piargs[1]+*rvset);
    double distance_kpc      = 1.;//ygets_d(piargs[2]+*rvset);
    string unt               = "geometrical";//ygets_q(piargs[3]+*rvset);
    SmartPointer<Screen>* sc = NULL; //yget_Screen(piargs[4]+*rvset);

    if (yarg_number(piargs[myarg]+*rvset)) {
      t1= ygets_d(piargs[myarg++]+*rvset);
      if (yarg_number(piargs[myarg]+*rvset)) {
	mass_sun= ygets_d(piargs[myarg++]+*rvset);
	if (yarg_number(piargs[myarg]+*rvset)) {
	  distance_kpc= ygets_d(piargs[myarg++]+*rvset);
	}
      }
    }
    if (yarg_string(piargs[myarg]+*rvset))
      unt = ygets_q(piargs[myarg++]+*rvset);
    if (yarg_Screen(piargs[myarg]+*rvset))
      sc = yget_Screen(piargs[myarg++]+*rvset);

    (*ph) -> save_txyz ( filename, t1, mass_sun, distance_kpc, unt,
			 (sc?*sc:SmartPointer<Screen>(NULL)) );

    if (debug()) cerr << filename << endl ;
  }

  YGYOTO_WORKER_XMLWRITE;

  // FUNCTION-LIKE //

  /* IS_HIT= */
  if ((iarg=kiargs[++k])>=0) { // is_hit=whatever
    if (debug()) cerr << "     is_hit=" << endl;
    if ((*rvset)++) y_error(rmsg);
    Astrobj::Properties junk;
    ypush_int((*ph)->hit(&junk));
  }

  /* GET_TXYZ */
  if ((iarg=kiargs[++k])>=0) { // get_txyz=
    if (debug()) cerr << "     get_txyz=1" << endl;
    if ((*rvset)++) y_error(rmsg);
    long nel =(*ph)->get_nelements();
      
    long dims[] = {2, nel, 4};
    double * data=ypush_d(dims);

    (*ph)->get_t(data);
    (*ph)->get_xyz(data+nel, data+2*nel, data+3*nel);
  }

  /* GET_COORD */
  if ((iarg=kiargs[++k])>=0) {
    if (debug()) cerr << "     get_coord=dates" << endl;
    if ((*rvset)++) y_error(rmsg);
    // Two cases : get_coord=1 or get_coord=dates. 
    // We will recognize the first case if the parameter is the
    // integer scalar 1, any other type or value is a date
    // specification.
    int argt = 1;yarg_number(iarg);
    long ntot[1] = {1};
    long dims[Y_DIMSIZE] = {0};
    double * dates = NULL;

    if (!yarg_nil(iarg)) {
      argt = yarg_number(iarg);
      dates=ygeta_d(iarg, ntot, dims);
    } else {
      dates = ypush_d(dims);
      dates[0]=1.;
    }

    if (debug())
      cerr << "DEBUG: gyoto_Photon(get_coord=array("
	   << (yarg_number(iarg)==1?"long":"double")
	   << ", " << *ntot 
	   << ") (rank=" << dims[0]
	   << ", dates[0]=" << dates[0] << ")" << endl;

    if (dims[0] == 0 && argt == 1 && dates[0] == 1) {
      if(debug())
	cerr << "DEBUG: retrieving all already computed coordinates" << endl;
      long nel =(*ph)->get_nelements();
      long ddims[] = {2, nel, 8};
      double * data = ypush_d(ddims);

      (*ph)->getCoord(data, data+nel, data+2*nel, data+3*nel);
      (*ph)->get_dot(data+4*nel, data+5*nel, data+6*nel, data+7*nel);
    } else {
      if(debug())
	cerr << "DEBUG: retrieving coordinates for specified date(s)"<< endl;
      if (dims[0] > Y_DIMSIZE-2)
	y_errorn("gyoto_Star(get_coord=dates): DATES must have at most %ld "
		 "dimensions", Y_DIMSIZE-2);
      dims[0] += 1;
      dims[dims[0]] = 7;
      double * data = ypush_d(dims);
      size_t nel = *ntot;
      (*ph) -> getCoord(dates, nel,
			data, data+ nel,data+2* nel,data+3* nel,
			data+4* nel, data+5* nel, data+6* nel);
    }
  }

  if ((iarg=kiargs[++k])>=0) { // get_cartesian
    if (debug()) cerr << "     get_cartesian=dates" << endl;
    if ((*rvset)++) y_error(rmsg);
    long ntot[1] = {1};
    long dims[Y_DIMSIZE];
    double* dates = ygeta_d(iarg, ntot, dims);
    if (dims[0] > Y_DIMSIZE-2)
      y_errorn("gyoto_Star(get_cartesian=dates): DATES must have at most %ld "
	       "dimensions", Y_DIMSIZE-2);
    dims[0] += 1;
    dims[dims[0]] = 6;
    double * data = ypush_d(dims);
    size_t nel = *ntot;
    (*ph) -> getCartesian(dates, nel,
			  data,         data +   nel, data + 2*nel,
			  data + 3*nel, data + 4*nel, data + 5*nel);
    
  }

  YGYOTO_WORKER_CLONE(Photon);
  YGYOTO_WORKER_HELP;

  if (debug()) cerr << endl;
}

YGYOTO_YUSEROBJ(Photon, Photon)

extern "C" {
  void gyoto_Photon_eval(void *obj, int argc) {
    // If no parameters, return pointer
    if (argc==1 && yarg_nil(0)) {
      ypush_long((long)((gyoto_Photon*)obj)->smptr());
      return;
    }

    SmartPointer<Photon> *ph_ = &(((gyoto_Photon*)obj)->smptr);

    static char const * knames[]={
      YGYOTO_PHOTON_GENERIC_KW, 0
    };

    YGYOTO_WORKER_INIT1(Photon, Photon, knames, YGYOTO_PHOTON_GENERIC_KW_N);
   
    ygyoto_Photon_generic_eval(ph_, kiargs, piargs, rvset, paUsed);

  }

  // Generic constructor/accessor
  void
  Y_gyoto_Photon(int argc)
  {
    YGYOTO_CONSTRUCTOR_INIT2(Photon, Photon, Photon, photon);
    gyoto_Photon_eval(OBJ, argc);
  }

}



