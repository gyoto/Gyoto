/*
    Copyright 2011 Thibaut Paumard, Frederic Vincent

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

#include "GyotoFactory.h"
#include "GyotoUtils.h"
#include "GyotoRegister.h"
#include "GyotoDefs.h"
#include <fitsio.h>
#include <csignal>

using namespace std;
using namespace Gyoto;

static char*     pixfile   = NULL;
static fitsfile* fptr      = NULL;
static int       status    = 0;
static long      fpixel[]  = {1,1,1};
static long      nelements = 0;
static double*   vect      = NULL;
static double*   impactcoords=NULL;
static SmartPointer<Astrobj::Properties> data = NULL;

void usage() {
  cout << "Usage:" << endl <<
    "    rayXML [--imin=i0 --imax=i1 --jmin=j0 --jmax=j1] input.xml output.dat" << endl;
}

void sigint_handler(int sig)
{
  if (sig!=SIGINT) cerr << "\n********GYOTO: sigint_handler trapping signal " << sig << ", this should not happen !" << endl;
  cerr << "GYOTO: SIGINT received: saving data to " << pixfile << "... ";
  signal(SIGINT, SIG_DFL);

  fits_write_pix(fptr, TDOUBLE, fpixel, nelements, vect, &status);
  fits_close_file(fptr, &status);
  fits_report_error(stderr, status);

  cerr << "Killing self." << endl;
  kill(getpid(), SIGINT);
}

static std::string curmsg = "";
static int curretval = 1;

void gyotoErrorHandler( const char *msg ) {
  cerr << curmsg << msg << endl;
  if (debug()) abort(); // to keep stack for debugger
  exit (curretval);
}

int main(int argc, char** argv) {

  /*
    This program aims at computing the null geodesics of photons from
    an observation screen to an astrophysical object (star orbit,
    fixed star or disk).  The final result is a list of illuminated
    pixels.
  */
  //	For debug output
  debug(0);
  // verbose(1);
  
  char * parfile=NULL;
  string ipctfile="";
  string param;

  size_t imin=1, imax=1000000000, jmin=1, jmax=1000000000;
  //  double tobs, tmin, fov, dist, paln, incl, arg;
  double tobs, fov, dist, paln, incl, arg;
  size_t res, nthreads;
  //  bool  xtobs=0, xtmin=0, xfov=0, xres=0, xdist=0, xpaln=0, xincl=0, xarg=0;
  bool  xtobs=0, xfov=0, xres=0, xdist=0, xpaln=0, xincl=0, xarg=0, xnthreads;
  bool  ipct=0;
  long  ipctdims[3]={0, 0, 0};
  double ipcttime;

  int i, stop=0;
  for (i=1;i<argc;++i) {
    param=argv[i];
    if (param.substr(0,1)=="-" && !stop) {
      if (param=="--") stop=1;
      else if (param.substr(0,10)=="--verbose=") verbose(atoi(param.substr(10).c_str()));
      else if (param.substr(0,8)=="--silent") verbose(0);
      else if (param.substr(0,7)=="--quiet") verbose(GYOTO_QUIET_VERBOSITY);
      else if (param.substr(0,9)=="--verbose") verbose(10);
      else if (param.substr(0,7)=="--debug") debug(1);
      else if (param.substr(0,7)=="--imin=") imin=atoi(param.substr(7).c_str());
      else if (param.substr(0,7)=="--imax=") imax=atoi(param.substr(7).c_str());
      else if (param.substr(0,7)=="--jmin=") jmin=atoi(param.substr(7).c_str());
      else if (param.substr(0,7)=="--jmax=") jmax=atoi(param.substr(7).c_str());
      else if (param.substr(0,15)=="--impact-coords")  {
	if (param.size() > 16 && param.substr(15,1)=="=")
	  ipctfile=param.substr(16);
	else ipct=1;
      }
      else if (param.substr(0,7)=="--time=") {
	tobs=atof(param.substr(7).c_str());
	xtobs=1;
	/*} else if (param.substr(0,7)=="--tmin=") {
	tmin=atof(param.substr(7).c_str());
	xtmin=1;*/
      } else if (param.substr(0,6)=="--fov=") {
	fov=atof(param.substr(6).c_str());
	xfov=1;
      } else if (param.substr(0,13)=="--resolution=") {
	res=atoi(param.substr(13).c_str());
	xres=1;
      } else if (param.substr(0,11)=="--distance=") {
	dist=atof(param.substr(11).c_str());
	xdist=1;
      } else if (param.substr(0,7)=="--paln=") {
	paln=atof(param.substr(7).c_str());
	xpaln=1;
      } else if (param.substr(0,14)=="--inclination=") {
	incl=atof(param.substr(14).c_str());
	xincl=1;
      } else if (param.substr(0,11)=="--argument=") {
	arg=atof(param.substr(11).c_str());
	xarg=1;
      }  else if (param.substr(0,11)=="--nthreads=") {
	nthreads=atoi(param.substr(11).c_str());
	xnthreads=1;
      }
      else {
	usage();
	return 1;
      }
    } else {
      if (!parfile) parfile=argv[i];
      else if (!pixfile) pixfile=argv[i];
      else {
	usage();
	return 1;
      }
    }
  }

  if (!pixfile) {
    usage();
    return 1;
  }

  // State copyright
  if (verbose() >= GYOTO_QUIET_VERBOSITY)
    cout << " Copyright (c) 2011 Frederic Vincent & Thibaut Paumard\n"
	 << " GYOTO is distributed under the terms of the GPL v. 3 license.\n"
	 << " We request that use of Gyoto in scientific publications be "
	 << " properly \n acknowledged. Please cite:\n"
	 << "  GYOTO: a new general relativistic ray-tracing code,\n"
	 << "  F. H. Vincent, T. Paumard, E. Gourgoulhon & G. Perrin 2011,\n"
	 << "  Classical and Quantum Gravity, accepted. [arXiv:1109.4769]"
	 << endl << endl;

  // set-up error reporter
  Gyoto::setErrorHandler ( &gyotoErrorHandler );

  curmsg = "In gyoto.C: Error initializing libgyoto: ";
  curretval = 1;
  Gyoto::Register::init();

  Factory *factory ;
  if (verbose() >= GYOTO_QUIET_VERBOSITY) cout << "Reading parameter file: " << parfile << endl;
  curmsg = "In gyoto.C: Error in Factory creation: ";
  curretval = 1;
  factory = new Factory(parfile);

  curmsg = "In gyoto.C: Error getting Kind: ";
  const string kind = factory->getKind();

  if (!kind.compare("Scenery")) {
    curmsg = "In gyoto.C: Error initializing ray-tracing: ";
    curretval = 2;
    SmartPointer<Scenery> scenery = factory -> getScenery();
    SmartPointer<Screen>  screen = scenery->getScreen();
    SmartPointer<Astrobj::Generic> object = scenery->getAstrobj();

    if (xtobs) screen -> setTime        ( tobs );
    else tobs= screen -> getTime();
    //      if (xtmin) screen -> setMinimumTime ( tmin );
    if (xres)  screen -> setResolution  ( res  );
    else res = screen -> getResolution();
    if (xfov)  screen -> setFieldOfView ( fov  );
    if (xdist) screen -> setDistance    ( dist );
    if (xincl) screen -> setInclination ( incl );
    if (xpaln) screen -> setPALN        ( paln );
    if (xarg)  screen -> setArgument    ( arg  );
    if (xnthreads)  scenery -> setNThreads    ( nthreads  );

    if (ipctfile != "") {
      //	  if (verbose() >= GYOTO_QUIET_VERBOSITY)
      size_t ipctnelt=0;
      cout << "Reading precomputed impact coordinates from " << ipctfile <<endl;
      fits_open_file(&fptr, ipctfile.c_str(), 0, &status);
      fits_movnam_hdu(fptr, ANY_HDU,
		      const_cast<char*>("Gyoto Impact Coordinates"),
		      0, &status);
      fits_read_key(fptr, TDOUBLE, "Gyoto Observing Date", &ipcttime,
		    NULL, &status);
      fits_get_img_size(fptr, 3, ipctdims, &status);
      fits_report_error(stderr, status);
      if (status) return status;

      if (ipctdims[0]==16 &&
	  size_t(ipctdims[1]) == res &&
	  size_t(ipctdims[2]) == res) {
	impactcoords = new double[(ipctnelt=16*res*res)];
      } else {
	cerr<<"ERROR: bad dimensions for precomputed impact coordinates\n";
	return 1;
      }
	  
      fits_read_subset(fptr, TDOUBLE, fpixel, ipctdims, fpixel,
		       0, impactcoords, NULL, &status);
      fits_close_file(fptr, &status);
      fptr=NULL;
      fits_report_error(stderr, status);
      if (status) return status;

      double dt = tobs * GYOTO_C / scenery -> getMetric() -> unitLength()
	- ipcttime;
      for (i=0; i < ipctnelt; i+=8)
	if (impactcoords[i] != DBL_MAX) impactcoords[i] += dt;
      ipcttime = tobs * GYOTO_C / scenery -> getMetric() -> unitLength();
    }

    Quantity_t quantities = scenery -> getRequestedQuantities();
    if (debug()) cerr << "DEBUG: Gyoto.C: Requested Quantities: "
		      << quantities <<endl;

    size_t nbnuobs=0;
    if (quantities & (GYOTO_QUANTITY_SPECTRUM | GYOTO_QUANTITY_BINSPECTRUM)) {
      SmartPointer<Spectrometer> spr = screen -> getSpectrometer();
      if (!spr) throwError("Spectral quantity requested but "
			   "no spectrometer specified!");
      nbnuobs = spr -> getNSamples();
    }
               //nb of frames that will be used for spectral cube
    size_t nbdata= scenery->getScalarQuantitiesCount();
               //nb of frames used for diverse interesting outputs
               //(obs flux, impact time, redshift..)
    size_t nelt=res*res*(nbdata+nbnuobs);
    vect = new double[nelt];

    // First check whether we can open file
    long naxis=3; 
    long naxes[] = {res, res, nbdata+nbnuobs};
    nelements=nelt; 

    fits_create_file(&fptr, pixfile, &status);
    fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
    fits_report_error(stderr, status);
    if (status) return status;

    // Allocate space for the output data
    data = new Astrobj::Properties();

    size_t curquant=0;
    size_t offset=res*res;

    if (debug()) {
      cerr << "DEBUG: gyoto.C: flag_radtransf = ";
      cerr << scenery -> getAstrobj() -> getFlag_radtransf() << endl;
      cerr << "DEBUG: gyoto.C: Requested quantities: ";
      cerr << scenery -> getRequestedQuantitiesString() << endl;
    }

    char keyname[FLEN_KEYWORD];
    char const * fmt="QUANT_%lu";
    char * CNULL=NULL;

    if (quantities & GYOTO_QUANTITY_INTENSITY) {
      data->intensity=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("Intensity"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_EMISSIONTIME) {
      data->time=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("EmissionTime"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_MIN_DISTANCE) {
      data->distance=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("MinDistance"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_FIRST_DMIN) {
      data->first_dmin=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("FirstDistMin"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_REDSHIFT) {
      if (debug())
	cerr << "DEBUG: gyoto.C: REDSHIFT requested\n";
      data->redshift=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("Redshift"),
		     CNULL, &status);
    }
    if ((quantities & GYOTO_QUANTITY_IMPACTCOORDS || ipct) && !ipctdims[0] ) {
      // Allocate if requested AND not provided
      cerr << "gyoto.C: allocating data->impactcoords" << endl;
      data->impactcoords = impactcoords = new double [res*res*16];
      ipcttime = tobs * GYOTO_C / scenery -> getMetric() -> unitLength();
    }
    if (quantities & GYOTO_QUANTITY_USER1) {
      data->user1=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("User1"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_USER2) {
      data->user2=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("User2"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_USER3) {
      data->user3=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("User3"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_USER4) {
      data->user4=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("User4"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_USER5) {
      data->user5=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("User5"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_SPECTRUM) {
      data->spectrum=vect+offset*(curquant++);
      data->offset=offset;
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("Spectrum"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_BINSPECTRUM) {
      data->binspectrum=vect+offset*(curquant++);
      data->offset=offset;
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("BinSpectrum"),
		     CNULL, &status);
    }
    
    signal(SIGINT, sigint_handler);

    curmsg = "In gyoto.C: Error during ray-tracing: ";
    scenery -> rayTrace(imin, imax, jmin, jmax, data,
			ipctdims[0]?impactcoords:NULL);

    curmsg = "In gyoto.C: Error while saving: ";
    if (verbose() >= GYOTO_QUIET_VERBOSITY)
      cout << "\nSaving to file: " << pixfile << endl;
    signal(SIGINT, SIG_DFL);


    // Save to fits file
    fits_write_pix(fptr, TDOUBLE, fpixel, nelements, vect, &status);

    if (quantities & GYOTO_QUANTITY_IMPACTCOORDS || ipct) {
      // Save if requested, copying if provided
      cout << "Saving precomputed impact coordinates" << endl;
      long naxes_ipct[] = {16, res, res};
      fits_create_img(fptr, DOUBLE_IMG, naxis, naxes_ipct, &status);
      fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
		     const_cast<char*>("Gyoto Impact Coordinates"),
		     CNULL, &status);
      fits_write_key(fptr, TDOUBLE, const_cast<char*>("Gyoto Observing Date"),
		     &ipcttime, "Geometrical units", &status);

      fits_write_pix(fptr, TDOUBLE, fpixel, res*res*16, impactcoords, &status);

      fits_report_error(stderr, status);
      if (status) return status;
    }

    fits_close_file(fptr, &status);
    fits_report_error(stderr, status);
    if (debug()) cerr << "DEBUG: gyoto.C: FITS file closed, cleaning" << endl;

    curmsg = "In gyoto.C: Error while cleaning (file saved already): ";

    if (debug()) cerr << "DEBUG: gyoto.C: delete [] vect" << endl;
    delete [] vect;

    if (impactcoords) {
      if (debug()) cerr << "gyoto.C: delete [] data->impact" << endl;
      delete [] impactcoords;
    }

    if (debug()) cerr << "DEBUG: gyoto.C: scenery==NULL" << endl;
    scenery = NULL;


    if (status) return status;

  } else {
    cerr << "Unknown kind for root element in XML file" << endl;
    return 1;
  }

  delete factory;
  return 0;
}
