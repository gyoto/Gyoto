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
static SmartPointer<AstrobjProperties> data = NULL;

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
  string param;

  size_t imin=1, imax=1000000000, jmin=1, jmax=1000000000;
  int save=0;
  //  double tobs, tmin, fov, dist, paln, incl, arg;
  double tobs, fov, dist, paln, incl, arg;
  size_t res;
  //  bool  xtobs=0, xtmin=0, xfov=0, xres=0, xdist=0, xpaln=0, xincl=0, xarg=0;
  bool  xtobs=0, xfov=0, xres=0, xdist=0, xpaln=0, xincl=0, xarg=0;


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
      else if (param.substr(0,6)=="--save")  save=1;
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
    SmartPointer<Astrobj> object = scenery->getAstrobj();

    if (xtobs) screen -> setTime        ( tobs );
    //      if (xtmin) screen -> setMinimumTime ( tmin );
    if (xres)  screen -> setResolution  ( res  );
    if (xfov)  screen -> setFieldOfView ( fov  );
    if (xdist) screen -> setDistance    ( dist );
    if (xincl) screen -> setInclination ( incl );
    if (xpaln) screen -> setPALN        ( paln );
    if (xarg)  screen -> setArgument    ( arg  );

    res = screen -> getResolution();

    if (debug()) cout << "DEBUG: gyoto.C: Nb of pixels= " << res << endl;

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
    data = new AstrobjProperties();

    size_t curquant=0;
    size_t offset=res*res;

    if (debug()) {
      cerr << "DEBUG: gyoto.C: flag_radtransf = ";
      cerr << scenery -> getAstrobj() -> getFlag_radtransf() << endl;
      cerr << "DEBUG: gyoto.C: Requested quantities: ";
      cerr << scenery -> getRequestedQuantitiesString() << endl;
    }

    char keyname[FLEN_KEYWORD];
    char * fmt="QUANT_%lu";
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
    if (quantities & GYOTO_QUANTITY_IMPACT_R) {
      data->rimpact=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("ImpactR"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_IMPACT_X) {
      data->x=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("ImpactX"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_IMPACT_Y) {
      data->y=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("ImpactY"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_IMPACT_Z) {
      data->z=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("ImpactZ"),
		     CNULL, &status);
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

    // scenery -> rayTrace(imin, imax, jmin, jmax, &data, save);
    curmsg = "In gyoto.C: Error during ray-tracing: ";
    scenery -> rayTrace(imin, imax, jmin, jmax, data, save);
    scenery = NULL;

    curmsg = "In gyoto.C: Error while saving: ";
    if (verbose() >= GYOTO_QUIET_VERBOSITY) cout << "\nSaving to file: " << pixfile << endl;
    signal(SIGINT, SIG_DFL);


    // Save to fits file
    fits_write_pix(fptr, TDOUBLE, fpixel, nelements, vect, &status);
    fits_close_file(fptr, &status);
    fits_report_error(stderr, status);

    curmsg = "In gyoto.C: Error while cleaning (file saved already): ";
    delete [] vect;

    if (status) return status;

  } else {
    cerr << "Unknown kind for root element in XML file" << endl;
    return 1;
  }

  delete factory;
  return 0;
}
