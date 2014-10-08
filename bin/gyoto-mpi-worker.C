/*
    Copyright 2011, 2013 Thibaut Paumard, Frederic Vincent

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

// Gyoto
#include "GyotoDefs.h"
#include "GyotoFactory.h"
#include "GyotoUtils.h"
#include "GyotoRegister.h"
#include "GyotoScenery.h"

// FITS I/O
#include <fitsio.h>

// signal()
#include <csignal>

// getpid()
#include <sys/types.h>
#include <unistd.h>
 
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

static SmartPointer<Scenery> sc = NULL;

#ifdef HAVE_MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/communicator.hpp>
#include <string>
#include <boost/serialization/string.hpp>
namespace mpi = boost::mpi;
#endif

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

void gyotoErrorHandler( const Gyoto::Error e ) {
  cerr << curmsg << e << endl;
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

  mpi::environment env;


  MPI_Comm parent_c;
  MPI_Comm_get_parent(&parent_c);

  mpi::communicator manager(parent_c,mpi::comm_take_ownership);
  mpi::communicator world;

  string pluglist= getenv("GYOTO_PLUGINS")?
    getenv("GYOTO_PLUGINS"):
    GYOTO_DEFAULT_PLUGINS;
  Gyoto::Error::setHandler ( &gyotoErrorHandler );
  curmsg = "In gyoto.C: Error initializing libgyoto: ";
  curretval = 1;
  Gyoto::Register::init(pluglist.c_str());

  int rk=world.rank();

  //Gyoto::debug(1);

  Scenery::mpi_tag task=Scenery::give_task;
  Scenery::is_worker=true;
  while (task != Scenery::terminate) {
    //    std::cerr << "Worker #" << rk << " waiting for task " << endl;
    manager.recv(0, Scenery::give_task, task);
    //    std::cerr << "Worker #" << rk << " received task " << task << endl;
    switch (task) {
    case Scenery::read_scenery: {
      std::string parfile;
      manager.recv(0, task, parfile);
      std::cerr << "Worker #" << rk << " reading \""<<parfile<<"\""<<std::endl;
      curmsg = "In gyoto.C: Error in Factory creation: ";
      curretval = 1;
      sc = Factory(const_cast<char*>(parfile.c_str())).getScenery();
      cerr <<"Maxiter=="<<sc->maxiter()<<endl;
      break;
    }
    case Scenery::raytrace: {
      size_t ij[2]={0, 0};
      //      std::cerr << "Worker #" << rk << " receiving ij for raytracing" <<std::endl;
      manager.recv(0, task, ij, 2);
      std::cerr << "Worker #" << rk << " raytracing i="<<ij[0]<<", j="<<ij[1]<<std::endl;


      // initialize AstrobjProperties
      size_t nbnuobs=0;
      Quantity_t quantities = sc -> getRequestedQuantities();
      if (quantities & (GYOTO_QUANTITY_SPECTRUM | GYOTO_QUANTITY_BINSPECTRUM)) {
	SmartPointer<Spectrometer::Generic> spr = sc -> screen() -> spectrometer();
	if (!spr) throwError("Spectral quantity requested but "
			     "no spectrometer specified!");
	nbnuobs = spr -> nSamples();
      }
      size_t nbdata= sc->getScalarQuantitiesCount();
      size_t nelt=(nbdata+nbnuobs);
      vect = new double[nelt];
      data = new Astrobj::Properties();
      size_t offset=1;
      size_t curquant=0;

      if (quantities & GYOTO_QUANTITY_INTENSITY) {
	data->intensity=vect+offset*(curquant++);
      }
      if (quantities & GYOTO_QUANTITY_EMISSIONTIME) {
	data->time=vect+offset*(curquant++);
      }
      if (quantities & GYOTO_QUANTITY_MIN_DISTANCE) {
	data->distance=vect+offset*(curquant++);
      }
      if (quantities & GYOTO_QUANTITY_FIRST_DMIN) {
	data->first_dmin=vect+offset*(curquant++);
      }
      if (quantities & GYOTO_QUANTITY_REDSHIFT) {
	data->redshift=vect+offset*(curquant++);
      }
      // if ((quantities & GYOTO_QUANTITY_IMPACTCOORDS || ipct) && !ipctdims[0] ) {
      // 	data->impactcoords = impactcoords = new double [res*res*16];
      // 	ipcttime = tobs * GYOTO_C / scenery -> metric() -> unitLength();
      // }
      if (quantities & GYOTO_QUANTITY_USER1) {
	data->user1=vect+offset*(curquant++);
      }
      if (quantities & GYOTO_QUANTITY_USER2) {
	data->user2=vect+offset*(curquant++);
      }
      if (quantities & GYOTO_QUANTITY_USER3) {
	data->user3=vect+offset*(curquant++);
      }
      if (quantities & GYOTO_QUANTITY_USER4) {
	data->user4=vect+offset*(curquant++);
      }
      if (quantities & GYOTO_QUANTITY_USER5) {
	data->user5=vect+offset*(curquant++);
      }
      if (quantities & GYOTO_QUANTITY_SPECTRUM) {
	data->spectrum=vect+offset*(curquant++);
	data->offset=int(offset);
      }
      if (quantities & GYOTO_QUANTITY_BINSPECTRUM) {
	data->binspectrum=vect+offset*(curquant++);
	data->offset=int(offset);
      }

      data->init(nbnuobs);
      sc->screen()->computeBaseVectors();
      sc->setPropertyConverters(data);

      // need to do something about impactcoords

      (*sc)(ij[0], ij[1], data, NULL, NULL);

      manager.send(0, Scenery::raytrace_done, rk);
      manager.send(0, Scenery::raytrace_done, ij, 2);
      manager.send(0, Scenery::raytrace_done, nelt);
      manager.send(0, Scenery::raytrace_done, vect, nelt);

      delete[] vect;
      delete data;
      //      std::cerr << "Worker #" << rk << " done raytracing i="<<ij[0]<<", j="<<ij[1]<<std::endl;

    }
      break;
    case Scenery::terminate:
      std::cerr << "Worker #" << rk << " terminating"<<std::endl;

      break;
    default:
      std::cerr << "unknown task" << endl;
    }
  }

  std::cerr << "Worker #" << rk << " done, exiting" << endl;

  return 0;
}
