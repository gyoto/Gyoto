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
#include "GyotoScreen.h"

// FITS I/O
#include <fitsio.h>

// signal()
#include <csignal>

// getpid()
#include <sys/types.h>
#include <unistd.h>
 
using namespace std;
using namespace Gyoto;

static SmartPointer<Scenery> sc = NULL;

#ifdef HAVE_MPI
#include <boost/mpi/environment.hpp>
#include <boost/mpi/intercommunicator.hpp>
#include <boost/mpi/collectives.hpp>
#include <string>
#include <boost/serialization/string.hpp>
namespace mpi = boost::mpi;
#endif


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
  verbose(GYOTO_SEVERE_VERBOSITY);

  mpi::environment env(argc, argv);


  MPI_Comm parent_c;
  MPI_Comm_get_parent(&parent_c);

  mpi::communicator team=mpi::intercommunicator(parent_c,mpi::comm_take_ownership).merge(true);

  string pluglist= getenv("GYOTO_PLUGINS")?
    getenv("GYOTO_PLUGINS"):
    GYOTO_DEFAULT_PLUGINS;
  Gyoto::Error::setHandler ( &gyotoErrorHandler );
  curmsg = "In gyoto-mpi-worker.C: Error initializing libgyoto: ";
  curretval = 1;
  Gyoto::Register::init(pluglist.c_str());

  sc = new Scenery();
  sc -> mpi_team_    = &team;

  Screen::Empty empty;

  Scenery::mpi_tag task=Scenery::give_task;
  Scenery::am_worker=true;
  while (task != Scenery::terminate) {
    sc->mpiTask(task);
    switch (task) {
    case Scenery::read_scenery: {
      std::string parfile;
      broadcast(team, parfile, 0);
      curmsg = "In gyoto-mpi-worker.C: Error in Factory creation: ";
      curretval = 1;
      sc = Factory(const_cast<char*>(parfile.c_str())).getScenery();
      sc -> mpi_team_    = &team;
     break;
    }
    case Scenery::raytrace:
      curmsg = "In gyoto-mpi-worker.C: Error in ray-tracing: ";
      curretval = 1;
      sc -> rayTrace(empty, NULL, NULL);
      break;
    case Scenery::terminate:
      sc = NULL;
      break;
    default:
      std::cerr << "unknown task " << task << endl;
    }
  }

  team.barrier();
  return 0;
}
