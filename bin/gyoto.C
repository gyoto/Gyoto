/*
    Copyright 2011-2016, 2018-2021 Thibaut Paumard, Frederic Vincent

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
#include "optionparser.h"

// feenableexcept()
#if defined HAVE_FENV_H
# include <fenv.h>
#endif

// FITS I/O
#include <fitsio.h>

// signal()
#include <csignal>

// getpid(), execlp
#include <sys/types.h>
#include <unistd.h>

// MPI
#if defined HAVE_MPI
# include <mpi.h>
# include <boost/mpi/communicator.hpp>
#endif

// ULONG_MAX
#include <climits>

// dlsym
#include <dlfcn.h>

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

namespace Gyoto {
  struct Arg;
}

struct Gyoto::Arg: public option::Arg
{
  static void printError(const char* msg1, const option::Option& opt, const char* msg2)
  {
    fprintf(stderr, "ERROR: %s", msg1);
    fwrite(opt.name, opt.namelen, 1, stderr);
    fprintf(stderr, "%s", msg2);
  }
  static option::ArgStatus Required(const option::Option& option, bool msg)
  {
    if (option.arg != 0)
      return option::ARG_OK;

    if (msg) printError("Option '", option, "' requires an argument\n");
    return option::ARG_ILLEGAL;
  }
};
enum  optionType { DEBUG, QUIET, VERBOSE, SILENT,
		   IMIN, IMAX, JMIN, JMAX, ISTEP, JSTEP, ISPEC, JSPEC};
enum  optionIndex { UNKNOWN, HELP, PLUGINS, LIST, VERSION, VERBOSITY, NOSIGFPE, RANGE,
		    BOUNDARIES, STEPS, IPCT, TIME, TMIN, FOV, RESOLUTION,
		    DISTANCE, PALN, INCLINATION, ARGUMENT, NTHREADS, NPROCESSES,
		    SETPARAMETER, UNIT, XMLWRITE};
const option::Descriptor usage[] =
{
 {UNKNOWN, 0, "", "",option::Arg::None, "\nUSAGE: gyoto [options] input.xml output.fits\t\n\n"
                                        "Generic options:\t\n  -- \tStop option processing." },
 {HELP, 0,"h","help",option::Arg::Optional, "  --help[=<c>, -h<c>  \tWithout argument, print usage and exit. With argument, document class <c> (e.g. \"Screen\", \"Astrobj::Star\") and exit." },
 {VERSION, 0, "V", "version", option::Arg::None, "  --version, -V  \tPrint the Gyoto version."},
 {LIST, 0,"l","list",option::Arg::None, "  --list, -l  \tPrint the Gyoto register of Astrobj, Metrics etc." },
 {NOSIGFPE, 0, "", "no-sigfpe",option::Arg::None, "  --no-sigfpe \tDo not enable SIGFPE."
#if !defined HAVE_FENV_H
  " (noop: this Gyoto lacks fenv.h support)."
#endif
 },
 {PLUGINS, 0,"p","plugins",option::Arg::Optional, "  --plugins=<l>, -p<l>  \tList of plug-ins to load instead of $GYOTO_PLUGINS." },
 {NTHREADS, 0, "T", "nthreads", Gyoto::Arg::Required, "  --nthreads=<n>, -T<n> \tNumber of parallel threads to use."},
 {NPROCESSES, 0, "P", "nprocesses", Gyoto::Arg::Required, "  --nprocesses=<n>, -P<n> \tNumber of MPI parallel processes to use."},
 {IPCT, 0, "", "impact-coords", option::Arg::Optional, "  --impact-coords[=<f>] \tRead impact coordinates from file <f> or store in output.fits."},
 {XMLWRITE, 0, "X", "xmlwrite", Gyoto::Arg::Required, "  --xmlwrite=<f>, -X<f> \tWrite back scenery to XML file <f>. Useful to see default values and check the effect of --parameter, see below."},
 {UNKNOWN, 0, "", "",option::Arg::None, "\nVerbosity level:" },
 {VERBOSITY, SILENT, "s", "silent", option::Arg::None, "  --silent, -s \tBe silent." },
 {VERBOSITY, QUIET, "q", "quiet", option::Arg::None, "  --quiet, -q \tBe quiet." },
 {VERBOSITY, VERBOSE, "v", "verbose", option::Arg::Optional, "  --verbose[=<l>], -v[<l>] \tBe verbose. Optional parameter: verbosity level." },
 {VERBOSITY, DEBUG, "d", "debug", option::Arg::None, "  --debug, -d \tEnable debug output." },
 {UNKNOWN, 0, "", "",option::Arg::None, "\nField selection:" },
 {BOUNDARIES, IMIN, "", "imin", Gyoto::Arg::Required, "  --imin=<arg>  \tFirst column (1)."},
 {BOUNDARIES, IMAX, "", "imax", Gyoto::Arg::Required, "  --imax=<arg>  \tLast column (ULONG_MAX)."},
 {STEPS, ISTEP, "", "di", Gyoto::Arg::Required, "  --di=<arg>  \tColumn step (1)."},
 {BOUNDARIES, JMIN, "", "jmin", Gyoto::Arg::Required, "  --jmin=<arg>  \tFirst line (1)."},
 {BOUNDARIES, JMAX, "", "jmax", Gyoto::Arg::Required, "  --jmax=<arg>  \tLast line (ULONG_MAX)."},
 {STEPS, JSTEP, "", "dj", Gyoto::Arg::Required, "  --dj=<arg>  \tLine step (1)."},
 {RANGE, ISPEC, "i", "ispec", Gyoto::Arg::Required, "  --ispec=<arg>, -i<arg>  \tColumn specification (imin[:imax[:di]])."},
 {RANGE, JSPEC, "j", "jspec", Gyoto::Arg::Required, "  --jspec=<arg>, -j<arg>  \tLine specification (jmin[:jmax[:dj]])."},
 {UNKNOWN, 0, "", "",option::Arg::None, "\nScreen parameters:" },
 {TIME, 0, "", "time", Gyoto::Arg::Required, "  --time=<arg> \tObserving date."},
 {TMIN, 0, "", "tmin", Gyoto::Arg::Required, "  --tmin=<arg> \tMinimum time."},
 {FOV, 0, "", "fov", Gyoto::Arg::Required, "  --fov=<arg> \tField-of-view."},
 {RESOLUTION, 0, "r", "resolution", Gyoto::Arg::Required, "  --resolution=<n>, -r<n> \tField size in pixels (on each side)."},
 {DISTANCE, 0, "", "distance", Gyoto::Arg::Required, "  --distance=<arg> \tDistance from observer."},
 {PALN, 0, "", "paln", Gyoto::Arg::Required, "  --paln=<arg> \tPosition angle of the line of nodes."},
 {INCLINATION, 0, "", "inclination", Gyoto::Arg::Required, "  --inclination=<arg> \tInclination."},
 {ARGUMENT, 0, "", "argument", Gyoto::Arg::Required, "  --argument=<arg> \tArgument of the x axis."},
 {UNKNOWN, 0, "", "",option::Arg::None, "\nArbitrary parameters:" },
 {UNIT, 0, "u", "unit", Gyoto::Arg::Optional, "  --unit[=<u>], -u[<u>] \tUnit for following parameters (until next instance of this option)."},
 {SETPARAMETER, 0, "E", "parameter", Gyoto::Arg::Required,
  "  --parameter=<Name>[=<value>],      -E<Name>[=<value>]"
  "\tSet arbitrary parameter by name. Optional value is expressed in unit previously set with --unit/-u. Examples: -ENThreads=5, -EAstrobj::Spectrum::Temperature=100."},
 {0,0,0,0,0,0}
};

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

#define ERROR_GENERIC          1
#define ERROR_INITIALIZING     2
#define ERROR_READING_SCENERY  3
#define ERROR_RAYTRACING       4
#define ERROR_MK_VIDEO 5

static std::string curmsg = "";
static int curretval = 1;

void gyotoErrorHandler( const Gyoto::Error e ) {
  cerr << curmsg << e << endl;
  if (debug()) abort(); // to keep stack for debugger
  exit (curretval);
}

static void gyotoVersion() {
  cout << " Copyright (c) 2011-2019 Frédéric Vincent, Thibaut Paumard,\n"
       << "                         Odele Straub and Frédéric Lamy.\n"
       << " GYOTO is distributed under the terms of the GPL v. 3 license.\n"
       << " We request that use of Gyoto in scientific publications be "
       << " properly \n acknowledged. Please cite:\n"
       << "  GYOTO: a new general relativistic ray-tracing code,\n"
       << "  F. H. Vincent, T. Paumard, E. Gourgoulhon & G. Perrin 2011,\n"
       << "  Classical and Quantum Gravity 28, 225011 (2011) "
       << "[arXiv:1109.4769]"
       << endl << endl;
}

int main(int argc, char** argv) {

  // Set-up error reporter
  Gyoto::Error::setHandler ( &gyotoErrorHandler );

  // No debug out put by default
  debug(0);
  
  char * parfile=NULL;
  string ipctfile="";
  string param;

  size_t imin=1, imax=ULONG_MAX, jmin=1, jmax=ULONG_MAX;
  ptrdiff_t di=1, dj=1;
  bool  ipct=0;
  long  ipctdims[3]={0, 0, 0};
  double ipcttime;

  string pluglist= getenv("GYOTO_PLUGINS")?
    getenv("GYOTO_PLUGINS"):
    GYOTO_DEFAULT_PLUGINS;

  argc-=(argc>0); argv+=(argc>0); // skip program name argv[0] if present

  // if first argument is exactly mk-video, make a video!
  if ( (argc>0) && (!strcmp(argv[0], "mk-video")) ) {
    curmsg = "In gyoto.C: in mk-video: ";
    curretval = ERROR_MK_VIDEO;
    GYOTO_DEBUG << "trying to load python plugin\n";
    void* handle=NULL;
    int (*mk_video)(int, char**) = NULL;
    Gyoto::Register::init(NULL);
    std::vector< std::string > plugnames = {"python3",
					    "python3.7", "python3.6", "python3.5",
					    "python2.7",
					    "python3.8", "python3.9"};
    for (size_t k=0; k < plugnames.size(); ++k) {
      GYOTO_DEBUG << "trying to load plug-in " << plugnames[k] << endl;
      handle = loadPlugin(plugnames[k].c_str(), 2);
      if (handle) {
	GYOTO_DEBUG << "trying find symbol \"mk_video\" in plug-in " << plugnames[k] << endl;
	mk_video = (int (*)(int, char**)) dlsym(handle, "mk_video");
	if (mk_video) break;
      }
    }
    if (!mk_video) GYOTO_ERROR("No Python plug-in containing mk_video() found");

    return mk_video(argc, argv);
  }

  // else parse arguments
  option::Stats  stats(true, usage, argc, argv);
  option::Option* options = new option::Option[stats.options_max];
  option::Option* buffer  = new option::Option[stats.buffer_max];
  option::Parser parse(true, usage, argc, argv, options, buffer, 1);

  if (parse.error())
    return 1;

  // Check whether to output usage string
  if (options[VERSION]
      || (!options[LIST] && (argc == 0 || parse.nonOptionsCount() != 2 || options[UNKNOWN])) || options[HELP] ) {
    cout << GYOTO_STRINGIFY(PACKAGE_STRING) << endl
	 << "ABI compatibility version: " << GYOTO_SOVERS << endl;
    if (!options[VERSION] && !options[HELP].arg) option::printUsage(std::cout, usage);
    if (!options[VERSION] && !options[LIST] && !options[HELP].arg) {
      if (options[HELP]) return 0;
      return 1;
    }
  }

  // Process options setting verbosity or debug level
  for (option::Option* opt = options[VERBOSITY]; opt; opt = opt->next()) {
    switch (opt->type()) {
    case DEBUG: debug(1); break;
    case SILENT: verbose(0); break;
    case QUIET: verbose(GYOTO_QUIET_VERBOSITY); break;
    case VERBOSE:
      if (opt->arg) verbose(atoi(opt->arg));
      else verbose(10);
      break;
    default:
      cerr << "gyoto: error parsing command line, unknown verbosity type" << endl;
      return 1;
    }
  }

  // Retrieve file names
  if (parse.nonOptionsCount() > 0) parfile=strdup(parse.nonOptions()[0]);
  if (parse.nonOptionsCount() > 1) pixfile=strdup(parse.nonOptions()[1]);

  if (options[VERSION]) {
    gyotoVersion();
    return 0;
  }

  if (options[PLUGINS])
    pluglist = options[PLUGINS].last()->arg?options[PLUGINS].last()->arg:"";
  curmsg = "In gyoto.C: Error initializing libgyoto: ";
  curretval = ERROR_INITIALIZING;
  Gyoto::Register::init(pluglist.c_str());

  SmartPointer<Scenery> scenery = NULL;
  if (parfile) {
    Factory *factory =NULL;
    if (verbose() >= GYOTO_QUIET_VERBOSITY) cout << "Reading parameter file: " << parfile << endl;
    curmsg = "In gyoto.C: Error in Factory creation: ";
    curretval = ERROR_READING_SCENERY;
    factory = new Factory(parfile);

    curmsg = "In gyoto.C: Error getting Kind: ";
    const string kind = factory->kind();

    if (kind.compare("Scenery")) {
      cerr << "Unknown kind for root element in XML file" << endl;
      return 1;
    }

    curmsg = "In gyoto.C: Error getting Scenery: ";
    scenery = factory -> scenery();

    curmsg = "In gyoto.C: Error deleting Scenery: ";
    delete factory;
  }

  if (options[HELP] && options[HELP].arg) {
    help(options[HELP].arg);
    if (!options[LIST]) return 0;
  }

  if (options[LIST]) {
    Gyoto::Register::list();
    return 0;
  }

  curmsg = "In gyoto.C: Error initializing ray-tracing: ";
  curretval = ERROR_RAYTRACING;
  SmartPointer<Screen>  screen = scenery->screen();
  string unit="";

  for (int i = 0; i < parse.optionsCount(); ++i) {
    option::Option& opt = buffer[i];
    switch(opt.index()) {
    case NOSIGFPE:
#if defined HAVE_FENV_H
      GYOTO_DEBUG << "enabling SIGFPE delivery\n";
      feenableexcept(FE_DIVBYZERO | FE_OVERFLOW | FE_INVALID);
#endif
      break;
    case BOUNDARIES:
      {
	double valtest=Gyoto::atof(opt.arg);
	if (valtest<=0){
	  cerr << "In gyoto.C: screen indices should be >0" << endl;
	  return 1;
	}
	size_t val=atoi(opt.arg);
	switch (opt.type()) {
	case IMIN: imin=val; break;
	case IMAX: imax=val; break;
	case JMIN: jmin=val; break;
	case JMAX: jmax=val; break;
	default:
	  cerr << "Gyoto BUG: unknown type of screen boundary" << endl;
	  return 1;
	}
      }
      break;
    case STEPS:
      {
	ptrdiff_t val=atoi(opt.arg);
	if (val==0) val=1;
	switch(opt.type()) {
	case ISTEP: di=val; break;
	case JSTEP: dj=val; break;
	default:
	  cerr << "Gyoto BUG: unknown type of screen step" << endl;
	  return 1;
	}
      }
      break;
    case RANGE:
      {
	string spec=opt.arg;
	size_t pos=spec.find(":"), pos2=string::npos;
	string sub=spec.substr(0, pos);
	size_t nmin=1, nmax=ULONG_MAX;
	ptrdiff_t dn=1;
	if (sub.length()) nmin=atoi(sub.c_str());
	if (pos==string::npos) {
	  nmax=nmin;
	} else {
	  pos2=spec.find(":", pos+1);
	  sub=spec.substr(pos+1, pos2-pos-1);
	  if (sub.length()) nmax=atoi(sub.c_str());
	  if (pos2 != string::npos) {
	    sub=spec.substr(pos2+1);
	    if (sub.length()) dn=atoi(sub.c_str());
	  }
	}
	GYOTO_DEBUG << "nmin="<<nmin<<", nmax="<<nmax<<", dn="<<dn<<endl;
	switch (opt.type()) {
	case ISPEC: imin=nmin; imax=nmax; di=dn; break;
	case JSPEC: jmin=nmin; jmax=nmax; dj=dn; break;
	default:
	  cerr << "Gyoto BUG: unknown type of pixel specification\n";
	  return 1;
	}
      }
      break;
    case IPCT:
      if (opt.arg) ipctfile=opt.arg;
      else ipct=1;
      break;
    case TIME:        screen -> time       (Gyoto::atof(opt.arg)); break;
    case TMIN:       scenery -> tMin       (Gyoto::atof(opt.arg)); break;
    case FOV:         screen -> fieldOfView(Gyoto::atof(opt.arg)); break;
    case RESOLUTION:  screen -> resolution (       atoi(opt.arg)); break;
    case DISTANCE:    screen -> distance   (Gyoto::atof(opt.arg)); break;
    case PALN:        screen -> PALN       (Gyoto::atof(opt.arg)); break;
    case INCLINATION: screen -> inclination(Gyoto::atof(opt.arg)); break;
    case ARGUMENT:    screen -> argument   (Gyoto::atof(opt.arg)); break;
    case NTHREADS:   scenery -> nThreads   (       atoi(opt.arg)); break;
    case NPROCESSES: scenery -> nProcesses (       atoi(opt.arg)); break;
    case UNIT: unit=opt.arg?opt.arg:""; break;
    case SETPARAMETER:
      {
	string arg=opt.arg;
	size_t pos=arg.find("=");
	string name=arg.substr(0, pos);
	string val=(pos==string::npos)?"":arg.substr(pos+1);
	GYOTO_DEBUG << "Setting parameter \"" << name << "\" to value \"" << val << "\" using unit \"" << unit << "\".\n";
	if(scenery -> setParameter(name, val, unit))
	  throwError("Unknown parameter");
      }
      break;
    case XMLWRITE: Factory(scenery).write(opt.arg); break;
    default: break;
    }
  }

#if defined HAVE_MPI
  if (scenery -> nProcesses() || getenv("OMPI_COMM_WORLD_SIZE")) {
    int status = MPI_Init(&argc, &argv);
    if (status) {
      cerr << "error initializing MPI"<< endl;
      return 2;
    }
    // If WORLD size is more than 1, use processes from WORLD
    // instead of spawning new processes
    int wsize=0, rank=0;
    MPI_Comm_size(MPI_COMM_WORLD, &wsize);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (wsize > 1) {
      if (rank==0) {
	GYOTO_INFO << "Process rank " << rank << " becoming manager\n";
	scenery -> nProcesses(-1);
      } else {
	GYOTO_INFO << "Process rank " << rank << " becoming worker\n";
	curmsg = "In gyoto.C: error in MPI worker: ";
	Scenery::mpiWorker();
	int started , stopped , error ;
	error = MPI_Initialized ( & started ) ;
	error = MPI_Finalized ( & stopped ) ;
	if (started && !stopped) MPI_Finalize();
	return 0;
      }
    }
  }
#endif

  // State copyright
  if (verbose() >= GYOTO_QUIET_VERBOSITY) {
    gyotoVersion();
  }

  {
    double tobs= screen -> time();
    size_t res = screen -> resolution();

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

      double dt = tobs * GYOTO_C / scenery -> metric() -> unitLength()
	- ipcttime;
      for (size_t i=0; i < ipctnelt; i+=8)
	if (impactcoords[i] != DBL_MAX) impactcoords[i] += dt;
      ipcttime = tobs * GYOTO_C / scenery -> metric() -> unitLength();
    }

    Quantity_t quantities = scenery -> getRequestedQuantities();
    if (debug()) cerr << "DEBUG: Gyoto.C: Requested Quantities: "
		      << quantities <<endl;

    size_t nbnuobs=0;
    if (quantities & GYOTO_QUANTITY_SPECTRAL) {
      SmartPointer<Spectrometer::Generic> spr = screen -> spectrometer();
      if (!spr) throwError("Spectral quantity requested but "
			   "no spectrometer specified!");
      nbnuobs = spr -> nSamples();
    }
               //nb of frames that will be used for spectral cube
    size_t nbdata= scenery->getScalarQuantitiesCount()
      +scenery->getSpectralQuantitiesCount()*nbnuobs;
               //nb of frames used for diverse interesting outputs
               //(obs flux, impact time, redshift..)
    size_t nelt=res*res*nbdata;
    vect = new double[nelt];

    // First check whether we can open file
    int naxis=3; 
    long naxes[] = {long(res), long(res), long(nbdata)};
    nelements=nelt; 

    fits_create_file(&fptr, pixfile, &status);
    fits_create_img(fptr, DOUBLE_IMG, naxis, naxes, &status);
    fits_report_error(stderr, status);
    if (status) return status;

    // Allocate space for the output data
    ::data = new Astrobj::Properties();
    ::data->alloc=true;

    size_t curquant=0;
    size_t offset=res*res;

    if (debug()) {
      cerr << "DEBUG: gyoto.C: flag_radtransf = ";
      cerr << scenery -> astrobj() -> opticallyThin() << endl;
      cerr << "DEBUG: gyoto.C: Requested quantities: ";
      cerr << scenery -> requestedQuantitiesString() << endl;
    }

    char keyname[FLEN_KEYWORD];
    char const * fmt="QUANT_%lu";
    char * CNULL=NULL;

    if (quantities & GYOTO_QUANTITY_INTENSITY) {
      ::data->intensity=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("Intensity"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_EMISSIONTIME) {
      ::data->time=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("EmissionTime"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_MIN_DISTANCE) {
      ::data->distance=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("MinDistance"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_FIRST_DMIN) {
      ::data->first_dmin=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("FirstDistMin"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_REDSHIFT) {
      if (debug())
	cerr << "DEBUG: gyoto.C: REDSHIFT requested\n";
      ::data->redshift=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("Redshift"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_NBCROSSEQPLANE) {
      if (debug())
	cerr << "DEBUG: gyoto.C: NBCROSSEQPLANE requested\n";
      ::data->nbcrosseqplane=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("NbCrossEqPlane"),
		     CNULL, &status);
    }
    if ((quantities & GYOTO_QUANTITY_IMPACTCOORDS || ipct) && !ipctdims[0] ) {
      // Allocate if requested AND not provided
      cerr << "gyoto.C: allocating data->impactcoords" << endl;
      ::data->impactcoords = impactcoords = new double [res*res*16];
      ipcttime = tobs * GYOTO_C / scenery -> metric() -> unitLength();
    }
    if (quantities & GYOTO_QUANTITY_USER1) {
      ::data->user1=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("User1"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_USER2) {
      ::data->user2=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("User2"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_USER3) {
      ::data->user3=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("User3"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_USER4) {
      ::data->user4=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("User4"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_USER5) {
      ::data->user5=vect+offset*(curquant++);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("User5"),
		     CNULL, &status);
    }
    double * curvect=vect+offset*curquant;
    if (quantities & GYOTO_QUANTITY_SPECTRUM) {
      ::data->spectrum=curvect;
      curvect += offset*nbnuobs;
      ++curquant;
      ::data->offset=int(offset);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("Spectrum"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_SPECTRUM_STOKES_Q) {
      ::data->stokesQ=curvect;
      curvect += offset*nbnuobs;
      ++curquant;
      ::data->offset=int(offset);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("SpectrumStokesQ"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_SPECTRUM_STOKES_U) {
      ::data->stokesU=curvect;
      curvect += offset*nbnuobs;
      ++curquant;
      ::data->offset=int(offset);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("SpectrumStokesU"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_SPECTRUM_STOKES_V) {
      ::data->stokesV=curvect;
      curvect += offset*nbnuobs;
      ++curquant;
      ::data->offset=int(offset);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("SpectrumStokesV"),
		     CNULL, &status);
    }
    if (quantities & GYOTO_QUANTITY_BINSPECTRUM) {
      ::data->binspectrum=curvect;
      curvect += offset*nbnuobs;
      ++curquant;
      ::data->offset=int(offset);
      sprintf(keyname, fmt, curquant);
      fits_write_key(fptr, TSTRING, keyname,
		     const_cast<char*>("BinSpectrum"),
		     CNULL, &status);
    }
    
    signal(SIGINT, sigint_handler);

    curmsg = "In gyoto.C: Error during ray-tracing: ";

    if (imax>res) imax=res;
    if (jmax>res) jmax=res;

    Screen::Range irange(imin, imax, di);
    Screen::Range jrange(jmin, jmax, dj);
    Screen::Grid  grid(irange, jrange, "\rj = ");

    if (verbose() >= GYOTO_QUIET_VERBOSITY)
      cout << "j = " << 1 << "/" << (jmax-jmin)/dj+1 << flush;

    scenery -> rayTrace(grid, ::data, ipctdims[0]?impactcoords:NULL);

    curmsg = "In gyoto.C: Error while saving: ";
    if (verbose() >= GYOTO_QUIET_VERBOSITY)
      cout << "\nSaving to file: " << pixfile << endl;
    signal(SIGINT, SIG_DFL);


    // Save to fits file
    fits_write_pix(fptr, TDOUBLE, fpixel, nelements, vect, &status);

    if (quantities & GYOTO_QUANTITY_IMPACTCOORDS || ipct) {
      // Save if requested, copying if provided
      cout << "Saving precomputed impact coordinates" << endl;
      long naxes_ipct[] = {16, long(res), long(res)};
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

    GYOTO_DEBUG << "screen = NULL\n";
    screen = NULL;

#if defined HAVE_MPI
    GYOTO_DEBUG << "scenery->mpiTerminate()\n";
    scenery->mpiTerminate();
#endif

    GYOTO_DEBUG << "scenery = NULL\n";
    scenery = NULL;

    if (status) return status;

  }

#if defined HAVE_MPI
  int started , stopped , error ;
  error = MPI_Initialized ( & started ) ;
  error = MPI_Finalized ( & stopped ) ;
  if (started && !stopped) MPI_Finalize();
#endif

  return 0;
}
