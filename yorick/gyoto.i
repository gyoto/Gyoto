plug_in, "gyoto";
/*
    Copyright 2011-2016 Thibaut Paumard

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

#include "graphk.i"

local GYOTO_PLUGINS;
/* DOCUMENT GYOTO_PLUGINS

    Comma-separated list of plug-ins to load when initializing gyoto
    (more plug-ins can be loaded afterwards). Must be set before
    #including gyoto.i. If this Yorick variable is not set, then the
    environment variable of the same name is used. If it is not set
    either, then the built-in default (usually
    "stdplug,nofail:lorene") is used. Each plug-in may be prefixed
    with "nofail:" to tell gyoto that it is not an error if this
    plug-oin cannot be loaded.

   EXAMPLES:
    To load gyoto with no plug-in at all, and load "myplugin" afterwards:
    GYOTO_PLUGINS="";
    GYOTO_NO_STD=1;
    #include "gyoto.i"
    noop, gyoto.requirePlugin("myplugin");

   SEE ASLO: gyoto.requirePlugin, GYOTO_NO_STD
 */

local GYOTO_PLUGINS;
/* DOCUMENT GYOTO_NO_STD=1

    By default, #include "gyoto.i" will also #include "gyoto_std.i"
    (which in turn loads the stdplug Gyoto plug-in).  To inhibit this
    behavios, set this variable (to something that evaluates as true)
    before #including gyoto.i.

   EXAMPLES:
    To load gyoto with no plug-in at all, and load only "myplugin":
    GYOTO_PLUGINS="";
    GYOTO_NO_STD=1;
    #include "gyoto.i"
    noop, gyoto.requirePlugin("myplugin");

   SEE ASLO: gyoto.requirePlugin, GYOTO_PLUGINS
 */

extern gyoto_haveXerces;
/* DOCUMENT have_xerces = gyoto.haveXerces()
    Tell whether GYOTO was compiled with Xerces support (XML i/o)
   OUTPUT:
    HAVE_XERCES=1 if compiled with Xerces, else 0.
*/

extern gyoto_haveCFITSIO;
/* DOCUMENT have_cfitsio = gyoto.haveCFITSIO()
    Tell whether GYOTO was compiled with CFITSIO support (FITS i/o)
   OUTPUT:
    HAVE_CFITSIO=1 if compiled with CFITSIO, else 0.
*/

extern gyoto_haveBoost;
/* DOCUMENT have_boost = gyoto.haveBoost()
    Tell whether GYOTO was compiled with Boost support (advanced integrators)
   OUTPUT:
    HAVE_BOOST=1 if compiled with Boost, else 0.
*/

extern gyoto_haveUDUNITS;
/* DOCUMENT have_udunits = gyoto.haveUNUITS()
    Tell whether GYOTO was compiled with UDUNITS support (advanced units)
   OUTPUT:
    HAVE_UDUNITS=1 if compiled with UDUNITS, else 0.
*/

extern gyoto_havePTHREAD;
/* DOCUMENT have_pthread = gyoto.havePTHREAD()
    Tell whether GYOTO was compiled with POSIX thread support
    (multi-thread parallel computing)
   OUTPUT:
    HAVE_PTHREAD=1 if compiled with PTHREAD, else 0.
*/

extern gyoto_haveMPI;
/* DOCUMENT have_mpi = gyoto.haveMPI()
    Tell whether GYOTO was compiled with MPI support (multi-process
    parallel computing)
    
   OUTPUT:
    HAVE_MPI=1 if compiled with MPI, else 0.

   SEE ALSO: gyoto.MPI_Init, .MPI_Finalize, .MPI_Initialized, .MPI_Finalized
*/

extern gyoto_MPI_Init;
/* DOCUMENT status = gyoto.MPI_Init([argv]);
    Initialize MPI.
    
    Wrapper around MPI_Init(numberof(argv), argv). argv may be
    modified upon completion.

    Returns 0 for success, no-zer otherwise.

    If MPI is not compiled-in, returns 1.

   SEE ALSO: gyoto.haveMPI, .MPI_Finalize, .MPI_Initialized, .MPI_Finalized
*/

extern gyoto_MPI_Finalize;
/* DOCUMENT gyoto.MPI_Finalize;
    Finalize MPI.
    
    Wrapper around MPI_Finalize(). Returns 0 for success, non-zero otherwise.
    
    If MPI is not compiled-in, returns 1.

   SEE ALSO: gyoto.haveMPI, .MPI_Init, .MPI_Finalized, .MPI_Initialized
*/

extern gyoto_MPI_Finalized;
/* DOCUMENT is_finalized=gyoto.MPI_Finalized();
    Tell whether some implemention of MPI_Finalize() was already called.

    If MPI support is not present, return 0.

   SEE ALSO: gyoto.haveMPI, .MPI_Finalize, .MPI_Init, .MPI_Initialzed
 */

extern gyoto_MPI_Initialized;
/* DOCUMENT is_initialized=gyoto.MPI_Initialized();
    Tell whether some implemention of MPI_Init() was already called.

    If MPI support is not present, return 0.

   SEE ALSO: gyoto.haveMPI, .MPI_Init, .MPI_Finalize, .MPI_Finalized
 */

extern __gyoto_setErrorHandler;
/* xDOCUMENT __gyoto_setErrorHandler
   Must be called once to attach the GYOTO error handler to Yorick's one
*/
__gyoto_setErrorHandler;

extern gyoto_requirePlugin;
extern gyoto_loadPlugin;
extern gyoto_havePlugin;
/* DOCUMENT gyoto.requirePlugin, plugin[, plugin2[, plugin3]] [, nofail=1]
            gyoto.loadPlugin, plugin[, plugin2[, plugin3]] [, nofail=1]
            gyoto.havePlugin(plugin]

   Load Gyoto plug-ins.

  INPUTS:
   gyoto.requirePlugin() and gyoto.loadPlugin() accept an aribtrary
   number of positional arguments, each a string or string array
   naming individual Gyoto plug-ins. For instance, all of the
   following attempt to load the plug-ins, stdplug, lorene and myplug:

    gyoto.requirePlugin, "stdplug", "lorene", myplug"
    gyoto.requirePlugin, ["stdplug", "lorene", myplug"]
    gyoto.requirePlugin, "stdplug", ["lorene", myplug"]

   In the first (preferred) form, each plug-in is loaded only if it is
   not already. In the second form, each plug-in is forcibly loaded,
   which has the effect that the classes it implements will be
   registered twice. The only use of this feature is to put the
   classes of this plug-in first in the registers to make sure they
   are the version that is used when instanciating a class without
   specifying it's plug-in. In short, this should not be used.

   gyoto.havePlugin("myplug") returns 1 if myplug has already been
   loaded, 0 otherwise.

  KEYWORDS:
   nofail= if set and true, failure to load a plug-in will not trigger
           an error. It applies to _all_ plug-ins in the list.

  EXAMPLE:
   gyoto.requirePlugin, "stdplug"
   gyoto.requirePlugin, "lorene", nofail=1
 */

extern __gyoto_initRegister;
/* xDOCUMENT __gyoto_initRegister[, pluglist]

   Must be called once to initialize the GYOTO plug-in register.

   PLUGLIST: list of plug-ins to load. If absent or nil, load the
   default list (either from the GYOTO_PLUGINS environment variable,
   or built-in default).
*/
__gyoto_initRegister, GYOTO_PLUGINS;

extern __gyoto_exportSupplier;
/* xDOCUMENT __gyoto_exportSupplier
   To be called by client plug-ins like gyoto_std.i
*/


require, "pl3d.i";
#include "gyoto_constants.i"

local gyoto;
/* DOCUMENT GYOTO -- General relativtY Orbit Tracer of Observatoire de paris

    GYOTO is a Yorick plug-in, based on an eponym C++ library, which
    allows computing geodesics in curved space-time. The two main uses
    are:
      -- computing the orbits of stars around relativistic objects;
      -- producing ray-traced images of compact objects and their
         surrondings.

    CONCEPTS AND SYNTAX
    ===================

    The Yorick plug-in implementation is object-oriented, reflecting
    that of the underlying C++ library. For the most part, it adheres
    to the Yorick uses, but some details can be surprising to a Yorick
    user.

    The gyoto namespace
    -------------------

    Most of the Gyoto functions are put into the "gyoto" namespace
    (actually an oxy group, for the interested reader). You can
      restore, gyoto;
    and be done with the "gyoto." prefix. A few (rare) functions have
    the "gyoto_" prefix, they are not in the namespace.


    Creating a GYOTO object
    -----------------------

    To create a GYOTO object, one calls one of the object creators:
    for instance
       gg = gyoto.KerrBL() ;
    results in GG being an opaque object containing a GYOTO Metric (a
    Kerr Metric in this case, using the Boyer-Lindquist coordinate
    system).

    The creators accept a filename as positional argument (see FILE
    FORMAT below): gg = gyoto.KerrBL( "my_KerrBL_descritpion.xml" );

    Arbitrary derived classes can also be instantiated using the base
    class constructor:
       gg = gyoto.Metric("KerrBL", [plugin]);
    The PLUGIN argument is optional and has two effects:
       - it ensures the plug-in PLUGIN is loaded before attempting to
         instantiate the object;
       - it ensures the class of the newly created object comes from
         this plug-in and that another class by the same name in
         another plug-in is not inadvertently used.

    Finally, it is possible to instantiate an object using a buffer
    containing XML data:
       #include "pathfun.i"
       gg = gyoto.KerrBL();
       xmldata = pathform(print(gg), delim="\n");
       gg2 =gyoto.Metric(xmldata).
    This should be equivalent to gg2 == gg.clone(), see COPYING
    vs. CLONIG below.

    GYOTO OBJECTS BEHAVE LIKE FUNCTIONS
    -----------------------------------

    Once a GYOTO object has been created, it can be accessed in an
    object-like fashion (the object then ressembles a Yorick function)
    or through the same function which was used to create it. For
    instance, to change the spin parameter of the GG Metric above, the
    two following are exactly synonymous:
       gg, spin=0.8 ;
       gyoto.KerrBL, gg, spin=0.8;

    Likewise, the value of the spin member can be retrieved like this:
       a = gg(spin=);
       a = gyoto.KerrBL(gg, spin=);

    Beware: in the second form, the positional parameter GG must be
    given before any keyword parameter. This is one of the few points
    were the GYOTO syntax does not fully adheres the Yorick uses.

    SYNTACTIC SUGAR: the dot operator
    ---------------------------------

    The normal way to access members or methods in a Gyoto object (see
    below) is:
       object, member1=val1, member2=val2...;   // set member1 and member2
       retval = object(member1=, member2=val2); // retrieve member1, set member2
    In order to get a more usual look-and-feel, it is also possible to
    write this instead:
       noop, object.member(val[,other_args]); // set member
       retval = object.member([,other_args]); // get member
    For instance, to retrieve the mass in a Metric object, these two
    syntaxes work (and will continue working in future releases, until
    a design change is necessary):
       m = gg(mass=, unit="sunmass");
       m = gg.mass(unit="sunmass");

    Incidently, it is currently possible to process several keywords
    also in this way:
       noop, object.member1(value1, member2=value2);
    is the same as
       object, member1=value1, member2=value2;
    This is not intentional and your code should not rely on this
    behaviour. For the sake of forward compatibility, when using the
    .member() notation, use only additional keywords if they are
    needed to process .member. Currently, this concerns only the
    "unit=" keyword.

    Additional sugar: when chaining dot extraction, intermediate "()"
    are implied, so that
       scenery.screen.resolution();
    is a shorthand for:
       scenery.screen().resolution();

    A note on optimisation: the dot operator (currently) works by
    calling the more basic keyword= syntax. It is therefore very
    slightly slower.

    In the following, we give the alternative dot notation alternative
    as a comment (following "//").
    
    MEMBERS
    -------
    
    All the GYOTO objects have properties called "members" such as the
    spin of a Kerr metric. Members can be set or retrieved using the
    corresponding keywords, either upon object instantiation or at any
    later point in time:
       gg = gyoto.KerrBL( spin=0.995, mass=4e6 ) ;
       gg, spin=0.5;
       spin_param = gg(spin=); // or spin_param = gg.spin()
       
    Note how giving a member keyword without any value, as in the last
    example above, allows _retrieving_ the previously set value.

    When setting member, it is also possible to call the object as a
    function. In that case, the return value will be the object
    itself, allowing to call it again as a function:
       spin_val = gg(spin=0.5)(spin=)
       // or spin_val = gg.spin(0.5).spin()
    Although the above example is trivial, this is useful in many
    slightly more complex situations:
       Getting the resolution in a Screen attached to a Scenery:
         res_val = scenery(screen=)(resolution=);
         // or res_val = scenery.screen.resolution();
       Setting resolution:
         noop, scenery.screen(resolution=res_val);
         // or noop, scenery.screen.resolution(res_val);
       (The noop above is not mandatory but avoid spurious display).
       
    Some member keywords accept more than one parameter, separated by
    comas. This is the second exception GYOTO makes to the Yorick
    syntax. For example, for setting the initial position and velocity
    of a star, one can use:
       st = gyoto.Star( initcoord=pos,vel );
       // or st = gyoto.Star().initcoord(pos, vel);
    Only one such keyword can be set at any time because it would be
    exceedingly difficult to parse them otherwise.
    
       
    COPYING vs. CLONING
    -------------------
       
    Be careful that GYOTO objects behave like pointers: if you make a
    copy of an object, any change in one of the instances will affect
    both. For instance:
       gg2 = gg;
       gg2, spin=0.5;
    In the above, the spin parameter of GG will be the same as that of
    GG2 (0.5). In fact, in memory, GG and GG2 point to the same C++
    object. This is on purpose: this way, if you affect a Metric to a
    Photon and a Star, then change the spin of the Metric, the change
    will be reflected in all of the objects linked with this Metric:
       gg = gyoto.KerrBL();
       st = gyoto.Star( metric=gg );
       ph = gyoto.Photon( metric=gg );
       gg, spin=0.5;
    In the above, the spin change affect the Photon and the Star as
    well.

    If what you want is a detached copy of a Gyoto object, you can get
    it with the CLONE method. Following the above example:
    
       gg_copy = gg;   // This is copying: gg2 and gg are the same
                       // object
       
       gg_clone = gg(clone=);        // This is cloning: gg_clone is a
       // or gg_clone = gg.clone();  // detached copy

       gg_copy, spin=0.2;
       gg_clone, spin=0.7;

       gg.spin();  // the spin in gg is the same as in gg_copy, not
                   // gg_clone.

    
    METHODS
    -------

    Some keywords are "subroutine-like" methods in that they can have
    side-effects, but do not return a specific value. One notable
    example is the XMLWRITE keyword which most objects accept for
    dumping a description of themselves to an XML file:
       gg, xmlwrite="filename.xml";
       // or noop, gg.xmlwrite("filename.xml")
    Another notable example is the HELP keyword, which lists the
    properties that can be set or retrieved in any object:
       gg, help=;
       // or noop, gg.help();

    Other methods keywords are function-like: they usually take one or
    several parameter and return a value. Only one value-returning
    keyword can be set at a time (be it a function-like method or
    member keyword set for retrieving a value, such as "spin="):
       coor = gg(makecoord=yinit, cst);
       // or coor = gg.makecoord(init, cst);
       txyz = star(get_cartesian=dates);
       // or txyz = star.get_cartesian(dates);

    Some objects have a default methods that is called when no keyword
    is present. This is the case for metrics, which return the metric
    coefficients at a specific location, and of sceneries, which
    perform ray-tracing:
       coefficient = gg(coordinates, mu, nu);
       // or coefficients = scenery.metric()(coordinates, mu, nu);
       image = scenery();
       
       
    GYOTO OBJECT TYPES
    ==================

    gyoto.Metric: the general relativity metric in which objects move...
          Notable sub-classes: gyoto.KerrBL, gyoto.KerrKS
    
    gyoto.Photon: massless particules used for ray-tracing
    
    gyoto.Astrobj: astrophysical objects.
          Notable sub-classes: gyoto.Star, gyoto.FixedStar,
          gyoto.PolishDoughnut, gyoto.ThinInfiniteDisk

    gyoto.Screen: the camera for ray-tracing

    gyoto.Spectrometer: the spectral capabilites of a gyoto.Screen
          (also found in gyoto.Photon)

    gyoto.Spectrum: a spectrum, only used in gyoto.Star so far
    
    gyoto.Scenery: relationship between all of the above

    
    FILE FORMAT
    ===========

    Most objects can be described in XML format. That's what you see
    when you print an object. You can also read/write an object
    from/to and XML file:
       sc = gyoto.Scenery( "description.xml" );
       print, sc;
       sc, xmlwrite="backup.xml";

       
    GYOTOY
    ======

    GYOTO comes with a little GUI program: gyotoy. It allows
    interactively selecting the initial conditions for displaying the
    trajectory of a star or of a photon around a Kerr black-hole.

    
    EXAMPLES
    ========

    Compute and display the orbit of a star:
       data = gyoto.Star(radius=0.5,
               metric=(gg=gyoto.KerrBL(spin=0.995)),
               initcoord=[0, 10.791, 1.5708, 0], [0, 0, 0.0166637],
               xfill=800
               )(
                 get_skypos=gyoto.Screen(metric=gg)
                 );
       plg, data(,2), data(,1);

    Ray-trace a scenery (orbiting star near a Kerr black-hole):
       sc = gyoto.Scenery(
             metric=(gg=gyoto.KerrBL()),
             screen=gyoto.Screen(metric=gg,
                                 observerpos=[1000, 100., 0.15, 0.],
                                 fov=pi/10.,
                                 resolution=128),
             astrobj=gyoto.Star(metric=gg,
                                 radius=0.5,
                                 initcoord=[600, 9, pi/2, 0], [0, 0, 0.037037])
            )
       pli, sc(,,"Intensity");

     Trace the trajectory of a photon in the secondary image of the above:
       ph = gyoto.Photon(initcoord=sc,77,45, xfill=870.623);
       txyz = ph.get_txyz();
       plg, txyz(,3), txyz(,2);
       limits, square=1;

       
  SEE ALSO: gyoto.Metric, gyoto.Astrobj, gyoto.Photon, gyoto.Screen,
            gyoto.Scenery, gyoto.Spectrum, gyoto.Spectrometer
            utilities: gyotoy, gyoto_plg3, gyoto.debug,
            gyoto_plgsky, gyoto_plmksky, gyoto_pltsky, gyoto_reticle,
            gyoto_orient3, gyoto_convert
    
 */

//////// SCENERY

extern gyoto_Scenery;
/* DOCUMENT scenery = gyoto.Scenery([filename,] [members=values ...])
             Create GYOTO Scenery object
         or scenery, [members=values]
             Set GYOTO Scenery member
         or res = scenery(member=) or res = scenery.member()
             Get GYOTO Scenery member
         or scenery, xmlwrite=filename
             Save Scenery description to XML
         or data = scenery([ irange, jrange [, quantities ]])
             Ray-trace scenery

   PURPOSE:
    Create and use GYOTO Scenery objects.

    A GYOTO Scenery is used to render an image by relativistic
    ray-tracing. A GYOTO Scenery contains references to:
     - a GYOTO Metric ("the stage");
     - a GYOTO Astrobj ("the actors");
     - a GYOTO Screen ("the camera").

    See GYOTO for basic concepts and syntax.


   MEMBERS:

    metric=  see gyoto.Metric(): what "straight" means for light
             travel;
             
    screen=  see gyoto.Screen(), it specifies where the obseerver is
             located and the obseerving time;
             
    astrobj= see gyoto.Astrobj(): where the light comes from;
    
    delta=   a double scalar, the initial integration step for the
             Photons laucnched during ray-tracing;
             
    adaptive= whether integration uses adaptive step

    maxiter= maximum number of iterations in integration (per Photon)
             
    quantities= an array of strings giving the names of physical
             quantities that should be retrieved during
             ray-tracing. Beware that not all Astrobj kinds support
             all of those quantities. The order in which the
             quantities are listed when setting this member is not
             relevant, the following are equivalent:
                 data = scenery(quantities=["EmissionTime", "Intensity"])();
             and
                 data = scenery(quantities=["Intensity", "EmissionTime"])();
             Setting quantities here is not mandatory as the third
             positional argument used for ray-tracing permits to
             override it in an ordered fashion.

             Recognized quantities:
               "Intensity": apparent intensity of the Astrobj;
               "EmissionTime": time at which each photonreaching the
                    screen was emitted;
               "MinDistance": minimum distance ever reached between
                    each photon (whether comingfrom the object or not)
                    and the Astrobj);
               "FirstDistMin": First local minimum in the
                    Astrobj/Photon distance;
               "Redshift": ;
               "ImpactR", "ImpactX", "ImpactY" and "ImpactZ": R
                    (specrical), X, Y and Z (Cartsesian) coordinates
                    at which a photon was emitted by the Astrobj;
               "Spectrum": Inu spectrum of the Astrobj on this pixel,
                    the spectrometer is specified in the Screen
                    object;
               "BinSpectrum": spectrum of the Astrobj on this pixel
                    (the spectrometer is specified in the Screen
                    object), as would be detected by a real
                    spectrometer: the value in each spectral channel
                    is the integral of Inu over the spectral channel;
               "User1" to "User5": other specific scalar quantities an
                    Astrobj may be able to compute, refer to the
                    documentation for the Astrobj kind of your choice.

    nthreads=number of parallel threads to use in
             gyoto.Scenery_rayTrace. This has no effect when
             ray-tracing using the "data = scenery()" syntax below.


    METHODS:
    xmlwrite, help

    MPI MULTI-PROCESSING:

     You can check whether MPI support is compiled-in using
     gyoto.haveMPI(). To benefit from it, you must call
     gyoto.MPI_Init() once, before using the following keywords. Before
     quitting Yorick, should should also call 'sc, mpispawn=0;
     gyoto.MPI_Finalize;'. Failure to do that may lead to helper
     processes lingering around.

     nprocesses=number of parallel processes to use. If you just set
            nprocesses this way, they will be launched (but not
            stopped) automatically the first time ray-tracing is
            attempted. Setting nprocesses after processes have been
            launched does not change the number of processes running.
             
     mpispawn=number of parallel MPI jobs (a.k.a. workers) to
            spawn. This also sets nprocesses, and actually launches
            the processes. use mpilaunch=0 to stop the processes (they
            consume CPU cycles even when idle). Always use mpiclone to
            send the current Scenery to the workers before calling
            gyoto.Scenery_rayTrace(). This can be done in the same
            call as mpispawn: sc, mpispawn=12, mpiclone=;

     mpiclone=[anything]: clone this Scenery into the MPI workers (see
             mpispawn). The workers do not stay synchonised
             automatically. mpiclone must always be called right
             before ray-tracing (i.e. after any modification to the
             Scenery or any of the objects it contains, and after
             having called mpispawn).
                    
    RAY-TRACING:
    
     Ray-traced data is retrieved calling the object like a function
     with no keyword:
        data = scenery ();
     or data = scenery (ispec, jspec, quant )
     or data = scenery (alpha, delta, quant )

    ISPEC and JSPEC are about anything that can be use to index an
    array (except range functions, wildcards and rubber indices):
    scalar integers, ranges (imin:imax:step), or arrays of
    integers. The result is that
       scenery() (ispec, jspec, ) == scenery(ispec, jspec, );
    Indices are in the range 1 to scenery.screen().resolution().

    ALPHA and DELTA are double scalars or arrays of scalars. If they
    are both arrays, they must have the same dimensions. They are used
    as angles to address the screen. Angles are normally in the range
    [-0.5,0.5]*scenery.screen().fov().

    QUANT is an array of Yorick strings where each element selects one
    quantity to retrieve (see the QUANTITIES member above). QUANT may
    be void to use the quantities already set in the Scenery or the
    default for the Astrobj. If specifyng QUANT, DATA will be a MxNxP
    double array, where P=numberof(QUANT) and DATA(,,i) will contain
    the value of the quantity specified by QUANT(i).

    QUANTITIES may also be a scalar to retrieve a single
    quantity.

    The "Spectrum" quantity is a bit peculiar since it take more than
    one plane in data.

   PARALLEL COMPUTING:

    Gyoto supports parallel computing using either multi-threading
    (pthreads) or multi-processing (MPI). Not all classes work well
    using multi-threading (in particular, Lorene metrics are not
    thread-safe). On the other hand, all of Gyoto should support
    multi-processing.

    If you want to use mutli-processing, you should take care of:
      - calling 'sc, mpispawn=<nprocs>' before sending ray-tracing;
      - calling 'sc, mpiclone=;' once the scenery is ready, before tracing;
      - before quitting yorick:
          * destroy all of your sceneries with sc[];
          * call gyoto.MPI_Finalize().

    libgyoto may automatically initialize MPI, in which case it will
    also automatically terminate it. To prevent this behaviour, use
    mpispawn at least once before ray-tracing.
    
   SEE ALSO:
     gyoto.Metric, gyoto.Screen, gyoto.Astrobj, gyoto.Photon,
     gyoto.Spectrometer, gyoto.Scenery_rayTrace
*/

extern gyoto_Scenery_rayTrace
/* DOCUMENT res = gyoto.Scenery_rayTrace(scenery, imin, imax, jmin, jmax,
                                         impactcoords)

     If IMPACTCOORDS is an unadorned, nil variable it is output. If it
     is an expression or non-nil, it is input.

     This function is a very thin wrapper around the C++ method
     Gyoto::Scenery::rayTrace(). It behaves exactly the same as the
     command-line utility gyoto, except:
     
     - the data are returned as a Yorick array instead of saved to a
       FITS file (if multiple quantities are requested, they are
       stored in successive planes in he output data cube);
     - (i|j)(min|max) default to sane values.
       
     Gyoto provides an alternative manner for ray-tracing, which is
     simply calling the Scenery object without any keyword
     argument. The two methods are not equivalent.
     
     Example using MPI multi-processing:
       #include "gyoto.i"
       restore, gyoto;
       sc = Scenery("../doc/examples/example-polish-doughnut.xml");
       data = Scenery_rayTrace(sc(mpispawn=12, mpiclone=));
       sc, mpispawn=0;
 */

func _gyoto_Scenery_adaptive_raytrace(sco, respmax, &computed) {
/* xDOCUMENT data = gyoto_Scenery_adaptive_raytrace(scenery, pmax, [computed])

BROKEN
   
     Ray-trace a GYOTO Scenery on an adaptive grid.

     For certain kinds of objects (in particular, Stars), this routine
     is much faster than the equivalent:
       data = scenery(resolution=3^pmax, raytrace=1);
     It is NOT guaranteed that the two methods yield the same result.

     The minimum distance between photon and object is first computed
     on a coarse grid which is then refined as required.

   SEE ALSO: gyoto.Scenery
 */
  write, format="%s\n",
    "WARNING: gyoto_Scenery_adaptive_raytrace() is under development";

  sc = sco(clone=); // don't modify calling object
  
  DBL_MAX=1e100;

  screen = sc(screen=);

  respmax=long(respmax);
  
  resp=1;
  resmax=3^respmax;
  step=long(3^(respmax-resp));
  first=step/2+1;

  data=array(double, resmax, resmax, 6);

  quantities = ["Intensity", "EmissionTime", "MinDistance", "ImpactX", "ImpactY", "ImpactZ"];
  
  screen, resolution=resmax;
  data(first::step, first::step, ) =
    sc(first::step, first::step, quantities);
  
  computed=array(long,resmax,resmax);

  nb=9;
  for (resp=2, res=3; resp<=respmax; ++resp) {

    // largest distance to neighbour
    dsub=data(first::step, first::step,);
    ind=where(dsub>=DBL_MAX);
    if (numberof(ind)) dsub(ind)=DBL_MAX;
    delta = array(double, res, res);
    d1=(dsub(dif,,4:6)^2)(,,sum);
    d2=(dsub(,dif,4:6)^2)(,,sum);
    delta(2:-1,2:-1)=[d1(:-1,2:-1), d1(2:,2:-1), d2(2:-1,:-1), d2(2:-1,2:)](,,max);
    delta(1,2:-1)=[d1(1,2:-1), d2(1,:-1), d2(1,2:)](,max);
    delta(0,2:-1)=[d1(0,2:-1), d2(0,:-1), d2(0,2:)](,max);
    delta(2:-1,1)=[d2(2:-1,1), d1(:-1,1), d1(2:,1)](,max);
    delta(2:-1,0)=[d2(2:-1,0), d1(:-1,0), d1(2:,0)](,max);
    delta(0,0)=[d1(0,0), d2(0,0)](max);
    delta(1,1)=[d1(1,1), d2(1,1)](max);
    delta(0,1)=[d1(0,1), d2(0,1)](max);
    delta(1,0)=[d1(1,0), d2(1,0)](max);


    // ! BEWARE : res is updated here
    res*=3;
    refine=array(int, res, res);
    refine(1::3, 1::3) =
    refine(2::3, 1::3) =
    refine(3::3, 1::3) =
    refine(1::3, 2::3) =
    refine(2::3, 2::3) =
    refine(3::3, 2::3) =
    refine(1::3, 3::3) =
    refine(2::3, 3::3) =
    refine(3::3, 3::3) =
      (dsub(,,3)<4*delta) | (dsub(,,3)<2);
    
    nstep=long(3^(respmax-resp));
    nfirst=first-nstep;
    data(nfirst     ::step,nfirst     ::step,)=
    data(first      ::step,nfirst     ::step,)=
    data(first+nstep::step,nfirst     ::step,)=
    data(nfirst     ::step,first      ::step,)=
    data(first+nstep::step,first      ::step,)=
    data(nfirst     ::step,first+nstep::step,)=
    data(first      ::step,first+nstep::step,)=
    data(first+nstep::step,first+nstep::step,)=
      dsub;

    step=nstep;
    first=nfirst;
    for (i=1; i<=res; ++i) {
      ibis=(i-1)*step+first;
      ind=where(refine(i,) & !computed(ibis,first::step));
      nb+=numberof(ind);
      if (numberof(ind)) {
        indbis=(ind-1)*step+first;
        data(ibis,indbis,)=sc(ibis, indbis, quantities);
        computed(ibis, indbis)=resp;
      }
    }

  }
  nb;
  return data;
}

func gyoto_matte_paint(set, paint, kind=, yaw=, pitch=, roll=)
/* DOCUMENT bg = gyoto.matte_paint(set, paint[, coordkind=coordkind])

   Vizualize lensing effects on a painted background. Only geometrical
   distorsion effects are rendered.

   gyoto.matte_paint() will ray-trace the Scenery SET, compute the
   origin direction of each Photon that do originate from infinite,
   and pass these directions to PAINT.

   INPUTS:
    set: a Gyoto Scenery or impact coordinates as computed with
      sc(,,impactcoords=);
    paint: a function-like object that will perform the actual
      painting. PAINT will be called as:
        bg = paint(theta, phi, mask=mask);
      where THETA and PHI denote the direction from which each Photon
      originates, as seen from the center of the coordinate system,
      and mask is 0 for Photons that do not exist (e.g. photons that
      would fall into the black-hole.

   OUTPUT:
    The rendered image.

 */
{
  restore, gyoto;
  coordkind=kind
  if (is_Scenery(set)) {
    data=sc(,,impactcoords=)(9:,,);
    coordkind=sc.metric().coordkind();
  } else {
    if (is_void(coordkind))
      error, "KIND must be provided if SET is not a Scenery";
    if (dimsof(set)(2)==8) data=set;
    else data=set(9:,,);
  }


  // Prepare rotation matrix to transform from compact object
  // coordinates to painter coordinates. Yaw, pitch and roll are the
  // angles defining the orientation of the painter relative to the
  // metric coordinate system.
  // First, R is identity.
  R=array(0., 3, 3);
  for (i=1; i<=3; ++i) R(i, i)=1.;
  // Then, apply each rotation in order.
  // yaw, about Z axis
  if (!is_void(yaw)) R = gyoto.rotation(3, yaw);
  // pitch, about rotated Y axis
  if (!is_void(pitch)) R = gyoto.rotation(2, pitch)(,+)*R(+,);
  // roll, about rotated X axis
  if (!is_void(roll)) R = gyoto.rotation(1, roll)(,+)*R(+,);
  
  // Transform the photon's 4-velocity into Cartesian coordinates
  if (coordkind==gyoto.coordkind.spherical) {

    t=data(1,,);
    r=data(2,,);
    theta=data(3,,);  
    phi=data(4,,);
    tdot=data(5,,);
    rp=rdot=data(6,,);
    thp=thetadot=data(7,,);
    php=phidot=data(8,,);

    ind=where(tdot);
    taup=tdot; taup(ind)=1./tdot(ind);
    rp =rdot    *taup;
    php=phidot  *taup;
    thp=thetadot*taup;

    sth=sin(theta);
    cth=cos(theta);
    sph=sin(phi);
    cph=cos(phi);

    ur= [sth*cph, sth*sph, cth];
    uth=[cth*cph, cth*sph, -sth];
    uph=[   -sph,     cph,   0];

    durdph=sth*uph;
    durdth=uth;

    durdt=php*durdph+thp*durdth;

    v=rp*ur + r*durdt;

    mask0=char(r<1e200);
  
  } else {
    tdot=data(5,,);
    xdot=data(6,,);
    ydot=data(7,,);
    zdot=data(8,,);
    v=[xdot, ydot, zdot]/tdot(,,-);
    mask0=char(tdot<1e200);
  }

  // Express this in the painter base:
  V=R(,+)*v(,,+);

  // transform to spherical coordinates
  vproj2=V(1,,)^2+V(2,,)^2;

  vr=sqrt(vproj2+V(3,,)^2);
  vph=atan(V(2,,), V(1,,));
  vth=atan(sqrt(vproj2), V(3,,));

  theta=pi-vth;
  phi=vph-pi;
  if (numberof((ind=where(phi<-pi)))) phi(ind)+=2.*pi;

  bg=paint(theta, phi);

  if (dimsof(bg)(1) == 3) {bg*=mask0(-,,);}
  else bg *= mask0;
  
  return bg;
}

if (is_func(use))
  painters_junk=save(painters_junk,
                     gyoto_painters_p_mode_eval,
                     gyoto_painters_mk_p_mode,
                     gyoto_painters_panorama_eval,
                     gyoto_painters_mk_panorama,
                     gyoto_painters_picture_eval,
                     gyoto_painters_mk_picture);
  
func gyoto_painters_p_mode_eval(theta, phi, mask=) {
  use, ntheta, nphi;
  bg=sin(ntheta*theta)*sin(nphi*phi)+1.;
  if (!is_void(mask)) bg*=mask;
  return bg;
}

func gyoto_painters_mk_p_mode(ntheta=, nphi=)
/* DOCUMENT painter = gyoto.painters.mk_p_mode([ntheta=ntheta, [nphi=nphi])

    Make a painter for gyoto.matte_paint(). In this painter, the sky
    is mainted in a p-mode spherical harmonic fashion, with m=nphi and
    l=nphi+ntheta.
   
   SEE ALSO: gyoto.matte_paint
*/
  {
  if (is_void(ntheta)) ntheta=1;
  if (is_void(nphi)) nphi=1;
  return closure(save(ntheta=ntheta, nphi=nphi,
                      painter_eval=gyoto.painters.p_mode_eval),
                 painter_eval);
}

func gyoto_rotation(axis, angle)
/* DOCUMENT R=gyoto.rotation(axis, angle)

     Yield rotation matrix about axis AXIS (Ox=1, Oy=2, Oz=3) by
     amount ANGLE.

     If vector V holds the coordinates of a point in original the
     coordinate base,
       R(,+)*v(+)
     holds coordinates of same point in rotated base while
       R(+,)*v(+)
     holds coordinates of rotated vector in original base.

     In both cases, if S and T are two rotations, R=T(,+)*S(+,) if the
     rotation obtained by first applying S, then T.
     
   SEE ALSO:
 */
{

  ca=cos(angle);
  sa=sin(angle);

  R2=[[ca, -sa], [sa, ca]];

  R3=array(0., 3, 3);
  for (i=1; i<=3; ++i) R3(i, i)=1.;

  plane=(axis+[1, 2])%3;

  if (numberof((ind=where(plane==0)))) plane(ind)=3;

  R3(plane, plane)=R2;

  return R3;
}

func gyoto_painters_picture_eval(theta_spherical, phi_spherical, mask=) {
  // see http://mathworld.wolfram.com/GnomonicProjection.html
  use, img, phi1, lambda0;
  dmatte=dimsof(img);
  ndims=dmatte(1);
  nx=dmatte(-1);
  ny=dmatte(0);


  phi=pi/2-theta_spherical ; // latitude
  lambda=-phi_spherical ; // longitude

  cos_c=sin(phi1)*sin(phi)+cos(phi1)*cos(phi)*cos(lambda-lambda0);
  x=cos(phi)*sin(lambda-lambda0)/cos_c;
  y=(cos(phi1)*sin(phi)-sin(phi1)*cos(phi)*cos(lambda-lambda0))/cos_c;

  // two extreme points in the long dimension
  phi_corners=phi1+fov*0.5*[-1.,1.]*(ny>nx);
  lambda_corners=lambda0+fov*0.5*[-1.,1.]*(nx>=ny);
  cos_c_corners=sin(phi1)*sin(phi_corners)+
    cos(phi1)*cos(phi_corners)*cos(lambda_corners-lambda0);
  x_corners=cos(phi_corners)*sin(lambda_corners-lambda0)/
    cos_c_corners;
  y_corners=(cos(phi1)*sin(phi_corners)-
             sin(phi1)*cos(phi_corners)*cos(lambda_corners-lambda0))/
    cos_c_corners;

  scale=(nx>=ny)?(nx/x_corners(dif)):(ny/y_corners(dif));

  // Below, +0.5 so that for a 4x4 grid, the center is at 2.5, and
  // another +0.5 so that the long() conversion yields the closest
  // approximation rather than the bottom truncation.
  i0=nx*0.5 +1.;
  j0=ny*0.5 +1.;

  xp=i0 + x*scale ;
  yp=j0 + y*scale ;

  if (numberof((ind=where(xp<=0)))) xp(ind)=1;
  if (numberof((ind=where(yp<=0)))) yp(ind)=1;
  if (numberof((ind=where(xp>nx)))) xp(ind)=nx;
  if (numberof((ind=where(yp>ny)))) yp(ind)=ny;

  dd=dimsof(phi);

  if (ndims==3) bg=array(structof(img), dmatte(2), dd(2), dd(3));
  if (ndims==2) bg=array(structof(img), dd(2), dd(3));
  for (i=1; i<=dd(2); ++i) {
    for (j=1; j<=dd(3); ++j) {
      bg(..,i,j)=img(..,long(xp(i, j)),long(yp(i, j)));
    }
  }

  return bg;
}


func gyoto_painters_mk_picture(img=, fov=, phi1=, lambda0=)
/* DOCUMENT painter = gyoto.painters.mk_picture(img=img)

    Make a painter for gyoto.matte_paint(). In this painter, only a
    fraction of the celestial sphere is available in rectilinear
    (=gnomonic) projection. In other words, as a normal photographic
    picture. Missing pixels are painted like the closest available
    pixel.

   Keywords:
    Mandatory:
    img: the image, suitable for displaying with pli.

    Optional:
    h_fov: horizontal angle of view of the image, in radians.
    v_fov: vertical field of view of the image.
    fov:   field-of-view in the largest dimension of the image.

    By default, FOV corresponds to 36mm film behind a 50mm
    focal-length lense (about 40°) and the pixel scale is assumed to
    be the same in the horizontal and vertical dimensions. All angles
    must be specified in radians.

   SEE ALSO: gyoto.matte_paint, .painters.mk_p_mode
*/
{
  if (is_void(img)) error, "img can't be void";
  d=dimsof(img);
  nx=d(-1);
  ny=d(0);
  if (is_void(fov)) fov=2.*atan(36./100.);
  if (is_void(phi1)) phi1=0.;
  if (is_void(lambda0)) lambda0=0.;
  return closure(save(img, h_fov, v_fov, phi1, lambda0,
                      painter_eval=gyoto.painters.picture_eval),
                 painter_eval);
}

func gyoto_painters_panorama_eval(theta, phi, mask=) {
  use, img, phi_fov, theta_fov;
  dmatte=dimsof(img);
  ndims=dmatte(1);
  nphi=dmatte(ndims);
  ntheta=dmatte(ndims+1);

  if (is_void(phi_fov)) phi_fov=2.*pi;
  if (is_void(theta_fov)) theta_fov=pi;


  // That's in pixel per radian.  Minus signs because theta increases
  // from top to bottom and phi from right to left, as seen from
  // inside the sphere.
  phi_scale=-nphi/phi_fov;
  theta_scale=-ntheta/theta_fov;

  // Below, +0.5 so that for a 4x4 grid, the center is at 2.5, and
  // another +0.5 so that the long() conversion yields the closest
  // approximation rather than the bottom truncation.
  i0=nphi*0.5 +1.;
  j0=ntheta*0.5 +1.;
  
  x=long(i0 + phi*phi_scale) ;
  y=long(j0 + (theta-0.5*pi)*theta_scale) ;

  if (numberof((ind=where(x<=0)))) x(ind)=1;
  if (numberof((ind=where(y<=0)))) y(ind)=1;
  if (numberof((ind=where(x>nphi)))) x(ind)=nphi;
  if (numberof((ind=where(y>ntheta)))) y(ind)=ntheta;
  
  dd=dimsof(phi);
  
  if (ndims==3) bg=array(structof(img), dmatte(2), dd(2), dd(3));
  if (ndims==2) bg=array(structof(img), dd(2), dd(3));
  for (i=1; i<=dd(2); ++i) {
    for (j=1; j<=dd(3); ++j) {
      bg(..,i,j)=img(..,x(i, j),y(i, j));
    }
  }

  return bg;
}


func gyoto_painters_mk_panorama(img=, phi_fov=, theta_fov=)
/* DOCUMENT painter = gyoto.painters.mk_panorama(img=img)

    Make a painter for gyoto.matte_paint(). In this painter, the sky
    is painted using an image. The image should be a 360°x180°
    panorama rendered in equirectangular projection. If the image does
    not cover the full celestial sphere, use phi_fov and theta_fov to
    specify the field of view in each direction. img should be an
    array(char, 3, nphi, ntheta).
   
   SEE ALSO: gyoto.matte_paint, .painters.mk_p_mode
*/
{
  if (is_void(phi_fov)) phi_fov=2.*pi;
  if (is_void(theta_fov)) theta_fov=pi;
  return closure(save(img, phi_fov, theta_fov,
                      painter_eval=gyoto.painters.panorama_eval),
                 painter_eval);
}

local gyoto_painters;
/* DOCUMENT namespace gyoto.painters
   
    Contains painters to be used with gyoto.matte_paint.

    EXAMPLE:
     sc=gyoto.Scenery("file.xml");
     painter=gyoto.painters.mk_p_mode(ntheta=80, nphi=80);
     bg=gyoto.matte_paint(sc, painter);

   SEE ALSO: gyoto.matte_paint, .painters.mk_p_mode
 */

if (is_func(use)) {
  gyoto_painters=save(p_mode_eval=gyoto_painters_p_mode_eval,
                      mk_p_mode=gyoto_painters_mk_p_mode,
                      picture_eval=gyoto_painters_picture_eval,
                      mk_picture=gyoto_painters_mk_picture,
                      panorama_eval=gyoto_painters_panorama_eval,
                      mk_panorama=gyoto_painters_mk_panorama);
  restore, painters_junk;
 }

// PHOTON CLASS
extern gyoto_Photon;
/* DOCUMENT photon = gyoto.Photon([filename], [members=values])
            photon, member=values
            value = photon(member=) or value = photon.member
            value = photon(function_method=params)
            photon, subroutine_method=params

   PURPOSE:

     Instantiate and use a single GYOTO photon.

     Photons are mass-less particles following light-like geodesics of
     a metric. For basic concepts, see GYOTO. For ray-tracing,
     gyoto.Scenery() is more appropriate.

   MEMBERS:

     Members can be set with the syntax "photon, member=value" and
     retrieved with the syntax "value=photon(member=)" or
     "value=photon.member":
     
        metric= a GYOTO Metric (see gyoto.Metric),
            initcoord=scenery,x,y also sets the metric.
        astrobj=a GYOTO Astrobj (see gyoto.Astroj), the target of
            ray-tracing.
                  
        initcoord= the initial coordinates (4-position & 4 velocity).
            There are many ways to specify this:

            initcoord=COORD8
                directly give the 8 coordinates;
            initcoord=POS4,VEL4
                as above, with a coma in-between;
            initcoord=POS4,VEL3 this time only 3 coordinates are given
                for the velocity vector. The light-ray will be tangent
                to this 3-vector.
            initcoord=SCREEN,DALPHA,DDELTA
                SCREEN is a gyoto.Screen, DALPHA and DDELTA specify
                the direction this photon comes from when it reaches
                the screen. DALPHA and DDELTA are in radians and must
                be floating-point values.
            initcoord=SCREEN,I,J
                As above, but I and J are integers specifying the
                pixel of the arrival SCREEN which the photon hits.
            initcoord=SCENERY,DALPHA,DDELTA
            initcoord=SCENERY,I,J
                As above, but specify a gyoto.Scenery instead of a
                gyoto.Screen. The Metric and Astrobj of the Senery
                will also be attached to the Photon.

            Those last ways of specifying the initial conditions are
            very useful to get the trajectory of a specific photon in
            a ray-traced scenery.

        spectro= a gyoto.Spectrometer

        delta= integration step (initial in case of adaptive, the
               default)

        adaptive= whether integration uses adaptive step

        maxiter= maximum number of iterations in integration

        integrator= "Legacy" | "runge_kutta_fehlberg78" |
            "runge_kutta_cash_karp54" |"runge_kutta_dopri5" |
            "runge_kutta_cash_karp54_classic"
            The integrator to use.

        deltamin=, deltamax=, deltamaxoverr=, abstol, reltol:
            numerical tuning parameters for the integrators other than
            Legacy (the tuning parameters for the Legacy integrator
            are set in the Metric object). You should normally need to
            tweak only abstol and reltol. If you find you need to
            change the others, try a higher order integrator, such as
            runge_kutta_fehlberg78.
     
   SUBROUTINE-LIKE METHODS:

     Several of these keywords can by specified whenever creating or
     accessing the object.

     xfill=TLIM Integrate the geodesic from the time specified with
            INITCOORD to tlim (the integrated geodesic remains stored
            in the PHOTON);

     save_txyz=FILENAME Dump the integrated geodesic in cartesian
            coordinates in ASCII file FILENAME.

     xmlwrite=filename as usual, save an XML description of this
            photon;

     help=; list properties;

     setparameter=name[,value] set any low-level parameter. Name and
            value are strings.

   FUNCTION-LIKE METHODS:

     The object PHOTON will return a value when called as a function
     with the following keywords set:

     is_hit=     Return 1 if this photon hits the Astrobj

     get_txyz=   Return the geodesic in Cartesian coordinates:
                     data = photon(xfill=tlim, get_txyz=)
            data will be a Nx4 double array where data(i,) contains
            the 4-position in Cartesian coordinates of the photon for
            all the dates computed by the integrator between
            INITCOORD[0] and TLIM.

     get_coord= Return the geodesic in Metric coordinatess: same as
            above, but in the prefered coordinate system for this
            metric, which may be Cartesian or spherical.

     get_coord=dates Same as above, but for the dates specified in
            double array DATES.

     get_cartesian=dates Get the 3-position and 3-velocity of the
            Photon in Cartesian coordinates for the specified dates.
     
   SEE ALSO: gyoto, gyoto.Metric, gyoto.Screen, gyoto.Scenery,
            gyoto.Astrobj
 */

/// METRIC

extern gyoto_Metric;
/* DOCUMENT gg = gyoto.Metric( filename, [members=values] )
            gg, members=values
            retval = gg(member=) or retval = gg.member;
            retval = gg(function_method=par1, par2...)
            gg, xmlwrite=filename
            coef = gg(coordinates, mu, nu)
            gg, help=;

   PURPOSE:
     Create and manipulate GYOTO Metric objects

   INTRODUCTION:

     See GYOTO for basic concepts.
   
     The GYOTO plug-in for Yorick introduces "Metric" objects (see
     GYOTO for an introduction). Such objects are created used for
     instance the gyoto.KerrBL() function. Any kind of metric (even if
     not explicitely exposed in this plug-in) can be loaded from an
     XML file using the FILENAME parameter. This XML file can by any
     GYOTO file containing a Metric section: the top-level can be a
     Metric, a Scenery, an Astrobj or whatever which contains a
     reference to a Metric.

     When printed, a Metric displays an XML description of itself,
     which can be dumped to a file using the XMLWRITE keyword.

     The HELP method prints the list of properties supported by this
     astrobj.

     Most GYOTO functions which accept a Metric as a parameter accept
     any kind of Metric. There are nevertheless specific
     functionsAstrobjs which make sense only in the framework of a
     specific kind of Metric, notably gyoto.KerrBL. This is the case
     for gyoto.PolishDoughnut for instance.


   MEMBER KEYWORDS:

     mass=, unitlength=, kind=, deltamin=, deltamax=, keplerian=
   
     All the Metric kinds have a mass that can be set and retrieved
     using the mass keyword:
        gg, mass=value;
        value = gg(mass=);
        
     Setting the mass gives the scale of the black hole in physical
     units. The unit length can be retrieve using the unitlength
     keyword:
        len = gg(unitlength=)
     where LEN is in meters if MASS was set in kilograms.

     The kind of the metric (e.g. "KerrBL") can be queried:
        kind_string = gg(kind=)

     The coordinate kind of this metric can be compared with members
     of the oxy group gyoto.coordkind:
        if (gg.coordkind() == gyoto.coordkind.cartesian) ...;

     The minimum and maximum value of the adaptive step in
     the integrator can be set using respectively deltamin and
     deltamax.

     KEPLERIAN can be set to 1 in oreder to get CIRCULARVELOCITY to
     return the Keplerian approximation.

   METHODS

     Without any keywords, the metric can output its coefficient at
     4-position POSITION:
        coefs = gg(position, mu, nu);
     where mu and nu are indices or ranges (1-4, 1 is for time).
        coefs = gg(position)
     returns the 16 metric coefficients at position.
   
     Additional function-like or subroutine like methods:
     
       coord=gg(prime2tdot=pos, vel): COORD is the 8-vector where
               COORD[1-4]==POS and COORD[5-8] is the 4-velocity
               corresponding to the 3-velocity VEL;
       
       vel4 =gg(nullifycoord=pos, vel)
               Return nullified (photon) coord tangent to vel at pos.

       prod = gg(scalarprod=pos, u1, u2)
               Return (4-dimensional) scalar product.

       vels = gg(circularvelocity=coords [, dir])
               On input, COORDS is an array of doubles yielding
               4-position vectors: x0=coords(1,..), x1=coords(2, ..),
               x3=coords(3, ..), x4=coords(4, ..). The return value
               VELS has the same dimensionality as coords(1:4, ..). It
               contains the circular velocity (on te equatorial plane,
               so COORDS is projected) corresponding to each position
               specified by coords.

               If the optional parameter DIR is -1, VELS corresponds
               to the counter-rotating circular velocity.

       Gamma= gg(christoffel=pos)
               Return the 64 Christoffel symbols at posistion POS.
               Beware of the index ordering:
                  Gamma(nu, mu, alpha) == Gamma^alpha_mu_nu

       gg, help=; print help on properties.

   SEE ALSO: gyoto, gyoto.KerrBL, gyoto.KerrKS
 */

extern gyoto_Astrobj;
/* DOCUMENT ao = gyoto.Astrobj( filename );
            ao, member1=val1, member2=val2...;
            val = ao(member=) or val = ao.member
            ao, xmlwrite=filename
            ao, help=;

     Generic class for describing an astronomical object (a star, an
     accretion disk...) in the context of GYOTO. See "help, gyoto" for
     an introduction.

     Several specific kinds of objects can be created using the
     functions listed in "SEE ALSO". Any kind of Astrobj can be
     loaded from an XML file using gyoto.Astrobj. This XML file can be
     any GYOTO file containing an Astrobj section (for instance, a
     Scenery file).
     
     When printed, an astrobj displays an XML description of
     itself. This description can be saved to disk using the XMLWRITE
     method keyword.

     The HELP method prints the list of properties supported by this
     astrobj.

   MEMBER KEYWORDS
   
     All the kinds of Astrobjs share a few members that can be
     accessed with the following keywords. To set the member, use "ao,
     member=value". To retrieve the current value of the member, use
     "retval = ao(member=)" or "retval = ao.member".

     rmax:        for optimization, tell the Scenery that this object
                  does not extend over RMAX from the center of the
                  coordinate system (RMAX is expressed in geometrical
                  units).
                  
     opticallythin: whether or not to compute radiative transfer
                  through the object for optically thin objects.

     kind:        read only. ao(kind=) return a string containing the
                  kind of this object (e.g. "Star", "PolishDoughnut").

   METHOD KEYWORDS

     Astrobj objects can also provide function-like or subroutine-like
     methods. All of them support the following (subroutine-like):
   
     xmlwite="filename.xml" dump a description of the object to an XML
                  file.

     setparameter="name","content" generic method to set a parameter
                  in an Astrobj object, even if is has not been
                  explicitely exposed in the Yorick plug-in. For
                  instance, if st is a gyoto.Star object, the two
                  following commands yield the same result, although
                  the former is faster than the latter:
                     st, radius=1.0;
                     st, setparameter="Radius","1.0";
                  The list of values that "name" can take is (or
                  should be) documented in the doxygen documentation
                  for the specific class. See
                  Gyoto::Astrobj::Generic::setParameter().
     
   EXAMPLES:
    The following implement specific objects, most require gyoto_std.i:
     gyoto.Star               A spherical object moving along a geodesic
     gyoto.FixedStar          A spherical object of constant coordinates
     gyoto.Torus              A simple torus (solid, Keplerian rotation)
     gyoto.ThinDisk           A geometrically thin disk
     gyoto.PatternDisk        As above, emission numerically provided
     gyoto.Disk3D             Thick disk, emission numerically provided

    It is also possible to instantiate an astrobj by kind name even if
    this kind is not explicitly implemented in the yorick plug-in:
     gyoto.Astrobj("PageThorneDisk")
                              A geometrically thin, optically thick disk
                              with Page & Thorne 1974 emission

   SEE ALSO: gyoto
     
 */

extern gyoto_ThinDisk;
/* DOCUMENT ao = gyoto.ThinDisk( filename );
            ao, member1=val1, member2=val2...;
            val = ao(member=) or val = ao.member
            ao, xmlwrite=filename
            ao, help=;

     A more specific version of the gyoto.Astrobj function. A very
     crude Astrobj can be instantiated using gyoto.ThinDisk. More
     elaborate derived classes also exist. gyoto.ThinDisk accepts a
     few keywords in addition to those processed by gyoto.Astrobj.
            
   MEMBER KEYWORDS

     innerradius:  inner radius of the disk.
                  
     outer radius: outer radius of the disk.

     dir:          1 if corotating (relative to the coordinate system),
                  -1 if coounter rotating.

   SEE ALSO: gyoto, gyoto.Astrobj
    There is one derived classe in gyoto_std.i:
     gyoto.PatternDisk        As above, emission numerically provided
    try also gyoto.Astrobj("PageThorneDisk")
     
 */

extern gyoto_Screen;
/* DOCUMENT screen = gyoto.Screen([keyword=value ...])
             
         or gyoto.Screen, screen, [keyword=value]
         or screen, [keyword=value]
             Set GYOTO Screen member
             
         or res = gyoto.Screen(screen, get_keyword=1)
         or res = screen(get_keyword=1)
             Get GYOTO Screen member

   PURPOSE:
    Create and use GYOTO Screen objects.

    A GYOTO Screen is in essence the virtual camera used for
    ray-tracing a GYOTO Scenery. It contains informations such as the
    resolution of the output images and the location of the camera
    relative to the Metric.
     
    A Screen is usually created using the gyoto.Screen function. It's
    members are then set or queried using various member keywords (see
    KEYWORDS below, "help, gyoto"). A description of the Screen can be
    written to an XML file using the xmlwrite= keyword. The help=
    keyword prints the list of poperties that the Screen class
    supports.

   INPUT:

    SCREEN: an opaque object referencing a GYOTO Screen instance. In
         the first form, SCENERY is instantiated with any member set
         accordingly to the KEYWORD=VALUE pairs. In the four other
         forms, SCREEN must have been instantiated by a previous call
         to gyoto.Screen().

   OUTPUT:
    SCREEN: in the first form, a new GYOTO Screen object is returned
         (actually an opaque Yorick object referencing a GYOTO
         object).

    RES: in the last two forms, the Screen object is queried thanks to
         a keyword starting with "get_" and the call returns the
         corresponding value.

   KEYWORDS:

     Member keywords are used to set the value of a member (screen,
     member=value) or to retrieve the value of a member
     (value=screen(member=) or value=screen.member). Many support
     "unit".:
       metric, time, tmin, fov (field-of-view), resolution (N pixels
       on each side), dangle1, dangle2, distance (meters), inclination,
       paln (position angle of the line of nodes), argument, freqobs,
       pojection (=[incl, paln, arg]), observerpos (alternative way to
       set time, dist, incl, and arg by giving the position of the
       camera. Here, dist is in geometrical units), fourvel, screen1,
       screen2, screen3, spectro (see gyoto.Spectrometer);

    Screens provide two function-like methods and the usual
    subroutine-like method xmlwrite plus two methods for specifying
    a mask:

     skycoord = [x0, x1, x2, x3] returns screen (=sky) coordinates
       d_alpha, d_delta, d_dist

     raycoord = [d_alpha, d_delta] returns 8-coordinate vector
       (position-velocity) of photon hitting screen with arrival
       direction d_alpha, d_delta.

     xmlwrite = filename to export screen to XML.

     help=nil list properties.

     mask= [nil | filename | image] use nil to retrieve the mask,
       filename to load mask from FITS file, image array to set mask
       from Yorick array.

     maskwrite=filename: save mask to FITS file.

   SEE ALSO:
     gyoto, gyoto.Metric, gyoto.Scenery, gyoto.Astrobj,
     gyoto.Spectrometer
*/


extern gyoto_Spectrum;
/* DOCUMENT sp = gyoto.Spectrum([filename, ][members=values]
         or sp, method=parameters
         or value = sp(method=parameters)
         
         or Inu = sp(frequency)
         or I_nu1_nu2 = sp(integrate=channels)

   A gyoto object with the usual gyoto syntax representing the
   spectrum of a source. So far, it can be attached only to a Star.

   The specificity of a spectrum is that it can be used as a function
   to retrieve the intensity of an object at a given wavelength.

   METHODS:
     operator():
       Without any keyword, SP is used as a function to evaluate Inu
       at a given FREQUENCY:
         Inu = sp(frequency)
       FREQUENCY may be an array if any shape, Inu will be
       conformable. FREQUENCY is in Hz.

     integrate=channels
       CHANNELS must be a 1-dimensional array of N doubles (N>=2). The
       returned value will be an array of doubles with N-1 elements:
         I_nu1_nu2 = sp(integrate=channels)
       I_nu1_nu2(i) is the integral of Inu from frequency nu1 to
       frequency nu2. This is the same than the BinSpectrum quantity
       in a Scenery.
     
   GENERIC MEMBERS:
     kind=   (read-only)

   GENERIC SUBROUTINE-LIKE METHODS
     xmlwrite=filename (as usual)
     help=nil (as usual)

   GENERIC FUNCTION-LIKE METHODS
     clone=  (as usual)
     
   SEE ALSO:  gyoto, gyoto.Star
 */

extern gyoto_debug;
/* DOCUMENT gyoto.debug, 1/0
    Turn GYOTO debug output on/off.
   SEE ALSO: gyoto
 */

extern gyoto_verbose;
/* DOCUMENT gyoto.verbose, level
         or level = gyoto.verbose()
    Set/get Gyoto verbosity level
   SEE ALSO: gyoto
 */

extern gyoto_Spectrometer;
/* DOCUMENT spectro = gyoto.Spectrometer([filename],[members=values])
         or spectro, xmlwrite=filename
         or var = spectro( channels= | midpoints= | widths= |
                           nsamples= | kind= | clone= )
         
     The spectrometric capabilities of a GYOTO Screen or Photon.

     The spectrometer attached to a Screen in a Scenery determines the
     frequencies for which the Quantities Spectrum and BinSpectrum are
     computed when ray-tracing.

     For basics, see GYOTO.

     In the first form, if FILENAME is specified, a new Spectrometer
     is instantiated from XML file FILENAME. FILENAME can also name a
     Spectrometer Kind, such as "Complex" or "wave".

   KEYWORDS:
       unit=     a string, affects methods called together:
                 channels_in_microns=spectro(channels=, unit="microm")

   METHODS:
       xmlwrite=filename, write XML desription of this spectrometer
                  to hard drive. E.g. spectro, wmlwrite="file.xml"

       help=nil, list properties
   
     The following keywords allow retrieving information on the
     spectrometer. They don't take a value but return one with this
     syntax:  retval=spectro(keyword=)

       kind=      returns the Spectrometer kind name, e.g. "wave",
                  "freq" or "Complex".
       nsamples=  returns the number of spectral channels.
       channels=  return an array(double, NSAMPLES, 2) yielding the
                  edges of each spectral channel in UNIT (default: Hz)
       midpoints= returns an array(double, NSAMPLES) yielding the central
                  frequency of each spectral channel in UNIT (default: Hz).
       widths=    returns the width of each spectral channel in UNIT
                  (default: Hz).
       clone=     returns a deep copy of this Spectrometer.
       
   SEE ALSO: gyoto, gyoto.Screen, gyoto.Scenery,
             gyoto.SpectroUniform, gyoto.SpetroComplex
 */

extern gyoto_SpectroUniform;
extern _gyoto_SpectroUniform_register_as_Spectro;
/* DOCUMENT spectro = gyoto.SpectroUniform([members=values])
         or spectro, keywords=values
         or retval = spectro( keyword= )
         
     gyoto.SpectroUniform() is a more specific implementation of
     gyoto.Spectrometer(). It allows creating an Spectrometer::Uniform
     object which models a spectrometer where NSAMPLES spectral
     channels are unifromly spaced in wavelength, frequency, or log10
     thereof.

     gyoto.SpectroUniform has a superset of the functionalities of
     gyoto.Spectrometer.

   MEMBERS:
     Members can be set with "spectro, member=value" and retrieved
     with "value = spectro(member=)" or "value=spectro.member".
     
       kind=     a string, one of "none", "wave", "freq", "wavelog",
                 "freqlog". KIND affects how BAND below is interpreted
                 and how CHANNELS and MIDPOINTS are returned. Since
                 these four "kinds" are actually based on the same
                 code, it is possible to convert a Spetrcometer
                 between them.

       nsamples= a long, the number of spectral channels in this
                 spectrometer;

       band=     a pair of doubles, the lower and upper boundaries
                 of the spectral band covered by this spectrometer.
                 The default unit in which BAND is expressed depends
                 on KIND:
                       KIND             UNIT
                       freq              Hz
                      freqlog         log10(Hz)
                       wave              m
                      wavelog         log10(m)
                 Keyword UNIT can be set to use a different unit than
                 the default. For kinds freqlog and wavelog,
                 log10(unit) is used.

   METHODS:
     The following methods from gyoto_Spectrometer are implemented:
     xmlwrite, help, clone, channels, midpoints, widths.
       
   SEE ALSO: gyoto.Spectrometer, gyoto.SpectroComplex
 */
_gyoto_SpectroUniform_register_as_Spectro;

extern gyoto_SpectroComplex;
extern _gyoto_SpCplx_register_as_Spectrometer;
/* DOCUMENT spectro = gyoto.SpectroComplex([members=values])
         or spectro, keywords=values
         or retval = spectro( keyword= )
         or subspectro = spectro(index)
         
     gyoto.SpectroComplex() is a more specific implementation of
     gyoto.Spectrometer(). It allows creating an Spectrometer::Complex
     object which models a spectrometer made of several simpler
     Spectrometer objects (sub-spectrometers).

     gyoto.SpectroComplex has a superset of the functionalities of
     gyoto.Spectrometer.

     A single sub-spectrometer can be retrieved by its INDEX. INDEX is
     1-based, as is customary in Yorick.

   EXAMPLE:
     sp1 = gyoto.SpectroUniform(kind="wave",
                                nsamples=10,
                                unit="microm",
                                band=[1, 2]);
     sp2 = gyoto.SpectroUniform(kind="freqlog",
                                nsamples=10,
                                unit="eV",
                                band=[10, 20]);
     sp = gyoto.SpectroComplex(append=sp1)(append=sp2);
     sp1bis = sp(1);
     sp2bis = sp(2);
     
   METHODS:
     The following methods from gyoto.Spectrometer are implemented:
     xmlwrite, help, clone, channels, midpoints, widths. In addition,
     gyoto.SpectroComplex accepts the following methods:
       append=subspectro add new sub-spectrometer
       remove=index      remove sub-spectrometer INDEX
       
   SEE ALSO: gyoto.Spectrometer, gyoto.SpectroUniform
 */
_gyoto_SpCplx_register_as_Spectrometer;

func gyoto_plg3(x,y,z, keywords=) {
/* DOCUMENT gyoto_plg3, x, y, z

    Plot (plg) a graph in 3D using the API provided by pl3d.i
     
   SEE ALSO: pl3d.i, set3_object, gyoto_plgsky
 */
  if (_draw3) {
    xyz= _nxt(x);
    keywords= _nxt(x);
    get3_xy, xyz, x, y;
    plgk, y, x, keywords=keywords;
    return;
  }
  set3_object, gyoto_plg3, _lst(transpose([x,y,z]), keywords);
}

func gyoto_plmk3(x,y,z, keywords=) {
/* DOCUMENT gyoto_plmk3, x, y, z

    Plot (plg) a graph in 3D using the API provided by pl3d.i
     
   SEE ALSO: pl3d.i, set3_object, gyoto_plgsky
 */
  if (_draw3) {
    xyz= _nxt(x);
    keywords= _nxt(x);
    get3_xy, xyz, x, y;
    plmkk, y, x, keywords=keywords;
    return;
  }
  set3_object, gyoto_plmk3, _lst(transpose([x,y,z]), keywords);
}

func gyoto_plgsky(alpha,delta, keywords=) {
/* DOCUMENT gyoto_plgsky, alpha, delta

    Plot (plg) data on the plane of the sky using the pl3d.i API. Not
    all transforamations are allowed (basically orient3 and rot3).

   SEE ALSO: gyoto_pltsky, gyoto_plmksky, gyoto_reticle, gyoto_plg3
 */
  if (_draw3) {
    list=alpha;
    alpha= _nxt(list);
    delta= _nxt(list);
    keywords= _nxt(list);
    scale= _getscl3()(1);
    plgk, delta*scale, -alpha*scale, keywords=keywords;
    return;
  }
  set3_object, gyoto_plgsky, _lst(alpha, delta, keywords);
}

func gyoto_plmksky(alpha,delta, keywords=) {
/* DOCUMENT gyoto_plmksky, alpha, delta

    Plot (plmk) data on the plane of the sky using the pl3d.i API. Not
    all transforamations are allowed (basically orient3 and rot3).

   SEE ALSO: gyoto_pltsky, gyoto_plgsky, gyoto_reticle
 */
  if (_draw3) {
    list=alpha;
    alpha= _nxt(list);
    delta= _nxt(list);
    keywords= _nxt(list);
    scale= _getscl3()(1);
    plmkk, delta*scale, -alpha*scale, keywords=keywords;
    return;
  }
  set3_object, gyoto_plmksky, _lst(alpha,delta, keywords);
}

func gyoto_pltsky(text, alpha, delta, keywords=) {
/* DOCUMENT gyoto_plgsky, text, alpha, delta, justify

    Write text (plt) on the plane of the sky using the pl3d.i API. Not
    all transforamations are allowed (basically orient3 and rot3).

   SEE ALSO: gyoto_plgsky, gyoto_plmksky, gyoto_reticle
 */
  if (_draw3) {
    list=text;
    text   = _nxt(list);
    alpha  = _nxt(list);
    delta  = _nxt(list);
    keywords= _nxt(list);
    scale= _getscl3()(1);
    if (is_void(keywords)) keywords=GraphK();
    keywords.tosys=&long(1);
    pltk, text, -alpha*scale, delta*scale, keywords=keywords;
    return;
  }
  set3_object, gyoto_pltsky, _lst(text, alpha,delta, keywords);
}

func gyoto_reticle(alphamax, deltamax, alphamin, deltamin, set=) {
/* DOCUMENT gyoto_reticle, alphamax [, deltamax, alphamin, deltamin]

    Draw alpha and delta axes on the plane of the sky using the pl3d.i
    API. Not all transforamations are allowed (basically orient3 and
    rot3).

   SEE ALSO: gyoto_plgsky, gyoto_plmksky, gyoto_pltsky
 */
  if (is_void(deltamax)) deltamax = alphamax;
  if (is_void(alphamin)) alphamin = -alphamax;
  if (is_void(deltamin)) deltamin = -deltamax;
  gyoto_plgsky,[alphamin, alphamax], [0,0];
  gyoto_plgsky,[0,0],[deltamin, deltamax];

  // Compute tick length
  len=max(alphamax-alphamin,deltamax-deltamin);

  _3ticks, alphamin, alphamax, _3nmajor;
  for (i=1; i<=numberof(xminor); ++i)
    gyoto_plgsky, xminor(i)(-:1:2), [-len*_3lminor, 0];
  for (i=1; i<=numberof(xmajor); ++i) {
    gyoto_plgsky, xmajor(i)(-:1:2), [-len*_3lmajor, 0];
    gyoto_pltsky, xlabel(i), xmajor(i), -len*_3lmajor,
      keywords=GraphK(justify=&string("CT"));
  }

  _3ticks, deltamin, deltamax, _3nmajor;
  for (i=1; i<=numberof(xminor); ++i)
    gyoto_plgsky, [len*_3lminor, 0], xminor(i)(-:1:2);
  for (i=1; i<=numberof(xmajor); ++i) {
    gyoto_plgsky, [len*_3lmajor, 0], xmajor(i)(-:1:2);
    gyoto_pltsky, xlabel(i), len*_3lmajor, xmajor(i),
      keywords=GraphK(justify=&string("RH"));
  }
}

func gyoto_orient3(incl, paln, phase) {
/* DOCUMENT gyoto_orient3, incl, paln, phase

     Call orient3 and rot3 after having converted the 3 (more or less)
     standard astronomical angles (expressed in degrees) into the
     corresponding ones in the pl3d.i API.

   SEE ALSO: orient3, rot3
 */

  local deg2rad;
  deg2rad=pi/180.;
  //  orient3,phase*deg2rad,incl*deg2rad-pi/2.;
  //  rot3,,,paln*deg2rad+pi/2.;
  //  restore3, [];

  // Reset orientation to something normal
  orient3, 0, 0;
  rot3, , pi, ;
  rot3, -pi/2,,;
  // start from here
  rot3,              , ,-phase*deg2rad;
  rot3, -incl*deg2rad, ,              ;
  rot3,              , ,-paln *deg2rad;
}

func gyoto_convert(&x, mass, distance, unit) {
/* DOCUMENT gyoto_convert, x, mass, distance, unit
         or x_unit = gyoto_convert(x, mass, distance, unit)

   INPUT:
    X: value in geometrical units
    MASS: mass defining geometrical units in solar masses
    DISTANCE: distance to observer in kiloparsecs
    UNIT: unit to convert to. One of "geometric", "m", "km", "sun
          radius", "arcsec", "mas" (milliarcseconds), "uas"
          (microarcseconds)

   OUTPUT:
    X (in the first form) or X_UNIT: the converted value.
          
   SEE ALSO:
 */

  m_sun = 1.98843e30;     // kg
  G     = 6.67428e-11;    // Gravitational constant, SI = m^3 * kg^-1 * s-2
  c     = 299792458.;     // m/s
  r_sun = 6.955e8;        // m
  kpc   = 3.08568025e19;  // m
  
  GMsunoC2 = G*m_sun/c^2;
  
  unit_d = GMsunoC2 * mass; // m

  if (unit=="geometrical")     x_unit=x;
  else if (unit=="m")          x_unit=x* unit_d;
  else if (unit=="km")         x_unit=x*(unit_d*1e-3);
  else if (unit=="sun radius") x_unit=x*(unit_d/r_sun);
  else if (unit=="rad")        x_unit=x*(unit_d/(distance*kpc));
  else if (unit=="degree")     x_unit=x*(unit_d/(distance*kpc*pi)*180);
  else if (unit=="arcmin")     x_unit=x*(unit_d/(distance*kpc*pi)*1.08e4);
  else if (unit=="arcsec")     x_unit=x*(unit_d/(distance*kpc*pi)*6.48e5);
  else if (unit=="mas")        x_unit=x*(unit_d/(distance*kpc*pi)*6.48e8);
  else if (unit=="uas")        x_unit=x*(unit_d/(distance*kpc*pi)*6.48e11);
  /* need double checking
  else if (unit=="s")          x_unit=x*(unit_d/c);
  else if (unit=="min")        x_unit=x*(unit_d/(c*60.));
  else if (unit=="hour")       x_unit=x*(unit_d/(c*3600.));
  else if (unit=="day")        x_unit=x*(unit_d/(c*86400.));
  else if (unit=="month")      x_unit=x*(unit_d/(c*2629743.8));
  else if (unit=="year")       x_unit=x*(unit_d/(c*31556926));
  */
  else error, "unknown unit "+unit;
  if (am_subroutine()) eq_nocopy, x, x_unit;
  else return x_unit;
}

extern is_gyoto_Photon;
extern is_gyoto_Astrobj;
extern is_gyoto_Metric;
extern is_gyoto_Spectrometer;
extern is_gyoto_Spectrum;
extern is_gyoto_Screen;
extern is_gyoto_Scenery;
/* DOCUMENT bool = gyoto.is_BASE(arg)
   
     BOOL is 1 if arg is a gyoto.BASE where base is Metric, Astrobj,
     Spectrum, Spectrometer, Scenery, Screen...

   SEE ALSO: gyoto
 */



extern gyoto_dontcatchSIGFPE;
extern gyoto_dontcatchSIGSEGV;
/* DOCUMENT gyoto.dontcatchSIGFPE
        and gyoto.dontcatchSIGSEGV

     By default, Yorick will catch SIGFPE and SIGSEGV and throw a
     Yorick error when they occur. Although this is arguably the right
     thing to do in the general case, this makes it difficult to debug
     compiled code that triggers these signals. These functions
     disables catching the corresponding signal, so that a debugger
     such as GDB can be used to track down their origin.

   SEE ALSO: gyoto.fedisableexcept
 */

extern gyoto_fedisableexcept;
extern gyoto_feenableexcept;
extern gyoto_haveFENV;
extern gyoto_FE;
/* DOCUMENT gyoto.fedisableexcept[(excepts)];
            gyoto.fedisableexcept[(excepts)];
            gyoto.haveFENV();
            except=gyoto.FE(except_name);

     Yorick tries very hard to set the floating-point environment such
     that e.g. divisions by zero will trigger SIGFPE. Gyoto tries to
     not allow SIGFPE to occur and this normally is a bug when it
     does. However, it is sometimes impossible to avoid it: depending
     on the compiler used to build Gyoto and its dependencies,
     optimizations may be done under the assumption that the code will
     always run with all exceptions disabled. See BUGS.md in the
     Gyoto source code for details.

     Gyoto therefore provides the following functions to control the
     floating-point environment from Yorick scripts. Parameters are
     discribed below.

      gyoto.fedisableexcept([int excepts]);
        Disable exceptions specified in EXCEPT. By default, all
        exceptions normally enabled by Yorick (FE_DIVBYZERO |
        FE_OVERFLOW | FE_INVALID). Values suitable for EXCEPT are
        those returned by gyoto.FE() and any ored combination of
        those.

      gyoto.feenableexcept([int excepts]);
        Enable exceptions specified in EXCEPT. By default, all
        exceptions normally enabled by Yorick.

      gyoto.haveFENV();
        Return 1 if support for setting the floating-point environment
        is built-in, else 0.

      gyoto.FE(string except_name);
        Get integer flag representation of a given exception, e.g.
        FE_DIVBYZERO=gyoto.FE("DIVBYZERO").  EXCEPT_NAME is one of:
        DIVBYZERO, INVALID, INEXACT, OVERFLOW, UNDERFLOW, ALL_EXCEPT.

     EXAMPLES:
      Disable all exceptions normally enabled by Yorick:
       gyoto.fedisableexcept();

      Disable absolutely all exceptions:
       gyoto.dedisableexcept(gyoto.FE("ALL_EXCEPT"));

      Specifically enable overflow and underflow exceptions, but not
      DIVBYZERO:
       gyoto.feenableexcept(gyoto.FE("OVERFLOW") | gyoto.FE("UNDEFLOW"));

      Enable all exceptions normally enabled by Yorick:
       gyoto.feenableexcept();

      Enable absolutely all exceptions:
       gyoto.feenableexcept(gyoto.FE("ALL_EXCEPT"))

   SEE ALSO: gyoto.dontcatchSIGFPE
 */

extern gyoto_listRegister;
/* DOCUMENT gyoto.listRegister
     List the register of known Astrobj, Metric, Spectrum and
     Spectrometer kinds.
   SEE ALSO: gyoto
*/

if (is_func(use)) {
  include, "gyoto_namespace.i", 1;
 } else {
  gyoto_namespace=noop;
 }

if (!GYOTO_NO_STD) include, "gyoto_std.i", 2;
