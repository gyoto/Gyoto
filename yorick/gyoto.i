plug_in, "gyoto";
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

#include "graphk.i"

extern gyoto_haveXerces;
/* DOCUMENT have_xerces = gyoto_haveXerces()
    Tell whether GYOTO was compiled with Xerces support (XML i/o)
   OUTPUT:
    HAVE_XERCES=1 if compiled with Xerces, else 0.
*/

extern __gyoto_setErrorHandler;
/* xDOCUMENT __gyoto_setErrorHandler
   Must be called once to attach the GYOTO error handler to Yorick's one
*/
__gyoto_setErrorHandler;

extern gyoto_loadPlugin;
/* DOCUMENT gyoto_loadPlugin, plugin[, plugin2[, plugin3]] [, nofail=1]

   Load Gyoto plug-ins.

  INPUTS:
   gyoto_loadPlugins() accepts an aribtrary number of positional
   arguments, each a string or string array naming individual Gyoto
   plugins. For instance, all of the following attempt to load the
   plug-ins, stdplug, lorene and myplug:
    gyoto_loadPlugin, "stdplug", "lorene", myplug"
    gyoto_loadPlugin, ["stdplug", "lorene", myplug"]
    gyoto_loadPlugin, "stdplug", ["lorene", myplug"]

  KEYWORDS:
   nofail= if set and true, failure to load a plug-in will not trigger
           an error. It applies to _all_ plug-ins in the list.

  EXAMPLE:
   gyoto_loadPlugin, "stdplug"
   gyoto_loadPlugin, "lorene", nofail=1
 */

extern __gyoto_initRegister;
/* xDOCUMENT __gyoto_initRegister
   Must be called once to initialize the GYOTO plug-in register
*/
__gyoto_initRegister;

extern __gyoto_exportSupplier;
/* xDOCUMENT __gyoto_exportSupplier
   To be called by client plug-ins like gyoto_std.i
*/


require, "pl3d.i";
#include "gyoto_Photon.i"
#include "gyoto_Scenery.i"
#include "gyoto_constants.i"

include, "gyoto_std.i", 2;


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


    Creating a GYOTO object
    -----------------------

    To create a GYOTO object, one calls one of the object creators:
    for instance
       gg = gyoto_KerrBL() ;
    results in GG being an opaque object containing a GYOTO Metric (a
    Kerr Metric in this case, using the Boyer-Lindquist coordinate
    system).

    The creators accept a filename as positional argument (see FILE
    FORMAT below): gg = gyoto_KerrBL( "my_KerrBL_descritpion.xml" );


    GYOTO OBJECTS BEHAVE LIKE FUNCTIONS
    -----------------------------------

    Once a GYOTO object has been created, it can be accessed in an
    object-like fashion (the object then ressembles a Yorick function)
    or through the same function which was used to create it. For
    instance, to change the spin parameter of the GG Metric above, the
    two following are exactly synonymous:
       gg, spin=0.8 ;
       gyoto_KerrBL, gg, spin=0.8;

    Beware: in the second form, the positional parameter GG must be
    given before any keyword parameter. This is one of the few points
    were the GYOTO syntax does not fully adheres the Yorick uses.

    
    MEMBERS
    -------
    
    All the GYOTO objects have properties called "members" such as the
    spin of a Kerr metric. Members can be set or retrieved using the
    correspinding keywords, either upon object instanciation or at any
    later point in time:
       gg = gyoto_KerrBL( spin=0.995, mass=4e6 ) ;
       gg, spin=0.5;
       spin_param = gg(spin=);
       
    Note how giving a member keyword without any value, as in the last
    example above, allows _retrieving_ the previously set value. When
    setting member, it is also possible to call the object as a
    function. In that case, the return value will be the object
    itself, allowing to call it again as a function:
       spin_val = gg(spin=0.5)(spin=);
    Although the above example is trivial, this is useful in some
    complex situations:
       Setting the resolution in a Screen attached to a Scenery:
         res_val = scenery(screen=)(resolution=);
       Setting resolution:
         noop, scenery(screen=)(resolution=res_val);
       (The noop above is not mandatory.) 

    Some member keywords accept more than one parameter, separated by
    comas. This is the first exception GYOTO makes to the Yorick
    syntax. For example, for setting the initial position and velocity
    of a star, one can use:
       st = gyoto_Star( initcoord=pos,vel );
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
       gg = gyoto_KerrBL();
       st = gyoto_Star( metric=gg );
       ph = gyoto_Photon( metric=gg );
       gg, spin=0.5;
    In the above, the spin change affect the Photon and the Star as
    well.

    If what you want is a detached copy of a Gyoto object, you can get
    it with the CLONE method. Following the above example:
    
       gg_copy = gg;   // This is copying: gg2 and gg are the same
                       // object
       
       gg_clone = gg(clone=);  // This is cloning: gg3 is a detached
                               // copy

       gg_copy, spin=0.2;
       gg_clone, spin=0.7;

       gg(spin=);  // the spin in gg is the same as in gg_copy, not
                   // gg_clone.

    
    METHODS
    -------

    Some keywords are "subroutine-like" methods in that they can have
    side-effects, but do not return a specific value. One notable
    example is the XMLWRITE keyword which most objects accept for
    dumping a description of themselves to an XML file:
       gg, xmlwrite="filename.xml";

    Other methods keywords are function-like: they usually take one or
    several parameter and return a value. Only one value-retruning
    keyword can be set at a time (be it a function-like method or
    member keyword set for retrieving a value, such as "spin="):
       coor = gg(makecoord=yinit, cst);
       txyz = star(get_cartesian=dates);

    Some objects have a default methods that is called when no keyword
    is present. This is the case for metrics, which return the metric
    coefficients at a specific location, and of sceneries, which
    perform ray-tracing:
       coefficient = gg(coordinates, mu, nu);
       image = scenery();
       
       
    GYOTO OBJECT TYPES
    ==================

    gyoto_Metric: the general relativity metric in which objects move...
          Notable sub-classes: gyoto_KerrBL, gyoto_KerrKS
    
    gyoto_Photon: massless particules used for ray-tracing
    
    gyoto_Astrobj: astrophysical objects.
          Notable sub-classes: gyoto_Star, gyoto_FixedStar,
          gyoto_PolishDoughnut, gyoto_ThinInfiniteDisk

    gyoto_Screen: the camera for ray-tracing

    gyoto_Spectrometer: the spectral capabilites of a gyoto_Screen
          (also found in gyoto_Photon)

    gyoto_Spectrum: a spectrum, only used in gyoto_Star so far
    
    gyoto_Scenery: relationship between all of the above

    
    FILE FORMAT
    ===========

    Most objects can be described in XML format. That's what you see
    when you print an object. You can also read/write an object
    from/to and XML file:
       sc = gyoto_Scenery( "description.xml" );
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
       data = gyoto_Star(radius=0.5,
               metric=(gg=gyoto_KerrBL(spin=0.995)),
               initcoord=[0, 10.791, 1.5708, 0], [0, 0, 0.0166637],
               xfill=800
               )(
                 get_skypos=gyoto_Screen(metric=gg)
                 );
       plg, data(,2), data(,1);

    Ray-trace a scenery (orbiting star near a Kerr black-hole):
       sc = gyoto_Scenery(
             metric=(gg=gyoto_KerrBL()),
             screen=gyoto_Screen(metric=gg,
                                 observerpos=[1000, 100., 0.15, 0.],
                                 fov=pi/10.,
                                 resolution=128),
             astrobj=gyoto_Star(metric=gg,
                                 radius=0.5,
                                 initcoord=[600, 9, pi/2, 0], [0, 0, 0.037037])
            )
       pli, sc(,,"Intensity");

     Trace the trajectory of a photon in the secondary image of the above:
       ph = gyoto_Photon(initcoord=sc,77,45, xfill=870.623);
       txyz = ph(get_txyz=1);
       plg, txyz(,3), txyz(,2);
       limits, square=1;

       
  SEE ALSO: gyoto_Metric, gyoto_Astrobj, gyoto_Photon, gyoto_Screen,
            gyoto_Scenery, gyoto_Spectrum, gyoto_Spectrometer
            utilities: gyotoy, gyoto_plg3, gyoto_debug,
            gyoto_plgsky, gyoto_plmksky, gyoto_pltsky, gyoto_reticle,
            gyoto_orient3, gyoto_convert
    
 */

extern gyoto_Metric;
/* DOCUMENT gg = gyoto_Metric( filename, [members=values] )
            gg, members=values
            retval = gg(member=);
            retval = gg(function_method=par1, par2...)
            gg, xmlwrite=filename
            coef = gg(coordinates, mu, nu)

   PURPOSE:
     Create and manipulate GYOTO Metric objects

   INTRODUCTION:

     See GYOTO for basic concepts.
   
     The GYOTO plug-in for Yorick introduces "Metric" objects (see
     GYOTO for an introduction). Such objects are created used for
     instance the gyoto_KerrBL() function. Any kind of metric (even if
     not explicitely exposed in this plug-in) can be loaded from an
     XML file using the FILENAME parameter. This XML file can by any
     GYOTO file containing a Metric section: the top-level can be a
     Metric, a Scenery, an Astrobj or whatever which contains a
     reference to a Metric.

     When printed, a Metric displays an XML description of itself,
     which can be dumped to a file using the XMLWRITE keyword.

     Most GYOTO functions which accept a Metric as a parameter accept
     any kind of Metric. There are nevertheless specific
     functionsAstrobjs which make sense only in the framework of a
     specific kind of Metric, notably gyoto_KerrBL. This is the case
     for gyoto_PolishDoughnut for instance.


   MEMBER KEYWORDS:

     mass=, unitlength=, kind=
   
     All the Metric kinds have a mass that can be set and retrieved
     using the mass keyword:
        gg, mass=value;
        value = gg(mass=);
        
     Setting the mass gives the scale of the black hole in physical
     units. The unit length can be retrieve using the unitlength
     keyword:
        len = gg(unitlength=)
     where LEN is in meters if MASS was set in kilograms.

     Finally, the kind of the metric (e.g. "KerrBL") can be queried:
        kind_string = gg(kind=)

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
       
       gg, nullifycoord=pos, vel    return nullified (photon) coord tangent
                                    to vel at pos.

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

   SET KEYWORDS:
     List of set-like keywords ("set" is never specified). Specific
     Metric kinds may recognize more:
       mass=new_mass                gyoto_Metric_setMass, gg, new_mass

   SEE ALSO: gyoto, gyoto_KerrBL, gyoto_KerrKS
 */

extern gyoto_Astrobj;
/* DOCUMENT ao = gyoto_Astrobj( filename );
            ao, member1=val1, member2=val2...;
            val = ao(member=)
            ao, xmlwrite=filename

     Generic class for describing an astronomical object (a star, an
     accretion disk...) in the context of GYOTO. See "help, gyoto" for
     an introduction.

     Several specific kinds of objects can be created using the
     functions listed in "SEE ALSO". Any kind of Astrobj can be
     loaded from an XML file using gyoto_Astrobj. This XML file can be
     any GYOTO file containing an Astrobj section (for instance, a
     Scenery file).
     
     When printed, an astrobj displays an XML description of
     itself. This description can be saved to disk using the XMLWRITE
     method keyword.

   MEMBER KEYWORDS
   
     All the kinds of Astrobjs share a few members that can be
     accessed with the following keywords. To set the member, use "ao,
     member=value". To retrieve the current value of the member, use
     "retval = ao(member=)".

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
                  instance, if st is a gyoto_Star object, the two
                  following commands yield the same result, although
                  the former is faster than the latter:
                     st, radius=1.0;
                     st, setparameter="Radius","1.0";
                  The list of values that "name" can take is (or
                  should be) documented in the doxygen documentation
                  for the specific class. See
                  Gyoto::Astrobj::Generic::setParameter().
     
   SEE ALSO: gyoto
    The following implement specific objects, most require gyoto_std.i:
     gyoto_Star               A spherical object moving along a geodesic
     gyoto_FixedStar          A spherical object of constant coordinates
     gyoto_Torus              A simple torus (solid, Keplerian rotation)
     gyoto_ThinDisk           A geometrically thin disk
     gyoto_PageThorneDisk     As above with Page & Thorne 1974 emission
     gyoto_PatternDisk        As above, emission numerically provided
     gyoto_Disk3D             Thick disk, emission numerically provided
     
 */

extern gyoto_ThinDisk;
/* DOCUMENT ao = gyoto_ThinDisk( filename );
            ao, member1=val1, member2=val2...;
            val = ao(member=)
            ao, xmlwrite=filename

     A more specific version of the gyoto_Astrobj function. A very
     crude Astrobj can be instanciated using gyoto_ThinDisk. More
     elaborate derived classes also exist. gyoto_ThinDisk accepts a
     few keywords in addition to those processed by gyoto_Astrobj.
            
   MEMBER KEYWORDS

     innerradius:  inner radius of the disk.
                  
     outer radius: outer radius of the disk.

     dir:          1 if corotating (relative to the coordinate system),
                  -1 if coounter rotating.

   SEE ALSO: gyoto, gyoto_Astrobj
    There are two derived classes in gyoto_std.i:
     gyoto_PageThorneDisk     As above with Page & Thorne 1974 emission
     gyoto_PatternDisk        As above, emission numerically provided
     
 */

extern gyoto_Screen;
/* DOCUMENT screen = gyoto_Screen([keyword=value ...])
             
         or gyoto_Screen, screen, [keyword=value]
         or screen, [keyword=value]
             Set GYOTO Screen member
             
         or res = gyoto_Screen(screen, get_keyword=1)
         or res = screen(get_keyword=1)
             Get GYOTO Screen member

   PURPOSE:
    Create and use GYOTO Screen objects.

    A GYOTO Screen is in essence the virtual camera used for
    ray-tracing a GYOTO Scenery. It contains informations such as the
    resolution of the output images and the location of the camera
    relative to the Metric.
     
    A Screen is usually created using the gyoto_Screen function. It's
    members are then set or queried using various member keywords (see
    KEYWORDS below, "help, gyoto"). A description of the Screen can be
    written to an XML file using the xmlwrite= keyword.

   INPUT:

    SCREEN: an opaque object referencing a GYOTO Screen instance. In
         the first form, SCENERY is instanciated with any member set
         accordingly to the KEYWORD=VALUE pairs. In the four other
         forms, SCREEN must have been instanciated by a previous call
         to gyoto_Screen().

   OUTPUT:
    SCREEN: in the first form, a new GYOTO Screen object is returned
         (actually an opaque Yorick object referencing a GYOTO
         object).

    RES: in the last two forms, the Screen object is queried thanks to
         a keyword starting with "get_" and the call returns the
         corresponding value.

   KEYWORDS:

     Member keywords are used to set the value of a member (screen,
     member=value) or to retrieve the value of a memeber
     (value=screen(member=)):
       metric, time, tmin, fov (field-of-view), resolution (N pixels
       on each side), distance (meters), inclination, paln (position
       angle of the line of nodes), argument, pojection (=[incl, paln,
       arg]), observerpos (alternative way to set time, dist, incl,
       and arg by giving the position of the camera. Here, dist is in
       geometrical units), spectro (see gyoto_Spectrometer);

    Screens provide two function-like methods and the usual
    subroutine-like method xmlwrite:

     skycoord = [x0, x1, x2, x3] returns screen (=sky) coordinates
       d_alpha, d_delta, d_dist

     raycoord = [d_alpha, d_delta] returns 8-coordinate vector
       (position-velocity) of photon hitting screen with arrival
       direction d_alpha, d_delta.

     xmlwrite = filename to export screen to XML.

   SEE ALSO:
     gyoto, gyoto_Metric, gyoto_Scenery, gyoto_Astrobj,
     gyoto_Spectrometer
*/


extern gyoto_Spectrum;
/* DOCUMENT sp = gyoto_Spectrum([filename, ][members=values]
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

   GENERIC FUNCTION-LIKE METHODS
     clone=  (as usual)
     
   SEE ALSO:  gyoto, gyoto_Star
 */

extern is_gyoto_Spectrum;

/* DOCUMENT bool = is_gyoto_Spectrum(arg)
     BOOL is 1 if arg is a gyoto_Spectrum
   SEE ALSO: gyoto, gyoto_Spectrum
 */

extern gyoto_debug;
/* DOCUMENT gyoto_debug, 1/0
    Turn GYOTO debug output on/off.
   SEE ALSO: gyoto
 */

extern gyoto_verbose;
/* DOCUMENT gyoto_verbose, level
         or level = gyoto_verbose()
    Set/get Gyoto verbosity level
   SEE ALSO: gyoto
 */

extern gyoto_Spectrometer;
extern gyoto_SpectroUniform;
/* DOCUMENT spectro = gyoto_Spectrometer([members=values])
         or spectro, members=values
         or table = spectro( channels= | midpoints= | widths= )
         
     The spectrometric capabilities of a GYOTO Screen or Photon.

     The spectrometer attached to a Screen in a Scenery determines the
     frequencies for which the Quantities Spectrum and BinSpectrum are
     computed when ray-tracing.

     For basics, see GYOTO.

   MEMBERS:
     Members can be set with "spectro, member=value" and retrieved
     with "value = spectro(member=)".
     
       kind=     a string, one of "none", "wave", "freq", "wavelog",
                 "freqlog". KIND affects how BAND below is interpreted
                 and how CHANNELS and MIDPOINTS are returned.

       nsamples= a long, the number of spectral channels in this
                 spectrometer;

       band=     a pair of doubles, the lower and upper boundaries
                 of the spectral band covered by this spectrometer.
                 The unit in which BAND is expressed depends on KIND:
                       KIND             UNIT
                       freq              Hz
                      freqlog          Log(Hz)
                       wave              m
                      wavelog          Log(m)

   METHODS:
     The following keywords allow retrieving the spectral band for
     each of the spectral channel in one of two forms:
       channels=  an array(double, NSAMPLES+1) yielding the edges
                  of each spectral channel, in the same unit as
                  BAND depending on KIND.
       midpoints= an array(double, NSAMPLES) yielding the central
                  _frequency_ of each spectral channel expressed in
                  the same unit as BAND depending on KIND.
       widths=    the width of each spectral channel in Hz.
       
       clone=     get a clone of this Spectrometer
       
   SEE ALSO: gyoto, gyoto_Screen, gyoto_Scenery
 */
extern _gyoto_SpectroUniform_register_as_Spectro;
_gyoto_SpectroUniform_register_as_Spectro;

extern gyoto_SpectroComplex;
extern _gyoto_SpCplx_register_as_Spectrometer;
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

func gyoto_warning(msg) {
  if (_gyoto_running) gyotoy_warning,msg;
  else error, msg;
}


extern is_gyoto_Astrobj;
/* DOCUMENT ret = is_gyoto_Astrobj( obj )

   Return 1 if OBJ is a GYOTO Astrobj

   SEE ALSO: gyoto gyoto_Astrobj
*/


extern gyoto_dontcatchSIGFPE;
extern gyoto_dontcatchSIGSEGV;

extern gyoto_listRegister;
/* DOCUMENT gyoto_listRegister
     list the register of known Astrobj, Metric and Spectrum kinds.
*/
