require, "gyoto.i";
gyoto_requirePlugin, "stdplug";
plug_in, "gyoto_std";
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

extern __set_GyotoStdPlugSupplier;
/* xxDOCUMENT __set_GyotoStdPlugSupplier
   Retrieve the functions from gyoto.so.
   Using a plug-in from another plu-gin is not easy.
*/
__set_GyotoStdPlugSupplier, __gyoto_exportSupplier();

extern _gyoto_KerrBL_register_as_Metric;
_gyoto_KerrBL_register_as_Metric;
extern _gyoto_KerrKS_register_as_Metric;
_gyoto_KerrKS_register_as_Metric;
/* xDOCUMENT _gyoto_KerrBL_register_as_Metric
    should be called once so that gg(get_spin=1...) works.
 */

extern gyoto_KerrBL;
extern gyoto_KerrKS;
/* DOCUMENT gg = gyoto_KerrBL([filename,][members=values])
         or gg = gyoto_KerrKS([filename,][members=values])
         
         or gg, members=values
         or value = gg(member=)

         or coord = gg(makecoord=yinit, cst) // KerrBL only

         or coefs = gg(position, mu, nu)

   PURPOSE

     Instanciate and use a GYOTO KerrBL (for Kerr Boyer-Lindquist) or
     Kerr (for Kerr Kerr-Schild) Metric.

     For basics, first read "help, gyoto" and "help, gyoto_Metric".

   SPECIFIC MEMBER
   
     In addition to the generic Metric members, KerrBL and KerrKS
     metrics have a spin that can be set and retrieved as follows:
       gg, spin=value;
       value = gg(spin=)

     Furthermore, KerrBL has two numerical tuning parameters: difftol=
     and deltamaxoverr=. Low values yield more accurate results at the
     expanse of computing time. Investigating the MinDistance map of a
     FixedStar helps finding the right hyper parameter set for
     investingating an object at a given location and a given
     resolution.
       
   SPECIFIC METHOD

     In addition to the generic methods provided by gyoto_Metric,
     KerrBL metrics provide the following:
     
       coord=gg(get_coord=yinit, cst)
                get 8-coordinate (4-position & 4-velocity) COORD
                corresponding to the 6 coordinate (4-position & 2
                momenta) YINIT and the 3 motion constants in CST.

   SEE ALSO: gyoto, gyoto_Metric
 */
gyoto_namespace, KerrKS=gyoto_KerrKS;
gyoto_namespace, KerrBL=gyoto_KerrBL;

/////////// STAR
extern _gyoto_Star_register_as_Astrobj;
_gyoto_Star_register_as_Astrobj
/* xDOCUMENT _gyoto_Star_register_as_Astrobj
      To be called exactly once ins gyoto_Star.i
*/
   
extern gyoto_Star;
/* DOCUMENT st = gyoto_Star([filename, ][members=values])

         or st, [members=values]
         or value = st(member=)

         or st, xfill=tlim;
         or data = st(function_method=parameters)
         
   PURPOSE:
   
     Create and manipulate GYOTO Star objects. Stars are a specific
     kind of Gyoto Astrobj which move in the metric following
     time-like geodesics and are sort of spherical (although what
     "spherical" means is coordinate-system-dependent and purely
     geometrical).

     There are two uses for GYOTO Stars:
        - for computing the trajectory of a massive particle in the metric;
        - as the light source for ray-tracing.

     For basic concepts, see gyoto and gyoto_Astrobj.
       

   MEMBERS:

     In addition to the generic members described in gyoto_Astrobj,
     Stars have the following members (which can be set with "st,
     member=value" and retrieved with "value=st(member=)"):
        radius= the radius of the sperical star, in geometrical units;
        metric= the metric of which the Stars follows the geodesics;
        initcoord=[x0, x1, x2, x3], [dx1/dx0, dx2/dx0, dx3/dx0]
     or initcoord=[x0, x1, x2, x3, dx1/dx0, dx2/dx0, dx3/dx0]
     or initcoord=[x0, x1, x2, x3, x0dot, x1dot, x2dot, x3dot]
              initial position of velocity, requires metric to have
              been set previously;
        delta=  integration step (initial if adaptive).
        adaptive= whether the integration uses adaptive step.
        maxiter= maximum number of iterations in integration.
        deltamaxoverradius, deltamaxoverdistance: numerical parameters

        integrator= "Legacy" | "runge_kutta_fehlberg78" |
            "runge_kutta_cash_karp54" |"runge_kutta_dopri5" |
            "runge_kutta_cash_karp54_classic"
            The integrator to use.

        deltamin=, deltamax=, deltamaxoverr=, abstol, reltol:
            numerical tuning parameters for the integrators other than Legacy
            (the tuning parameters for the Legacy integrator are set in the
            Metric object). You should normally need to
            tweak only abstol and reltol. If you find you need to
            change the others, try a higher order integrator, such as
            runge_kutta_fehlberg78.


   METHODS

     In addition to the generic Astrobj methods, stars provide the
     following:

     Subroutine-like (st, subroutine_keyword=parameters):
        reset=anything    (e.g.:   st, reset=;  )
              Forget the result of a previous integration (see XFILL
              below).
     
        xfill=tlim
              The orbit is computed between x0 (in INITCOORD) and
              tlim. XFILL requires INITIALCOOORD to have been
              set. It is not guaranteed that any of the computed dates
              is actually tlim, merely that the orbit is computed a
              least from T0 to TLIM.

              The result of the integretation is cached inside the
              star object and can be retrieved with the various get_*
              methods below. Beware that changing a parameter in the
              metric will ne be taken into account in the retrieved
              results. To start again the integration, do one of the
              following:
                   st, metric=new_metric;
                   st, initcoord=conditions;
                   st, reset=;
                   
        setparameter=name,string_value

     Function-like (retval=st(function_method=parameters)):

        startrace=tmin, tmax : make a StarTrace object. See
              gyoto_StarTrace.
     
      The following return the position of the star for all the dates
      that where evaluated by the integrated when XFILL, above, was
      called, in various coordinate systems:
        get_skypos=screen
              Retrieve RA, Dec and depth
                data = st(get_skypos=screen);
                da=data(,1), dd=data(,2), dD=data(,3);
              SCREEN is a gyoto_Screen;
              
              da, dd and dD are the offsets from the center of the
              Metric in RA, Dec and along the line-of-sight.
              
        get_txyz=anything
              Retrieve Cartesian coordinates
                data = st(get_txyz=);
                t=data(,1), x=data(,2), y=data(,3), z=data(,4);
              ANYTHING can be anything, including nil.

        get_coord=
              Retrieve 8-coordinates (in Metric-prefered coord system)
                 data = st(get_coord=);
                 data(i, )=[t_i, x1_i, x2_i, x3_i,
                            x0dot_i, x1dot_i, x2dot_i, x3dot_i]
                            
        get_prime=
              Retrieve 3-velocity coordinates (in Metric-prefered
              coord system)
                 data = st(get_prime=);
                 data(i, )=[x1dot_i/x0dot_i, x2dot_i/x0dot_i, x3dot_i/x0dot_i] 
              
      The following retrieve coordinates for specific dates. There is
      some interpolation involved which is done rather
      precisely. DATES is a scalar double or an array of doubles of
      any shape. All the output sub-arrays DATA(..,i) are conformable
      with DATES:

        get_coord=dates    Retrieve 7 coordinates:
                data = star (get_coord=dates);
              STAR is a Star and DATES is an array of double with N
              elements. DATA will be an array of Nx7 elements where
              each of the 7 subarrays DATA(,i) is conformable with
              DATES. In the following, i is the index in the usual
              physical notation (X0 = time):
                X0 = DATES
                for i in {1, 2, 3}, Xi = DATA(..,i);
                for i in {0, 1, 2, 3} Xidot = DATA(..,i+4).

        get_cartesian=dates    Retrieve 6-coordinates:
                data = star (get_cartesian=dates);
              STAR is a Star and DATES is an array of double with N
              elements. DATA will be an array of Nx6 elements where
              each of the 7 subarrays DATA(,i) is conformable with
              DATES.
                X = DATA(..,1); Y = DATA(..,2); Z = DATA(..,3);
                dX/dt = DATA(..,4); dY/dt = DATA(..,5); dZ/dt = DATA(..,6);

   EXAMPLE:
       data = gyoto_Star(get_txyz=1, initialpos=pos, v, xfill=tlim,
                         metric=gyoto_KerrBL(spin=a));
                
   SEE ALSO: gyoto_Metric, gyoto_Kerr, gyoto_StarTrace
 */
gyoto_namespace, Star=gyoto_Star;

/////////// STARTRACE
extern _gyoto_StarTrace_register_as_Astrobj;
_gyoto_StarTrace_register_as_Astrobj
/* xDOCUMENT _gyoto_StarTrace_register_as_Astrobj
      To be called exactly once ins gyoto_std.i
*/
   
extern gyoto_StarTrace;
/* DOCUMENT st = gyoto_StarTrace([filename, ][members=values])

     A StarTrace is like a Star that is at all time ubiquitous on its
     orbit. This allows efficiently precomputing a mask for
     ray-tracing many images of the star on its orbit. A StarTrace
     accepts all the same keywords as a Star, except the STARTRACE
     keyword itself. StarTrace accepts the STAR keyword instead to
     retrieve the underlying Star object.

   EXAMPLE:

     sc = gyoto_Scenery("file.xml");
     star = sc.astrobj;
     startrace = star(startrace=0, 100);
     startrace, adaptive=0, opticallythin=0, delta=1;
     sc, astrobj=startrace;
     mask = sc(,,"Intensity");
     noop, sc.screen(mask=mask);
     sc, astrobj=star;
     image = sc(,,);

   SEE ALSO: gyoto, gyoto_Star
   
 */
gyoto_namespace, StarTrace=gyoto_StarTrace;

///////// FIXEDSTAR

extern _gyoto_FixedStar_register_as_Astrobj;
_gyoto_FixedStar_register_as_Astrobj;

extern gyoto_FixedStar;
/* DOCUMENT fs = gyoto_FixedStar([filename, ][members=values])
            fs, members=values;
            retval = fs(member=)
            fs, xmlwrite=filename;

    A GYOTO astrophysicla object: a non-moving, coordinate-spherical
    blob.

    See GYOTO for basic concepts. This star is fixed in the coordinate
    system with 3D position POS3D, for any time. The time coordinate
    is missing in POS3D. RADIUS is the stars radius.

   
   MEMBER KEYWORDS:
   
    Members can be set with "fs, member=val" or retrieved with
    "relval=fs(member=)" where member is one of:
   
      metric=    any GYOTO Metric object
      position=  3D position of the star
      radius=    radius of the star, in geometrical units
      deltamaxoverradius, deltamaxoverdistance: numerical parameters
     
    In addition, the standard gyoto_Astrobj members and methods are
    supported.
   
   SEE ALSO: gyoto, gyoto_Astrobj, gyoto_Metric
   
 */

gyoto_namespace, FixedStar=gyoto_FixedStar;

//// TORUS
extern _gyoto_Torus_register_as_Astrobj;
_gyoto_Torus_register_as_Astrobj;

extern gyoto_Torus;
/* DOCUMENT torus = gyoto_Torus([filename, ][members=values]);
            torus, members=values;
         or value = torus(member=);
         or torus, xmlwrite=filename

     A simple Torus for use with GYOTO. See "help, gyoto" for basics.

     The Torus is in solid rotation at the circular velocity along the
     circle of diameter LARGERADIUS. SMALLRADIUS is the radius of a
     vertical cross-section. Torus emission is set using SPECTRUM and
     OPACITY.
     
     MEMBERS:
       largeradius, smallradius, spectrum, opacity

    SEE ALSO: gyoto, gyoto_Astrobj, gyoto_Star, gyoto_Spectrum
*/
gyoto_namespace, Torus=gyoto_Torus;

//// PATTERNDISK
extern _gyoto_PatternDisk_register_as_Astrobj;
_gyoto_PatternDisk_register_as_Astrobj;
extern gyoto_PatternDisk;
/* DOCUMENT disk = gyotoPatternDisk(...)
            disk, member=value...

    This is a subkind of gyoto_ThinDisk. The disk is "painted" with a
    pattern, hence the name. The grid for the pattern is set by three
    numbers: NPHI (number of grid points in the azimuthal direction),
    REPEATPHI if the pattern must be repeated several times in the
    azimuthal direction (i.e. if the angular periodicity of the
    pattern is a fraction of 2*pi), and NR (number of grid points in
    the radial direction). The disk extends from INNERRADIUS to
    OUTERRADIUS (see gyoto_ThinDisk) with a regular spacing along the
    radial direction, unless RADIUS is specified.

    The pattern is specified by the surface brightness EMISSION==Jnu
    at NNU frequencies going from NU0 to NU0*DNU*(NNU-1). The cube
    EMISSION is an array(double, NNU, NPHI, NR).

    An optional OPACITY cube with the same dimensions as EMISSION can
    be provided. This allows using PatternDisk for any solid object
    confined in the equatorial plane and orbiting the central object
    in circular motion.

    By default, the fluid is supposed to be corotating at the local
    circular velocity, but the fluid velocity field can be specified
    with VELOCITY==array(double, 2, NPHI, NR).
    VELOCITY(1,..)==dphi/dt; VELOCITY(2,..)==dr/dt.

    The fluid VELOCITY field must not be mistaken by the apparent
    pattern velocity. The pattern is is solid (apparent) rotation at
    angular velocity PATTERNVELOCITY.

   KEYWORDS:

    fitsread="filename.fits"  read pattern from FITS file.
    
    fitswrite="filename.fits" write pattern to FITS file.

    patternvelocity=double(value) set (or get) pattern angular
                       velocity.

    repeatphi=N the pattern angular periodicity is 2*pi/N.

    nu0=        first frequency (Hz)

    dnu=        frequencty spacing (Hz)

    copyintensity=EMISSION
                * if EMISSION is nil, retrieve the surface brightness
                  cube;
                * if EMISSION==0, free the cube;
                * if EMISSION is an array of NNU x NPHI x NR doubles,
                  attach (copy) this array into DISK as the surface
                  brightness cube. If this cube doesn't have the same
                  dimensions as the previously set one, the velocity
                  and radius arrays will also be freed (as they have
                  inconsistent dimensions).

     copyopacity=OPACITY
                same as above for the opacity cube.

     copyvelocity=VELOCITY
                same as COPYINTENSITY but to attach the fluid velocity
                field, a 2 x NPHI x NR array where
                VELOCITY(1,..)==dphi/dt and VELOCITY(2,..)==dr/dt.

     copygridradius=RADIUS
                same as above, but RADIUS is a NR element vector
                specifying the R coordinate of the grid points. If
                RADIUS is not attached (if set, it can be detached
                with copygridradius=0), the grid points are regularly
                spaced between INNERRADIUS and OUTERRADIUS (see
                gyoto_ThinDisk).

   SEE ALSO:
    gyoto_Astrobj, gyoto_ThinDisk
*/
gyoto_namespace, PatternDisk=gyoto_PatternDisk;

//// DIRECTIONAL DISK
extern _gyoto_DirectionalDisk_register_as_Astrobj;
_gyoto_DirectionalDisk_register_as_Astrobj;
extern gyoto_DirectionalDisk;
/* DOCUMENT disk = gyotoDirectionalDisk(...)
            disk, member=value...

    This is a subkind of gyoto_ThinDisk.
    The disk is 2D with emitted intensity depending
    only on the coordinate radius and
    on the direction of emission cos(i), where i
    is the angle between the direction of emission and
    the local normal. Intensity is given
    in an array EMISSION, array(double, NNU, NI, NR),
    with NNU the number of emitted frequencies,
    NI the number of direction cosines, NR the number
    of radii.

    The radii and direction cosine must be supplied
    in RADIUS and COSI.

   KEYWORDS:

    fitsread="filename.fits"  read pattern from FITS file.
    
    fitswrite="filename.fits" write pattern to FITS file.

    copyintensity=EMISSION
                * if EMISSION is nil, retrieve the surface brightness
                  cube;
                * if EMISSION==0, free the cube;
                * if EMISSION is an array of NNU x NPHI x NR doubles,
                  attach (copy) this array into DISK as the surface
                  brightness cube. If this cube doesn't have the same
                  dimensions as the previously set one, the velocity
                  and radius arrays will also be freed (as they have
                  inconsistent dimensions).

    copygridfreq=FREQ
                a NNU element vector containing the frequencies of emission

    copygridcosi=COSI
                a NI element vector containing the direction cosine

    copygridradius=RADIUS
                a NR element vector
                specifying the R coordinate of the grid points

   SEE ALSO:
    gyoto_PatternDisk
*/
gyoto_namespace, DirectionalDisk=gyoto_DirectionalDisk;

//// DISK3D
extern _gyoto_Disk3D_register_as_Astrobj;
_gyoto_Disk3D_register_as_Astrobj;
extern gyoto_Disk3D;
/* DOCUMENT disk = gyoto_Disk3D(...)
            disk, member=value...

    Geometrically thick disk. The grid for the pattern is set by 4
    numbers: NPHI (number of grid points in the azimuthal direction),
    REPEATPHI if the pattern must be repeated several times in the
    azimuthal direction (i.e. if the angular periodicity of the
    pattern is a fraction of 2*pi), NZ (number of grid points in
    the vertical direction), and NR (number of grid points in
    the radial direction). The disk extends from RIN to
    ROUT with a regular spacing along all directions.

    The pattern is specified by the surface brightness which can be
    computed from EMISSQUANT which is typically the temperature
    at NNU frequencies going from NU0 to NU0*DNU*(NNU-1). The cube
    EMISSQUANT is an array(double, NNU, NPHI, NZ, NR).

    Optional opacity field may be provided, with same dimensions as
    emissquant.

    The fluid velocity field must be specified
    with VELOCITY==array(double, 3, NPHI, NZ, NR).
    VELOCITY(1,..)==dphi/dt; VELOCITY(2,..)==dz/dt;
    VELOCITY(3,..)==dr/dt.

   KEYWORDS:

    fitsread="filename.fits"  read pattern from FITS file.
    
    fitswrite="filename.fits" write pattern to FITS file.

    repeatphi=N the pattern angular periodicity is 2*pi/N.

    nu0=        first frequency (Hz)

    dnu=        frequencty spacing (Hz)

    rin=        inner radius

    rout=       outer radius

    zmin=       smallest z value (if >=0 then the disk is assumed to be
                symmetric by z->-z transformation) 

    zmax=       biggest z value

    copyemissquant=EMISSQUANT
                * if EMISSQUANT is nil, retrieve the
                  cube;
                * if EMISSQUANT==0, free the cube;
                * if EMISSQUANT is an array of NNU x NPHI x NZ x NR doubles,
                  attach (copy) this array into DISK.
                  If this cube doesn't have the same
                  dimensions as the previously set one, the velocity
                  array will also be freed (as it has
                  inconsistent dimensions).

    copyopacity=OPACITY
                same as COPYEMISSQUANT for opacity field.
    
    copyvelocity=VELOCITY
                same as COPYEMISSQUANT but to attach the fluid velocity
                field, a 3 x NPHI x NZ x NR array where
                VELOCITY(1,..)==dphi/dt, VELOCITY(2,..)==dz/dt
                and VELOCITY(3,..)==dr/dt

   SEE ALSO:
    gyoto_Astrobj, gyoto_PatternDisk
*/
gyoto_namespace, Disk3D=gyoto_Disk3D;

//// POLISHDOUGHNUT

extern _gyoto_PolishDoughnut_register_as_Astrobj;
_gyoto_PolishDoughnut_register_as_Astrobj;
/* xDOCUMENT _gyoto_<KIND>_register_as_Astrobj
      To be called exactly once in gyoto_<KIND>.i
*/
   
extern gyoto_PolishDoughnut;
/* DOCUMENT pd = gyoto_PolishDoughnut([filename, ][members=values])

         or pd, members=values
         or value = pd(member=)

   PURPOSE:
     Create and manipulate GYOTO PolishDoughnut
     objects. PolishDoughnuts are toroidal accretion structures in
     KerrBL metric. They can be optically thin or optically thick and
     the spectrum can be computed.

     For basic syntax, see "help, gyoto" and "help, gyoto_Astrobj".

   MEMBERS:
     In addition to those listed in "help, gyoto_Astrobj":
     
        lambda= (a double) the lambda parameter
        tempratio, centraldensity, centraltempovervirial, beta
        metric= a gyoto_Metric object (KerrBL only)

   METHODS:
     Some characteristic values can be retrieved with
     value=pd(keyword=):
        l0, Wsurface, Wcentre, rcusp, rcentre

   SEE ALSO: gyoto, gyoto_Astrobj, gyoto_KerrBL
 */
gyoto_namespace, PolishDoughnut=gyoto_PolishDoughnut;

/////// SPECTRUM KIND ///////

/* PowerLaw */

extern gyoto_PowerLawSpectrum;
/* DOCUMENT spectrum = gyoto_PowerLawSpectrum([filename, ][members=values]);
            spectrum, members=values;
            value = spectrum(member=);

            Inu = spectrum(nu)

   SEE ALSO: gyoto, gyoto_Spectrum
*/
gyoto_namespace, PowerLawSpectrum=gyoto_PowerLawSpectrum;

extern _gyoto_PowerLawSpectrum_register_as_Metric;
_gyoto_PowerLawSpectrum_register_as_Metric;

/* BlackBody */

extern gyoto_BlackBodySpectrum;
/* DOCUMENT spectrum = gyoto_BlackBodySpectrum([filename, ][members=values]);
            spectrum, members=values;
            value = spectrum(member=);

            Inu = spectrum(nu)

   SEE ALSO: gyoto, gyoto_Spectrum
*/
gyoto_namespace, BlackBodySpectrum=gyoto_BlackBodySpectrum;

extern _gyoto_BlackBodySpectrum_register_as_Metric;
_gyoto_BlackBodySpectrum_register_as_Metric;

