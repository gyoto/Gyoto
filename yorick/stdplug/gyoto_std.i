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

#include "gyoto_Kerr.i"
#include "gyoto_Star.i"

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
     
    In addition, the standard gyoto_Astrobj members and methods are
    supported.
   
   SEE ALSO: gyoto, gyoto_Astrobj, gyoto_Metric
   
 */

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

//// PAGETHORNEDISK
extern _gyoto_PageThorneDisk_register_as_Astrobj;
_gyoto_PageThorneDisk_register_as_Astrobj;
extern gyoto_PageThorneDisk;
/* DOCUMENT disk = gyotoPageThorneDisk(...)
            disk, member=value...

    This is a subkind of gyoto_ThinDisk with an emission law based on
    Page & Thorne 1974. Works only in Kerr metric (KerrBL or
    KerrKS). The spin is cached when the metric is set: if you later
    change the spin in the metric, use the UPDATESPIN keyword to
    update it, i.e.:

      disk, metric=gyoto_KerrBL();
      noop, disk(metric=)(spin=0.5);
      disk, updatespin=

   KEYWORDS:

    updatespin= [] updated cached value of spin parameter.

   SEE ALSO:
    gyoto_Astrobj, gyoto_ThinDisk, gyoto_KerrBL, gyoto_KerrKS
*/

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
    gyoto_Astrobj, gyoto_ThinDisk, gyotoPageThorneDisk
*/

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

     copyvelocity=VELOCITY
                same as COPYEMISSQUANT but to attach the fluid velocity
                field, a 3 x NPHI x NZ x NR array where
                VELOCITY(1,..)==dphi/dt, VELOCITY(2,..)==dz/dt
                and VELOCITY(3,..)==dr/dt

   SEE ALSO:
    gyoto_Astrobj, gyoto_PatternDisk
*/

/////// SPECTRUM KIND ///////

/* PowerLaw */

extern gyoto_PowerLawSpectrum;
/* DOCUMENT spectrum = gyoto_PowerLawSpectrum([filename, ][members=values]);
            spectrum, members=values;
            value = spectrum(member=);

            Inu = spectrum(nu)

   SEE ALSO: gyoto, gyoto_Spectrum
*/

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

extern _gyoto_BlackBodySpectrum_register_as_Metric;
_gyoto_BlackBodySpectrum_register_as_Metric;
