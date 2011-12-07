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
