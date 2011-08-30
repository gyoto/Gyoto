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
     
       coord=prime2tdot=pos, vel: COORD is the 8-vector where
              COORD[1-4]==POS and COORD[5-8] is the 4-velocity
              corresponding to the 3-velocity VEL;

       
       nullifycoord=pos, vel       return nullified (photon) coord tangent
                                    to vel at pos.

   SET KEYWORDS:
     List of set-like keywords ("set" is never specified). Specific
     Metric kinds may recognize more:
       mass=new_mass                gyoto_Metric_setMass, gg, new_mass

   SEE ALSO: gyoto, gyoto_KerrBL, gyoto_KerrKS
 */


extern gyoto_Metric_setMass;
extern gyoto_Metric_getMass;
/* DOCUMENT        gyoto_Metric_setMass, gg, mass
         or        gg, mass=mass
         or mass = gyoto_Metric_getMass(gg)
         or mass = gg(get_mass=1)

     Set or get mass from a from Gyoto Metric object GG
     previously created using e.g. gyoto_Kerr_new.

   SEE ALSO: gyoto_Metric
*/

extern gyoto_Metric_gmunu;
func gyoto_Metric_g(met, x)
/* DOCUMENT coef = gyoto_Metric_gmunu(met, x, mu, nu)
         or coef = met(get_gmunu=x, mu, nu)
         or    g = gyoto_Metric_g(met, x)
   
     Get coefficient COEF of indices MU and NU of Gyoto Metric MET at
     position X.

     In the first form, a single coeeficient is computed and returned
     as a Yorick scalar.
     
     In the third form, all the coefficients are returned as a 4x4
     Yorick array.

     Warning: as usual, the time coordinate is given first. Since
     Yorick indices start with 1, mu and nu are between 1 and 4. For
     instance, the time coefficient of the metric (usually g00) is
     given by mu=nu=1, or g(1, 1).

     INPUTS:
      MET: an Gyoto Metric object, created using e.g. gyoto_Kerr_new();
      X: 4-(or more)-element vector, coordinates of the point where
         the metric is to be evaluated;
      MU, NU: indices of the metric to retrieve, each between 1 and 4;
   
   SEE ALSO: gyoto_Metric
 */
{
  g=array(double, 4, 4);
  for (mu=1;mu<=4;++mu)
    for (nu=1;nu<=4;++nu)
      g(mu,nu)=gyoto_Metric_gmunu(met, x, mu, nu);
  return g;
}

extern gyoto_Metric_SysPrimeToTdot;
/* DOCUMENT tdot=gyoto_Metric_SysPrimeToTdot(gg, coord, v)
         or tdot=gg(get_tdot=coord, vel)

    Compute tdot component of quadri-velocity corresponding to a given
    3-velocity v.

   INPUTS:
    GG: a Gyoto Metric object (e.g. from gyoto_Kerr_new())
    COORD: a quadruplet [x0, x1, x2, x3] giving the position
    V: a triplet [dx1/dt, dx2/dt, dx3/dt]

   OUTPUT:
    TDOT = dt/dtau
   
   SEE ALSO: gyoto_Metric
 */
extern gyoto_Metric_get_refCount;
// a debugging function to read the C++ reference counter
// useless from Yorick
