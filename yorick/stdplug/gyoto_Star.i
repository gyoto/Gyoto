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
               

     Function-like (retval=st(function_method=parameters)):

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

        get_coord=1
              Retrieve 8-coordinates (in Metric-prefered coord system)
                 data = st(get_coord=1);
                 data(i, )=[t_i, x1_i, x2_i, x3_i,
                            x0dot_i, x1dot_i, x2dot_i, x3dot_i] 
              
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
                
   SEE ALSO: gyoto_Metric, gyoto_Kerr
 */


extern gyoto_Star_xFill;
/* DOCUMENT gyoto_Star_xFill, st, tlim
         or st, xfill=tlim

    Integrate Star orbit from previously set initial condition to tlim.
    
   SEE ALSO: gyoto_Star
 */

//extern gyoto_Star_position;
/* xDOCUMENT pos=gyoto_Star_position(star, t)

     Return position of GYOTO Particule object STAR at time T.
     
   SEE ALSO: gyoto_Star_get_xyz, gyoto_Star_get_t
 */

extern gyoto_Star_getSkyPos;
/* DOCUMENT gyoto_Star_getSkyPos, star, dalpha, ddelta, dDistance
         or da_dd_dD = star(get_skypos=1)
         
     Return orbit projected on sky of GYOTO Star object STAR.
     
   SEE ALSO: gyoto_Star
 */


extern gyoto_Star_get_t;
extern gyoto_Star_get_xyz;
/* DOCUMENT gyoto_Star_get_xyz, star, x, y, z
         or t    = gyoto_Star_get_t(star)
         or txyz = star(get_txyz=1)
         
     Return orbit (cartesian coordinates x, y, z and/or time t) of
     GYOTO Star object STAR. Warning: these coordinates are a
     Cartesian expression of the Boyer-Lindquist coordinates, they do
     not match the Kerr-Schild coordinates.
     
   SEE ALSO: gyoto_Star
 */

extern gyoto_Star_get_coord;
extern gyoto_Star_get_dot;
extern gyoto_Star_get_prime;
/* DOCUMENT gyoto_Star_get_coord, star, t, r, theta, phi
         or gyoto_Star_get_dot, star, tdot, rdot, thetadot, phidot
         or gyoto_Star_get_prime, star, dr_dt, dtheta_dt, dphi_dt
         or coord_dot = star(get_coord=1)

     Get STAR's trajectory, quadri-velocity or classical velocity in
     Boyer-Lindquist coordinates.

     INPUTS:
      STAR: a Gyoto Star object as returned by gyoto_Star();

     OUTPUTS:
      t, r, theta, phi: four arrays of doubles giving the
         Boyer-Lindquist coordinates of the Star along its trajectory,
         as computed when initiated by gyoto_Star();
      tdot, rdot, thetadot, phidot: four arrays of doubles giving the
         Boyer-Lindquist coordinates of the Star's quadri-velocity
         (dx^alpha/dtau, where tau is the proper time);
      dr_dt, dtheta_dt, dphi_dt: rdot/tdot, thetadot/tdot,
         phidot/tdot, the three coordinates of the velocity in a
         coordinate system attached to the observer.
      coord_dot = [t, x0, x1, x2, tdot, x0dot, x1dot, x2dot]
      
   SEE ALSO: gyoto_Star
 */

extern is_gyoto_Star;
/* DOCUMENT ret = is_gyoto_Star( obj )

   Return 1 if OBJ is a GYOTO Star

   SEE ALSO: gyoto gyoto_Astrobj gyoto_Star is_gyoto_Astrobj
*/
