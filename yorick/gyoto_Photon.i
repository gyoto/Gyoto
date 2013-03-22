plug_in, "gyoto";
/*
    Copyright 2011, 2013 Thibaut Paumard

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


// PHOTON CLASS
extern gyoto_Photon;
/* DOCUMENT photon = gyoto_Photon([filename], [members=values])
            photon, member=values
            value = photon(member=)
            value = photon(function_method=params)
            photon, subroutine_method=params

   PURPOSE:

     Instanciate and use a single GYOTO photon.

     Photons are mass-less particles following light-like geodesics of
     a metric. For basic concepts, see GYOTO. For ray-tracing,
     gyoto_Scenery() is more appropriate.

   MEMBERS:

     Members can be set with the syntax "photon, member=value" and
     retrieved with the syntax "value=photon(member=)":
     
        metric= a GYOTO Metric (see gyoto_Metric),
            initcoord=scenery,x,y also sets the metric.
        astrobj=a GYOTO Astrobj (see gyoto_Astroj), the target of
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
                SCREEN is a gyoto_Screen, DALPHA and DDELTA specify
                the direction this photon comes from when it reaches
                the screen. DALPHA and DDELTA are in radians and must
                be floating-point values.
            initcoord=SCREEN,I,J
                As above, but I and J are integers specifying the
                pixel of the arrival SCREEN which the photon hits.
            initcoord=SCENERY,DALPHA,DDELTA
            initcoord=SCENERY,I,J
                As above, but specify a gyoto_Scenery instead of a
                gyoto_Screen. The Metric and Astrobj of the Senery
                will also be attached to the Photon.

            Those last ways of specifying the initial conditions are
            very useful to get the trajectory of a specific photon in
            a ray-traced scenery.

        spectro= a gyoto_Spectrometer

     
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
     
   SEE ALSO: gyoto, gyoto_Metric, gyoto_Screen, gyoto_Scenery,
            gyoto_Astrobj
 */

