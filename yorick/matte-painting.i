#!/usr/bin/env yorick
/*
    This script examplifies how to perform matte painting with Gyoto,
    i.e. vizualizing lensing effects on a painted background.

    Copyright 2014 Thibaut Paumard

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

#include "gyoto.i"
#include "jpeg.i"
restore, gyoto;

///// First, define an appropriate Scenery:

// Choose a metric, set its mass:
earth_mass=5.9722e24; // mass of the Earth, in kg
kerrbl=KerrBL(spin=0.999999, mass=1*earth_mass, unit="kg");
kerrks=KerrKS(mass=kerrbl.mass(), spin=kerrbl.spin());
minkowski=Metric("Minkowski");
minkowski, mass=kerrbl.mass(), setparameter="Cartesian";
metric=kerrks;

// Define the screen (=camera):
res=512; // resolution of the output image
scr=Screen(inclination=pi/2, paln=pi, argument=pi/2,
           metric=metric, resolution=res);
scr,distance=2, unit="m"; // camera is at 1m from the BH
scr, fov=1.0;

// Use an empty Astrobj, just set rmax to 10 meters:
ao=Astrobj("Complex", rmax=10./metric.unitlength());

// Build a Scenery with this Metric, Screen and Astrobj:
sc=Scenery(metric=metric, astrobj=ao, screen=scr, nthreads=8);

///// Second, choose a painter to paint the sky:
// To use a jpeg file containing a full-sky, 360°x180° panorama:
painter=gyoto.painters.mk_panorama(img=jpeg_read("456185645_e56abcc2cd_o.jpg")(,,::-1));
// To use a p-mode-like pattern:
//painter=gyoto.painters.mk_p_mode(ntheta=80, nphi=80);

///// Third, simply call the appropriate function:
ipct=sc(,,impactcoords=);
bg=matte_paint(ipct, painter, kind=metric.coordkind(), yaw=1);

///// Fourth, display the image:
fma;
pli, bg;

///// Fifth, possibly, save it as a JPEG:
//jpeg_write, "file_out.jpg", bg(,,::-1);
