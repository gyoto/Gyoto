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
kerrbl=KerrBL(spin=0., mass=1*earth_mass, unit="kg");
kerrks=KerrKS(mass=kerrbl.mass(), spin=kerrbl.spin());
minkowski=Metric("Minkowski");
minkowski, mass=kerrbl.mass(), setparameter="Spherical";
metric=kerrbl;

// Define the screen (=camera):
res=128; // resolution of the output image
scr=Screen(inclination=pi/2, paln=pi, argument=pi/2,
           metric=metric, resolution=res, anglekind="Rectilinear");
scr,distance=2, unit="m"; // camera is at 2m from the BH
scr, fov=pi/2.;
// scr, fov=2.*atan(36./100.); // ~40°, 36mm film behind a 50mm lense

// Use an empty Astrobj, just set rmax to 10 meters:
ao=Astrobj("Complex", rmax=10./metric.unitlength());

// Build a Scenery with this Metric, Screen and Astrobj:
sc=Scenery(metric=metric, astrobj=ao, screen=scr, nthreads=8);

///// Second, choose a painter to paint the sky:
// To use a jpeg file containing a full-sky, 360°x180° panorama:
// download file:
//   http://farm1.staticflickr.com/192/456185667_adde9d2f8e_o_d.jpg
//painter=gyoto.painters.mk_panorama(img=jpeg_read("456185645_e56abcc2cd_o.jpg")(,,::-1));
// Or a more conventional picture:
//   https://i.pinimg.com/236x/19/53/11/195311bdb5b029cb72214b95833e6dcf.jpg
//painter=gyoto.painters.mk_picture(img=jpeg_read("195311bdb5b029cb72214b95833e6dcf.jpg")(,,::-1), fov=2.*atan(36./100.))
// To use a p-mode-like pattern:
painter=gyoto.painters.mk_p_mode(ntheta=80, nphi=80);


///// Third, simply call the appropriate function:
// matte_paint() can be called directly on sc, but to specify only a
// region or to call matte_paint repeatedly (e.g. on distinct
// painters, of to fiddle with yaw, pitch & roll), we can precompute
// ipct:
nx=res
ny=res*10/16; // use a 16/10 aspect ratio
x1=(res-nx)/2+1;
x2=res+1-x1;
y1=(res-ny)/2+1;
y2=res+1-y1;
ipct=sc(x1:x2,y1:y2,impactcoords=);
//ipct=fits_read("ipct.fits");
//bg=matte_paint(ipct, painter, kind=metric.coordkind());
bg=matte_paint(ipct, painter, kind=metric.coordkind(), yaw=-1.62, roll=-0.32);

///// Fourth, display the image:
fma;
pli, bg;
limits, square=1;

///// Fifth, possibly, save it as a JPEG:
//jpeg_write, "file_out.jpg", bg(,,::-1);
