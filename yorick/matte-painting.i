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

matte=[];
// To use a jpeg file as matte:
//   - uncomment the following line;
//   - set the path to the file you want;
//   - set x0, y0 and fov0 below.
//matte=jpeg_read("file_in.jpg")(,,::-1);

// Pixel coordinates of the center of the zone of interest in the
// JPEG image:
x0=2000;
y0=1500;
// Size in pixel of the zone of the JPEG file that should map the
// screen field-of-view:
fov0=3000;

// Choose a metric, set its mass:
earth_mass=5.9722e24; // mass of the Earth, in kg
kerrbl=KerrBL(spin=0., mass=0.1*earth_mass, unit="kg");
kerrks=KerrKS(mass=kerrbl.mass(), spin=kerrbl.spin());
minkowski=Metric("Minkowski");
minkowski, mass=kerrbl.mass(), setparameter="Cartesian";
metric=kerrbl;

// Define the screen (=camera):
res=128; // resolution of the output image
scr=Screen(inclination=pi/2, paln=pi, metric=metric, resolution=res);
scr,distance=1, unit="m"; // camera is at 1m from the BH

// Use an empty Astrobj, just set rmax to 10 meters:
ao=Astrobj("Complex", rmax=10./metric.unitlength());

// Build a Scenery with this Metric, Screen and Astrobj:
sc=Scenery(metric=metric, astrobj=ao, screen=scr, nthreads=8);

// Compute the last coordinates of the "escaping" photons
// (Actually coming from infinity)
data=sc(,,impactcoords=)(9:,,);

// Transform the photon's 4-velocity into equatorial coordinates
if (sc.metric().coordkind()==gyoto.coordkind.spherical) {

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

  mask0=r<1e200;
  
 } else {
  tdot=data(5,,);
  xdot=data(6,,);
  ydot=data(7,,);
  zdot=data(8,,);
  v=[xdot, ydot, zdot]/tdot(,,-);
  mask0=tdot<1e200;
 }

vproj2=v(,,1)^2+v(,,2)^2

vr=sqrt(vproj2+v(,,3)^2);
vph=atan(v(,,2), v(,,1));
vth=atan(sqrt(vproj2), v(,,3));

ind=where(r>1e308);

// A simple analytical matte:
bg=sin(80*vph)*sin(80*vth);
if (numberof(ind)) bg(ind)=bg(min, min);

// If a JPEG was provided:
if (!is_void(matte)) {

  dmatte=dimsof(matte);

  phi0=(vph(0,0)+vph(1,1))*0.5;
  theta0=(vth(0,0)+vth(1,1))*0.5;
  Dphi=(vph(0,0)-vph(1,1));
  Dtheta=(vth(0,0)-vth(1,1));

  x=long((vph-phi0)/Dphi*fov0)+x0;
  y=long((vth-theta0)/Dtheta*fov0)+y0;

  if (anyof((mask=(x<1)))) x(where(mask))=1;
  if (anyof((mask=(y<1)))) y(where(mask))=1;
  if (anyof((mask=(x>dmatte(3))))) x(where(mask))=dmatte(3);
  if (anyof((mask=(y>dmatte(4))))) y(where(mask))=dmatte(4);

  bg=array(char, 3, res, res);
  for(k=1; k<=3; ++k) {
    for (i=1; i<=res; ++i) {
      for (j=1; j<=res; ++j) {
        bg(k,i,j)=matte(k,x(i, j),y(i, j));
      }
    }
    bg(k,,) *= mask0;
  }

 }

// Display the image:
fma;
pli, bg;


// Possibly, save it as a JPEG:
//jpeg_write, "file_out.jpg", bg(,,::-1);
