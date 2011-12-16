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

#include "gyoto.i"
#include "gyoto_std.i"

aa=0.;

gg=gyoto_KerrBL(spin=aa);

write, format="%s", "Creating Photon: ";
ph = gyoto_Photon();
if (ph())
  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";

write, format="%s", "Attaching metric: ";
ph,metric=gg;
if (ph())
  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";

// initial conditions :
ri=6;
thetai=pi/2.;
phii=0.;

ttravel=r0;
ti=0.;
pri=0.;//canonical momentum
pthetai=0.;
yinit=[ti, ri, thetai, phii, ti, pri, pthetai];
cst=[1,0.921103,2.,0.];//4 Kerr cst of motion a, E, L, Q

write, format="%s", "Checking gyoto_Kerr_MakeCoord: ";
coord=gg(makecoord=yinit, cst);
coord;
//if (abs((coord-[0,10.791,1.5708,0,1.12641,0,0,0.0187701]))(max)<1e-6)
if (abs((coord-[0,6,1.5708,0,1.38165,0,0,0.0555556]))(max)<1e-5)
  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";

pos=coord(:4);
v= coord(6:)/coord(5);
write, format="%s", "Checking gyoto_Star: ";
st=gyoto_Star(metric=gg, radius=1., initcoord=pos, v);
write, format="%s\n","done.\n";

write, format="%s\n", "Trying gyoto_Star_xFill";
st(xfill=770);
"done";
gyoto_Star_xFill(st,770.);

//Computing position of star at a given proper time :
//time=212.4034;//proper time
//write, format="%s\n", "Checking gyoto_Star_position: ";
//pos=gyoto_Part_position(st,time);
//pos;
//if (abs(pos-[10.5718661339679, 1.57079398752261, 59.5795847453848])(max)<1e-5)
//  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";

gyoto_Star_get_xyz,st,x,y;


gyoto_Star_get_coord, st, t, r, theta, phi;
gyoto_Star_get_dot, st, tdot, rdot, thetadot, phidot;
gyoto_Star_get_prime, st, rp, thetap, phip;

write, format="%s", "Checking gyoto_Metric_SysPrimeToTdot: ";
tdotbis=array(double,numberof(t));
for (n=1; n<= numberof(t); ++n)
  tdotbis(n)=gg(prime2tdot= [t(n), r(n), theta(n), phi(n)],
                [rp(n), thetap(n), phip(n)]
                );
if (max (abs( (tdot-tdotbis)/tdot ) ) < 2e-3)
  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";


write, format="%s", "Checking gyoto_Metric_g: ";
norm=array(double, numberof(t));
for (n=1; n<= numberof(t); ++n) {
  g=gg([t(n), r(n), theta(n), phi(n)]);
  qvel=[tdot(n), rdot(n), thetadot(n), phidot(n)];
  norm(n)=sum(g*qvel(,-)*qvel(-,));
 }
if (max(abs(norm+1)) < 3e-3)
  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";

// Ray tracing

write, format="%s", "Checking gyoto_Metric_setObserverPos: ";
screen=gyoto_Screen(metric=gg, observerpos=[1000., 100., 0.05, 0.]);
write, format="%s\n","done.\n";
write, format="%s", "Checking gyoto_Metric_setSpin: ";
gg,spin=0.;
write, format="%s\n","done.\n";
write, format="%s", "Checking gyoto_Star(): ";
orbit=gyoto_Star(metric=gg, radius=2,
                 initcoord=[600, 6, 1.57, 0], [0, 0, 0.068041381745])
write, format="%s\n","done.\n";

N=21;
delta=pi/(10.*N);
screen, fov=pi/10., resolution=21;
write, format="%s", "Checking gyoto_Photon_new: ";
ph=gyoto_Photon_new();
write, format="%s\n","done.\n";

i=15; j=9;
xscr=delta*(i-(N+1)/2.);
yscr=delta*(j-(N+1)/2.);
write, format="%s", "Checking gyoto_Photon_setInitialCondition: ";
gyoto_Photon_setInitialCondition, ph, gg, orbit, screen, -xscr, yscr;
ph2=gyoto_Photon(metric=gg, astrobj=orbit);
ph1=gyoto_Photon(metric=gg, astrobj=orbit);
ph1, initcoord=screen, -xscr, yscr;
ph2, initcoord=screen, i, j;
write, format="%s\n","done.\n";
write, format="%s", "Checking gyoto_Photon_setDelta: ";
gyoto_Photon_setDelta, ph, 1.;
write, format="%s\n","done.\n";
write, format="%s", "Checking gyoto_Photon_hit: ";
//if(gyoto_Photon_hit(ph, 0.))
if( ph(is_hit=1) && ph1(is_hit=1) && ph2(is_hit=1))
  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";

"_________________________";

hitmap2=hitmap1=hitmap=array(0, N, N);
for (i=1; i<=N; i++) {
  write , format="*** Column %i ***\n", i; 
  xscr=delta*(i-(N+1)/2.);
  for (j=1; j<=N; j++) {
    yscr=delta*(j-(N+1)/2.);
    gyoto_Photon_setInitialCondition, ph, gg, orbit, screen, -xscr, yscr;
    ph1, initcoord=screen, -xscr, yscr;
    ph2, initcoord=screen, i, j;
    gyoto_Photon_setDelta, ph, 1.;
    hitmap(i,j)=gyoto_Photon_hit(ph, 0.);
  }
 }
ph1=ph2=[];
"_________________________";

ph=[];
st=[];
orbit=[];
screen=[];

gg=[];

write, format= "%s\n"," ALL TESTS PASSED";
