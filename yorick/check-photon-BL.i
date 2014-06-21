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

#include "check-helpers.i"

begin_section, "Photon in KerrBL metric";

aa=0.;

gg=gyoto_KerrBL(spin=aa);
gg2=gg.clone;

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

doing, "Setting metric spin";
gg, spin=0.95;
done;
doing, "Computing initial coordinate";
coord=gg.makecoord(yinit, cst);
done;

pos=coord(:4);
v= coord(6:)/coord(5);
doing, "Checking gyoto_Star";
st=gyoto_Star(metric=gg, radius=1., initcoord=pos, v);
done;

doing, "Trying gyoto_Star_xFill";
st,xfill=770.;
done;

//Computing position of star at a given proper time :
//time=212.4034;//proper time
//write, format="%s\n", "Checking gyoto_Star_position: ";
//pos=gyoto_Part_position(st,time);
//pos;
//if (abs(pos-[10.5718661339679, 1.57079398752261, 59.5795847453848])(max)<1e-5)
//  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";


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
gg, setparameter="GenericIntegrator";
i=15; j=9;
xscr=delta*(i-(N+1)/2.);
yscr=delta*(j-(N+1)/2.);
write, format="%s", "Checking gyoto_Photon_setInitialCondition: ";
ph2=gyoto_Photon(metric=gg, astrobj=orbit);
ph1=gyoto_Photon(metric=gg, astrobj=orbit);
ph1, initcoord=screen, -xscr, yscr;
ph2, initcoord=screen, i, j;
write, format="%s\n","done.\n";
write, format="%s", "Checking gyoto_Photon(delta=1): ";
write, format="%s\n","done.\n";
write, format="%s", "Checking gyoto_Photon(is_hit=1): ";
if( ph1(is_hit=1) && ph2(is_hit=1))
  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";

"_________________________";

hitmap2=hitmap1=hitmap=array(0, N, N);
ph=gyoto_Photon(metric=gg, astrobj=orbit);
for (i=1; i<=N; i++) {
  write , format="*** Column %i ***\n", i; 
  xscr=delta*(i-(N+1)/2.);
  for (j=1; j<=N; j++) {
    yscr=delta*(j-(N+1)/2.);
    ph, initcoord=screen, -xscr, yscr;
    ph, tmin=0.;
    ph1, initcoord=screen, -xscr, yscr;
    ph2, initcoord=screen, i, j;
    ph, delta=1.;
    hitmap(i,j)=ph(is_hit=1);
  }
 }
ph2=[];


// Check that changing spin can be done on attached metric
ph1, metric=ph1.metric();
ph1, xfill=0.;
txyz=ph1.get_txyz();
  plg, txyz(2,), txyz(1,);
ph2=ph1.clone();
gg2=gg.clone();
write, format="%s", "Mutating metric spin... ";
gg, spin=0.5;
write, format="%s\n", "done.";
gg2, spin=0.5;
ph2, metric=gg2;


ph1=ph2=[];

// CLONES AND HOOKS

"_________________________";

ph=[];
st=[];
orbit=[];
screen=[];

gg=[];

end_section, "Photon in KerrBL metric";
