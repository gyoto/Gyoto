/*
    Copyright 2011, 2013-2015 Thibaut Paumard

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

// NODISPLAY implies batch mode
if (get_env("GYOTO_CHECK_NODISPLAY")) {
  batch, 1;
  __xytitles=xytitles; __fma=fma; __winkill=winkill; __pli=pli; __plg=plg;
  __pause=pause; __window=window;
  xytitles = fma = winkill = pli = plg = pause = window = noop;
 }

#include "check-helpers.i"

begin_section, "basic functionality";

aa=0.995;
write, format="%s", "Checking gyoto_Kerr_new: ";
gg=gyoto_KerrBL(spin=aa);
write, format="%s\n","done.\n";

// initial conditions :
ri=10.791;
thetai=1.5708;
phii=0.;
ti=0.;
pri=0.;//canonical momentum
pthetai=0.;
cst=[1,0.921103,2.,0.];//4 Kerr cst of motion mu, E, L, Q
yinit=[ti,ri,thetai,phii,-cst(2),pri,pthetai,cst(3)];

write, format="%s", "Checking gyoto_Kerr_MakeCoord: ";
coord=gg(makecoord= yinit, cst);
if (abs((coord-[0,10.791,1.5708,0,1.12641,0,0,0.0187701]))(max)<1e-6)
  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";

pos=coord(:4);
v= coord(6:)/coord(5);
write, format="%s", "Checking gyoto_Star(): ";
st=gyoto_Star(metric=gg, radius=1., initcoord=pos, v);
write, format="%s\n","done.\n";

write, format="%s\n", "Trying gyoto_Star_xFill";
st,xfill=770.;

//Computing position of star at a given proper time :
//time=212.4034;//proper time
//write, format="%s\n", "Checking gyoto_Star_position: ";
//pos=gyoto_Part_position(st,time);
//pos;
//if (abs(pos-[10.5718661339679, 1.57079398752261, 59.5795847453848])(max)<1e-5)
//  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";

txyz=st.get_txyz; dates=txyz(,1); x=txyz(,2); y=txyz(,3);
coords=st.get_coord();
primes=st.get_prime();

write, format="%s", "Checking gg(prime2tdot= pos, vel): ";
N=dimsof(coords)(2);
tdotbis=array(double,N);
for (n=1; n<= N; ++n)
  tdotbis(n)=gg (prime2tdot=coords(n,1:4), primes(n,));
tdot=coords(,5);
if (max (abs( (tdot-tdotbis)/tdot ) ) < 1e-4)
  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";


write, format="%s", "Checking gg(position): ";
norm=array(double, N);
for (n=1; n<= N; ++n) {
  g = gg (coords(n,));
  qvel=coords(n,5:);
  norm(n)=sum(g*qvel(,-)*qvel(-,));
 }
if (max(abs(norm+1)) < 1e-4)
  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";

// Ray tracing

write, format="%s", "Checking gyoto_Screen(observerpos=...): ";
screen=gyoto_Screen(metric=gg, observerpos=[1000., 100., 0.05, 0.]);
write, format="%s\n","done.\n";
write, format="%s", "Checking gyoto_Metric(spin=): ";
gg, spin=0.;
write, format="%s\n","done.\n";
write, format="%s", "Checking gyoto_Star(): ";
orbit=gyoto_Star(metric=gg, radius=2,
                 initcoord=[600, 6, 1.57, 0], [0, 0, 0.068041381745])
write, format="%s\n","done.\n";

N=51;
delta=pi/(10.*N);
write, format="%s", "Checking gyoto_Photon(): ";
ph=gyoto_Photon();
write, format="%s\n","done.\n";

i=35; j=19;
xscr=delta*(i-(N+1)/2.);
yscr=delta*(j-(N+1)/2.);

write, format="%s", "Checking gyoto_Photon_setInitialCondition: ";
ph, metric=gg, astrobj=orbit, initcoord=screen, -xscr, yscr;
write, format="%s\n","done.\n";

write, format="%s", "Checking gyoto_Photon_setDelta: ";
ph, delta=1.;
write, format="%s\n","done.\n";
write, format="%s", "Checking gyoto_Photon_hit: ";
if(ph(is_hit=1))
  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";

"_________________________";

hitmap=array(0, N, N);
ph, tmin=0.;
for (i=1; i<=N; i++) {
  write , format="*** Column %i ***\n", i; 
  xscr=delta*(i-(N+1)/2.);
  for (j=1; j<=N; j++) {
    yscr=delta*(j-(N+1)/2.);
    ph, metric=gg, astrobj=orbit, initcoord=screen, -xscr, yscr;
    // gyoto_Photon_setDelta, ph, 1.;
    hitmap(i,j)=ph(is_hit=1);
  }
 }
"_________________________";

screen=[];
ph=[];
st=[];
orbit=[];

gg=[];

end_section, "basic functionality";

#include "check-photon-BL.i"
#include "check-scenery.i"
#include "check-kerrbl.i"
#include "check-kerrks.i"
#include "check-minkowski.i"
#include "check-star.i"
#include "check-startrace.i"
#include "check-patterndisk.i"
#include "check-directionaldisk.i"
#include "check-disk3d.i"
#include "check-polish-doughnut.i"
// Don't run check-mpi.i automatically,
// It may hang the computer if not plugged to the network
//#include "check-mpi.i"

write, format="\n\n%s\n%s\n%s\n%s\n\n",
  "  ********************************************",
  "  *             ALL TESTS PASSED             *",
  "  ********************************************",
  "  (You may want to run 'make check-mpi' still)";

if (anyof(get_argv() == "check.i")) quit;
//if (batch()) quit;
