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

aa=0.995;
write, format="%s", "Checking gyoto_KerrBL: ";
gg=gyoto_KerrBL(spin=aa);
write, format="%s\n","done.\n";

write, format="%s", "Checking gyoto_Metric_get_RefCount: ";
if (gyoto_Metric_get_refCount(gg)==1)
  write, format="%s\n","done.\n"; else error, "PREVIOUS CHECK FAILED";

// initial conditions :
ri=10.791;
thetai=1.5708;
phii=0.;
ti=0.;
pri=0.;//canonical momentum
pthetai=0.;
cst=[1,0.921103,2.,0.];//4 Kerr cst of motion mu, E, L, Q
yinit=[ti,ri,thetai,phii,-cst(2),pri,pthetai,cst(3)];

write, format="%s", "Checking makecoord: ";
coord=gg(makecoord=yinit, cst);
if (abs((coord-[0,10.791,1.5708,0,1.12641,0,0,0.0187701]))(max)<1e-6)
  write, format="%s\n","done."; else error, "PREVIOUS CHECK FAILED";

write, format="%s\n", "Creating metric using gyoto_Kerr(spin=0.7)";
gg=gyoto_KerrBL(spin=0.7);
if (gyoto_KerrBL(gg, spin=)!=0.7) error, "CHECK FAILED";
write, format="%s\n","done.";

write, format="%s", "Setting spin... ";
gg, spin=0.9;
write, format="%s\n", "done.";
  
write, format="%s", "Getting spin... ";
if (gg(spin=)!=0.9) error, "CHECK FAILED";
write, format="%s\n", "done.";

//write, format="%s", "Pointer to this Metric object: gg()==";
gg();

write, format="%s\n", "Printing object. \"gg\" yields: ";
gg;


write, format="%s", "Checking gmunu... ";
for (i=1; i<=4; ++i)
  for (j=1; j<=4; ++j)
    if (gg([0, 5, 0, pi/2], i, j) !=
        gyoto_Metric_gmunu( gg, [0, 5, 0, pi/2], i, j))
      error, "CHECK FAILED";
write, format="%s\n", "done.";

write, format="%s", "Checking prime2tdot... ";
if (gg(prime2tdot=[0, 6, pi/2, pi/2], [0., 0.01, 0.01]) !=
    gyoto_Metric_SysPrimeToTdot(gg, [0, 6, pi/2, pi/2], [0., 0.01, 0.01]))
  error, "CHECK FAILED";
write, format="%s\n", "done.";

// Free memory for testing with valgrind
gg=[];
aa=[];
ri=[];
thetai=[];
phii=[];
ti=[];
pri=[];
pthetai=[];
cst=[];
yinit=[];


write, format="\n%s\n", "ALL TESTS PASSED";
