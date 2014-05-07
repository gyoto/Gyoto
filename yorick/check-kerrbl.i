/*
    Copyright 2011, 2014 Thibaut Paumard

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

begin_section, "KerrBL metric", "Kerr in Boyer-Lindquist coordinates";

aa=0.995;
write, format="%s", "Checking gyoto_KerrBL: ";
gg=gyoto_KerrBL(spin=aa);
write, format="%s\n","done.\n";

check_gmunu, gg;
check_gmunu_up, gg;
check_christoffels, gg;

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
if (gg.spin!=0.9) error, "CHECK FAILED";
write, format="%s\n", "done.";

//write, format="%s", "Pointer to this Metric object: gg()==";
gg();

write, format="%s\n", "Printing object. \"gg\" yields: ";
gg;

write, format="%s", "Setting deltamin... ";
gg, deltamin=40;
if (gg.deltamin!=40) error, "CHECK FAILED";
write, format="%s\n", "done.";

write, format="%s", "Setting deltamax... ";
gg, deltamax=400;
if (gg.deltamax!=400) error, "CHECK FAILED";
write, format="%s\n", "done.";

write, format="%s", "Setting difftol... ";
gg, difftol=1e-3;
if (gg.difftol!=1e-3) error, "CHECK FAILED";
write, format="%s\n", "done.";

write, format="%s", "Setting deltamaxoverr... ";
gg, deltamaxoverr=2e-3;
if (gg.deltamaxoverr!=2e-3) error, "CHECK FAILED";
write, format="%s\n", "done.";

write, format="%s", "Cloning... ";
gg2=gg(clone=);
if (gg2.deltamax!=400 || gg2.deltamin!=40 || gg2.difftol!=1e-3)
  error, "CHECK FAILED";
write, format="%s\n", "done.";

// Free memory for testing with valgrind
gg2=[];
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

end_section, "KerrBL metric", "Kerr in Boyer-Lindquist coordinates";
