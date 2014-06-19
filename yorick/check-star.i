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

begin_section, "Star Astrobj", "in KerrBL metric";

write, format="%s\n", "Attempting star construction:";
st=gyoto_Star();
write, format="%s", "Printing star:";
st;

write, format="%s", "Cloning...";
st2=st.clone();
write, format="%s\n", "DONE.";

write, format="%s", "Printing clone:";
st2;

write, format="%s", "Gyoto::Star pointer is at address: "; st();

write, format="%s", "Setting radius... ";
st, radius=0.5;
write, format= "%s\n", "done.";

write, format="%s", "Getting radius... ";
if (st.radius()!=0.5) error, "CHECK FAILED";
write, format= "%s\n", "done.";

write, format="%s", "Setting metric... ";
st, metric=gyoto_KerrBL(spin=0.);
write, format= "%s\n", "done.";

write, format="%s", "Changing metric spin... ";
noop, st.metric.spin(0.7);
if (st.metric.spin()!=0.7) error, "CHECK FAILED";
gg=st.metric(); gg, spin=0.995;
if (st.metric.spin()!=0.995) error, "CHECK FAILED";
write, format= "%s\n", "done.";

write, format="%s", "Setting initial condition... ";
st, initcoord=[0, 10.791, 1.5708, 0], [0, 0, 0.0166637];
write, format="%s\n", "done.";

doing, "Setting integrator";
if (gyoto_haveBoost()) {
  st, integrator="runge_kutta_fehlberg78";
  st, deltamaxoverr=0.1;
 } else {
  noop, st(metric=)(deltamaxoverr=0.1);
 }
done;

write, format="%s", "Computing orbit... ";
st, xfill=800;
write, format="%s\n", "done.";

write, format="%s", "Instanciating Screen... ";
screen=gyoto_Screen(metric=st.metric());
write, format="%s\n", "done.";


write, format="%s", "Retrieving projected orbit... ";
data=st(get_skypos=screen);
write, format="%s\n", "done.";

write, format="%s\n", "Printing Star object:";
st;

if (!nodisplay) {
  write, format="%s\n", "Check it out (pausing for 1s)!";
  plg,data(,2), data(,1);
  limits;
  pause, 1000;
  winkill;
 }


write, format="%s", "All in one call... ";
data2=gyoto_Star(radius=0.5,
               metric=gyoto_KerrBL(spin=0.995),
               initcoord=[0, 10.791, 1.5708, 0], [0, 0, 0.0166637],
               xfill=800
               )(
                 get_skypos=(screen2=gyoto_Screen(metric=gyoto_KerrBL(spin=0.995)))
                 );
write, format="%s\n", "done.";

if (!nodisplay) {
  write, format="%s\n", "Check it out (pausing for 1s)!";
  plg,data2(,2), data2(,1);
  pause, 1000;
  winkill;
 }

// Free memroy to check with valgrind;
//data=[];
//st=[];

end_section, "Star Astrobj", "in KerrBL metric";
