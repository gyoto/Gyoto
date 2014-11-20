/*
    Copyright 2013-2014 Thibaut Paumard

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

begin_section, "StarTrace Astrobj";

if (gyoto_haveXerces()) {
  sc = gyoto.Scenery(GYOTO_EXAMPLES_DIR+"example-moving-star.xml");
  noop, sc.screen.mask(0); // make sure no mask is set yet
  st = sc.astrobj();
 } else {
  // No XML, build from scratch
  met=gyoto.KerrBL();
  st=gyoto.Star(metric=met, radius=2.,
                initcoord=[600., 9., 1.5707999999999999741, 0],  [0., 0., 0.037037]);
  screen=gyoto.Screen(metric=met,
                      observerpos=[1000., 100., 0.78, 0.],
                      time=1000.,
                      resolution=128,
                      fov=0.1*pi);
  sc=gyoto.Scenery(metric=met, screen=screen, astrobj=st, tmin=0.);
 }

write, format="%s", "Instanciating StarTrace from Star... ";
stt = st(startrace=600, 800);
write, format="%s\n", "done.";

write, format="%s", "Mutating StarTrace... ";
stt, adaptive=0, delta=1, opticallythin=0;
write, format="%s\n", "done.";

sc, astrobj=stt, nthreads=8;

write, format="%s\n", "Ray-tracing StarTrace... ";
tic;
mask = sc(,,"Intensity");
tac();
write, format="%s\n", "done.";
sc, astrobj=st;

write, format="%s\n", "Ray-tracing Star without mask... ";
tic;
im1=sc(,,);
tac();
write, format="%s\n", "done.";

noop, sc.screen.mask(mask);
write, format="%s\n", "Ray-tracing Star with mask... ";
tic;
im2=sc(,,);
tac();
write, format="%s\n", "done.";

end_section, "StarTrace Astrobj";
