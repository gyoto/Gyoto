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

// From yutils, for tic() and tac()
#include "util_fr.i"

if (get_env("GYOTO_CHECK_NODISPLAY"))
  window = pause = fma = winkill = pli = plg = noop;

write, format="%s", "New scenery... ";
sc=gyoto_Scenery();
write, format="%s\n", "done.";

write, format="%s\n", "Printing this Scenery:";
sc;

write, format="%s %i\n", "Pointer to this Scenery:", sc(get_pointer=1);

write, format="%s", "New scenery, setting only \"time\"... ";
sc2=gyoto_Scenery(screen=gyoto_Screen(time=1));
write, format="%s\n", "done.";

write, format="%s\n", "Printing this Scenery:";
sc2;
write, format="%s %i\n", "Pointer to this Scenery:", sc2(get_pointer=1);

write, format="%s", "Attaching metric to scenery... ";
sc, metric=gyoto_KerrBL();
write, format="%s\n", "done.";

write, format="%s\n", "Printing this Scenery:";
sc;
write, format="%s", "Pointer to this Scenery: ";
write, format="%i\n", sc(get_pointer=1);

write, format="%s", "Creating star... ";
ao=gyoto_Star(metric=sc(metric=), radius=0.5,
              initcoord=[0,6,pi/2.,0], [0,1e-3,0]);
write, format="%s\n", "done.";

write, format="%s", "Attaching astrobj to scenery... ";
sc, astrobj=ao;
write, format="%s\n", "done.";

write, format="%s", "Retrieving astrobj... ";
ao2=sc(astrobj=);
write, format="%s\n", "done.";

write, format="%s", "Setting time... ";
noop,sc(screen=)(time=10000);
write, format="%s\n", "done.";
write, format="%s %e\n", "Checking time:", (time=sc(screen=)(time=));
if (time!=10000) error, "CHECK FAILED";

// write, format="%s", "Setting tmin... ";
// rien=sc(screen=)(tmin=-10000);
// write, format="%s\n", "done.";
// write, format="%s %e\n", "Checking tmin:", (tmin=sc(screen=)(get_tmin=1));
// if (tmin!=-10000) error, "CHECK FAILED";


// write, format="%s", "Setting dtau... ";
// sc, dtau=0.1;
// write, format="%s\n", "done.";
// write, format="%s %e\n", "Checking dtau:", (dtau=sc(get_dtau=1));
// if (dtau!=0.1) error, "CHECK FAILED";


write, format="%s", "Setting field-of-view... ";
noop,sc(screen=)(fov=pi/4.);
write, format="%s\n", "done.";
write, format="%s %e\n", "Checking field-of-view:", (fov=sc(screen=)(fov=));
if (fov!=pi/4) error, "CHECK FAILED";

write, format="%s", "Setting resolution... ";
noop,sc(screen=)(resolution=16);
write, format="%s\n", "done.";
write, format="%s %i\n", "Checking resolution:", (res=sc(screen=)(resolution=));
if (res!=16) error, "CHECK FAILED";


write, format="%s", "Setting inclination... ";
noop,sc(screen=)(inclination=pi/3.);
write, format="%s\n", "done.";
write, format="%s %e\n", "Checking inclination:",
  (incl=sc(screen=)(inclination=));
if (incl!=pi/3) error, "CHECK FAILED";

if (gyoto_haveXerces()) {
  write, format="%s", "Writing XML description... ";
  sc,xmlwrite="test.xml";
  write, format="%s\n", "done.";
  remove, "test.xml";

  write , format="%s", "Reading Scenery from XML description... ";
  sc3=gyoto_Scenery("../doc/examples/example-moving-star.xml");
 } else {
  write , format="%s", "No Xerces, creating scenery from scratch... ";
  gg  = gyoto_KerrBL();
  spectro = gyoto_Spectrometer(kind="wave", nsamples=1, band=[2e-6, 2.4e-6]);
  scr = gyoto_Screen(metric=gg,
                     observerpos=[1000., 100., 0.78, 0. ],
                     time=1000.,
                     resolution=128,
                     spectro=spectro);
  ao = gyoto_Star(metric=gg, radius=2, rmax=50,
                  initcoord=[600,9,1.5707999999999999741,0],[0., 0., 0.037037])
  sc3=gyoto_Scenery(metric=gyoto_KerrBL(),
                    screen=scr,
                    astrobj=ao,
                    tmin=0.);
 }
noop, sc3(screen=)(resolution=32);
noop, sc3(astrobj=)(radius=2);
write, format="%s\n" , "done.";

write, format="%s", "Ray-tracing on 1 thread (sc())... \n";
sc3, nthreads=1;
tic;
im1 = sc3(,,"Intensity"); // raytrace
tac();
pli, im1;
pause, 1000;
write, format="%s\n" , "done.";

write, format="%s", "Ray-tracing on 2 threads (sc())... \n";
sc3, nthreads=2;
tic;
im1 = sc3(,,"Intensity"); // raytrace
tac();
window,1; pli, im1;
pause, 1000;
write, format="%s\n" , "done.";

write, format="%s", "Ray-tracing on 2 threads (gyoto_Scenery_rayTrace())... \n";
sc3, nthreads=2, quantities="Intensity";
//tic;
im1 = gyoto_Scenery_rayTrace(sc3);
window,2; pli, im1;
pause, 1000;
//tac();
write, format="%s\n" , "done.";

write, format="%s", "Ray-tracing on 1 thread... \n";
sc3, nthreads=1;
tic;
im1 = gyoto_Scenery_rayTrace(sc3);
tac();
write, format="%s\n" , "done.";


/* write, format="%s", "Ray-tracing on adaptive grid... ";
   data = gyoto_Scenery_adaptive_raytrace(sc3, 4);
   fma;
   pli, data(,,3), cmax=100;
   write, format="%s\n" , "done."; */

write, format="%s", "Cloning...";
sc4=sc3(clone=);
write, format="%s\n", "DONE.";

write, format="%s\n", "Printing clone:";
sc4;

ph = gyoto_Photon(initcoord=sc3, 6, 19);
ph(is_hit=);


if (batch()) {

  // Free memory for easier checking with valgrind
  data=[];
  sc4=[];
  sc3=[];
  sc2=[];
  sc=[];
  ao=ao2=[];
  pause, 1000;
  winkill;
 }

