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
restore, gyoto;

begin_section, "Scenery";

// From yutils, for tic() and tac()
#include "util_fr.i"

write, format="%s", "New scenery... ";
sc=gyoto_Scenery();
write, format="%s\n", "done.";

write, format="%s\n", "Printing this Scenery:";
sc;

write, format="%s %i\n", "Pointer to this Scenery:", sc(get_pointer=1);

// MaxIter is a property of type size_t. The implementation of such
// properties rely on undefined behaviour. Check here that it works in
// practice. If one of these test fail, we have to find a new
// implementation of GYOTO_PROPERTY_SIZE_T.

doing, "Using maxiter()";
max1=sc.maxiter();
done;

doing, "Using MaxIter()";
max2=sc.MaxIter();
done;

doing, "Comparing";
if (max1!=max2) error, "maxiter() and MaxIter() do not yield the same value";
done;

doing, "Checking MaxIter(val)";
if ((sc.MaxIter(25)).maxiter() != 25) error, "MaxIter and maxiter do not agree";
done;

doing, "Checking maxiter(val)";
if ((sc.maxiter(40)).MaxIter() != 40) error, "MaxIter and maxiter do not agree";
done;

// Set back the original value. Phew, GYOTO_PROPERTY_SIZE_T works as
// expected.
sc, MaxIter=max1;

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
ao=gyoto_Star(metric=sc.metric(), radius=0.5,
              initcoord=[0,6,pi/2.,0], [0,1e-3,0]);
write, format="%s\n", "done.";

write, format="%s", "Attaching astrobj to scenery... ";
sc, astrobj=ao;
write, format="%s\n", "done.";

write, format="%s", "Retrieving astrobj... ";
ao2=sc.astrobj();
write, format="%s\n", "done.";

write, format="%s", "Setting time... ";
noop,sc.screen.time(10000);
write, format="%s\n", "done.";
write, format="%s %e\n", "Checking time:", (time=sc.screen.time());
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
noop,sc.screen.fov(pi/4.);
write, format="%s\n", "done.";
write, format="%s %e\n", "Checking field-of-view:", (fov=sc.screen.fov());
if (fov!=pi/4) error, "CHECK FAILED";

write, format="%s", "Setting resolution... ";
noop,sc.screen.resolution(16);
write, format="%s\n", "done.";
write, format="%s %i\n", "Checking resolution:", (res=sc.screen.resolution());
if (res!=16) error, "CHECK FAILED";


write, format="%s", "Setting inclination... ";
noop,sc.screen.inclination(pi/3.);
write, format="%s\n", "done.";
write, format="%s %e\n", "Checking inclination:",
  (incl=sc.screen.inclination());
if (incl!=pi/3) error, "CHECK FAILED";

if (gyoto_haveXerces()) {
  write, format="%s", "Writing XML description... ";
  sc,xmlwrite="test.xml";
  write, format="%s\n", "done.";
  remove, "test.xml";

  write , format="%s", "Reading Scenery from XML description... ";
  sc3=gyoto_Scenery(GYOTO_EXAMPLES_DIR+"example-moving-star.xml");
 } else {
  write , format="%s", "No Xerces, creating scenery from scratch... ";
  gg  = gyoto_KerrBL();
  spectro = gyoto_Spectrometer("wave", nsamples=1, band=[2e-6, 2.4e-6]);
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
noop, sc3.screen.resolution(32);
noop, sc3.astrobj.radius(2);
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
mask = gyoto_Scenery_rayTrace(sc3);
tac();
write, format="%s\n" , "done.";

write, format="%s", "Setting mask... ";
noop, sc3.screen.mask(mask);
write, format="%s\n" , "done.";

write, format="%s", "Checking mask... ";
mask2 = sc3.screen.mask();
if (!allof(mask2==mask)) error, "CHECK FAILED!";
write, format="%s\n" , "done.";

write, format="%s", "Ray-tracing on 1 thread with mask... \n";
sc3, nthreads=1;
tic;
im1 = gyoto_Scenery_rayTrace(sc3);
tac();
if (!allof(im1==mask)) error, "CHECK FAILED!";
write, format="%s\n" , "done.";

write, format="%s", "Ray-tracing with mask (sc())... \n";
sc3, nthreads=1;
tic;
im1 = sc3(,,"Intensity"); // raytrace
if (!allof(im1==mask)) error, "CHECK FAILED!";
tac();
pli, im1;
pause, 1000;
write, format="%s\n" , "done.";

//noop,sc3.screen(maskwrite="toto.fits");
//noop,sc3.screen(xmlwrite="toto.xml");

/* write, format="%s", "Ray-tracing on adaptive grid... ";
   data = gyoto_Scenery_adaptive_raytrace(sc3, 4);
   fma;
   pli, data(,,3), cmax=100;
   write, format="%s\n" , "done."; */

write, format="%s", "Cloning...";
sc4=sc3.clone();
write, format="%s\n", "DONE.";

write, format="%s\n", "Printing clone:";
sc4;

ph = gyoto_Photon(initcoord=sc3, 6, 19);
ph.is_hit();

doing, "Reading Scenery...";
sc=Scenery(GYOTO_EXAMPLES_DIR+"example-complex-astrobj.xml");
done;

sc, nthreads=8, nprocesses=0, mpispawn=0;

doing, "Integrating whole field...";
data=sc();
done;

r1=8:25:4;
r2=2:-2:3;
v1=[1, 4, 16];
v2=[15, 20, 22];
s1=[[1, 2], [3, 4]];
s2=[[[12, 13], [14, 15]], [[1, 2], [3, 4]], [[10, 11], [16, 17]], [[7, 8], [20, 22]]];

doing, "Integrating subfield...";
data2=sc(r1, r2, );
done;

doing, "Comparing...";
if (anyof(data2 != data(r1, r2, ))) error, "result differ";
done;

doing, "Integrating subfield...";
data2=sc(v1, v2, );
done;

doing, "Comparing...";
if (anyof(data2 != data(v1, v2, ))) error, "result differ";
done;

doing, "Integrating subfield...";
data2=sc(s1, s2, );
done;

doing, "Comparing...";
if (anyof(data2 != data(s1, s2, ))) error, "result differ";
done;

doing, "Integrating subfield...";
data2=sc(r1, v2, );
done;

doing, "Comparing...";
if (anyof(data2 != data(r1, v2, ))) error, "result differ";
done;

doing, "Integrating subfield...";
data2=sc(v1, r2, );
done;

doing, "Comparing...";
if (anyof(data2 != data(v1, r2, ))) error, "result differ";
done;

fov=sc.screen().fov();
npix=sc.screen().resolution();
delta= fov/double(npix);

xx=((indgen(npix)-0.5)*delta-fov/2.)(, -:1:32);
yy=transpose(xx);

data2=array(double, npix, npix, 1);
doing, "integrating pix. by pix., specifying angles...\n";
verbosity=verbose(0);
for (j=1; j<=npix; ++j) {
  write, format="\rj = %d/%d", j, npix;
  for (i=1; i<=npix; ++i) {
    data2(i, j, 1) = sc(-xx(i, j), yy(i, j), );
  }
 }
write, format="%s\n", "";
done;
verbose, verbosity;

doing, "Comparing results";
diff=data-data2;
ind=where(data);
diff(ind)/=data(ind);
mdiff=max(abs(diff));
if (mdiff > 1e-6) error, "Results differ";
output, " OK (max rel. dif.: "+pr1(mdiff)+")";

doing, "integrating whole field, specifying angles...\n";
data2 = sc(-xx, yy, );
done;

doing, "Comparing results";
diff=data-data2;
ind=where(data);
diff(ind)/=data(ind);
mdiff=max(abs(diff));
if (mdiff > 1e-6) error, "Results differ";
output, " OK (max rel. dif.: "+pr1(mdiff)+")";

if (batch()) {

  // Free memory for easier checking with valgrind
  xx=yy=data2=data=ind=[];
  sc4=[];
  sc3=[];
  sc2=[];
  sc=[];
  ao=ao2=[];
  pause, 1000;
  winkill;
 }

end_section, "Scenery";
