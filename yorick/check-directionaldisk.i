/*
    Copyright 2014 Frederic Vincent, Thibaut Paumard

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

begin_section, "DirectionalDisk Astrobj";

nnu=2; ni=2; nr=10;
intensity=array(double, nnu, ni, nr);
intensity(,1,1::2)=10.;
intensity(,1,2::2)=7.;
intensity(,1,3::8)=4.;
intensity(,1,9::10)=1.;
intensity(,2,)=0.1*intensity(,1,);

freq=array(double,nnu);
cosi=array(double,ni);
radius=indgen(nr)*5.;

freqobs=1e18;
freq(1)+=freqobs*100.; freq(2)+=freqobs/100.; // in decreasing order
//cosi+=0.5;
cosi(1)=0.3;cosi(2)=0.6;

metric = gyoto_KerrBL(mass=4e6*GYOTO_SUN_MASS);

write, format="%s", "Creating DirectionalDisk...";
pd = gyoto_DirectionalDisk(copyintensity=intensity,
                           copygridfreq=freq, copygridcosi=cosi,
                           copygridradius=radius,
                           metric=metric, rmax=100,
                           innerradius=6, outerradius=50
                           );
write, format="%s\n", " done.";

write, format="%s\n", "Printing DirectionalDisk:";
pd;
write, format="%s\n", " done.";

screen = gyoto_Screen(metric=metric, resolution=64,
                      time=1000.*metric.unitlength()/GYOTO_C,
                      distance=100.*metric.unitlength(), fov=1.1,
                      inclination=110./180.*pi, paln=pi,
                      freqobs=freqobs);

write, format="%s", "Attaching DirectionalDisk to scenery...";
sc = gyoto_Scenery(metric=metric, screen=screen, astrobj=pd);
write, format="%s\n", " done.";

if (gyoto_haveXerces() && gyoto_haveCFITSIO()) {
  write, format="%s", "Saving data to fits file...";
  pd, fitswrite="!check-directionaldisk.fits.gz";
  write, format="%s\n", " done.";

  write, format="%s", "Saving scenery to XML file...";
  sc, xmlwrite="check-directionaldisk.xml";
  write, format="%s\n", " done.";

  write, format="%s", "Reading back scenery...";
  sc2 = gyoto_Scenery("check-directionaldisk.xml");
  write, format="%s\n", " done.";

  write, format="%s", "Removing temporary files...";
  remove, "check-directionaldisk.xml";
  remove, "check-directionaldisk.fits.gz";
  write, format="%s\n", " done.";
 } else {
  write, format="%s", "Cloning...";
  sc2 = sc.clone();
  write, format="%s\n", " done.";
 }
  
write, format="%s", "Getting DirectionalDisk...";
pd2 = sc2.astrobj();
write, format="%s\n", " done.";

write, format="%s", "Comparing intensity array...";
if (anyof(intensity != pd2.copyintensity())) error, "CHECK FAILED";
write, format="%s\n", " done.";

write, format="%s", "Comparing freq array...";
if (anyof(freq != pd2.copygridfreq())) error, "CHECK FAILED";
write, format="%s\n", " done.";

write, format="%s", "Comparing cosi array...";
if (anyof(cosi != pd2.copygridcosi())) error, "CHECK FAILED";
write, format="%s\n", " done.";

write, format="%s", "Comparing radius array...";
if (anyof(radius != pd2.copygridradius())) error, "CHECK FAILED";
write, format="%s\n", " done.";

write, format="%s", "Performing raytracing...\n";
im = sc();
write, format="%s\n", "done.";

write, format="%s", "Displaying image...";
window, style="nobox.gs";
pli, im;
write, format="%s\n", " done.";
pause, 1000;
if (batch()) winkill;

end_section, "DirectionalDisk Astrobj";
