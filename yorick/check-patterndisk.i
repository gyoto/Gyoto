/*
    Copyright 2011-2014, 2019 Thibaut Paumard

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

begin_section, "PatternDisk Astrobj";

opacity=array(double, 11, 3, 1);
opacity(1::2, 1::2, )=100.;
opacity(2::2, 2::2, )=100.;

intensity=opacity*0.+1.;

metric = gyoto_KerrBL(mass=4e6*GYOTO_SUN_MASS);

write, format="%s", "Creating PatternDisk...";
pd = gyoto_PatternDisk(copyintensity=intensity, copyopacity=opacity,
                       innerradius=3, outerradius=28, repeatphi=8,
                       metric=metric, rmax=50);
write, format="%s\n", " done.";

write, format="%s\n", "Printing PatternDisk:";
pd;
write, format="%s\n", " done.";

screen = gyoto_Screen(metric=metric, resolution=64,
                      time=1000.*metric.unitlength()/GYOTO_C,
                      distance=100.*metric.unitlength(), fov=30./100.,
                      inclination=110./180.*pi, paln=pi);

write, format="%s", "Attaching PatternDisk to scenery...";
sc = gyoto_Scenery(metric=metric, screen=screen, astrobj=pd);
write, format="%s\n", " done.";

if (gyoto_haveXerces() && gyoto_haveCFITSIO()) {
  write, format="%s", "Saving data to fits file...";
  pd, fitswrite="!check-patterndisk.fits.gz";
  write, format="%s\n", " done.";

  write, format="%s", "Saving scenery to XML file...";
  sc, xmlwrite="check-patterndisk.xml";
  write, format="%s\n", " done.";

  write, format="%s", "Reading back scenery...";
  sc2 = gyoto_Scenery("check-patterndisk.xml");
  write, format="%s\n", " done.";

  doing, "Checking dmax";
  if (sc2.screen().dmax() != sc.screen().dmax())
    error, "dmax was not conserved when writing and reading XML";
  done;

  doing, "Checking tmin";
  if (sc2.tmin() != sc.tmin())
    error, "tmin was not conserved when writing and reading XML";
  done;

  write, format="%s", "Removing temporary files...";
  remove, "check-patterndisk.xml";
  remove, "check-patterndisk.fits.gz";
  write, format="%s\n", " done.";
 } else {
  write, format="%s", "Cloning...";
  sc2 = sc.clone();
  write, format="%s\n", " done.";
 }
  
write, format="%s", "Getting PatternDisk...";
pd2 = sc2.astrobj();
write, format="%s\n", " done.";

write, format="%s", "Comparing intensity array...";
if (anyof(intensity != pd2.copyintensity())) error, "CHECK FAILED";
write, format="%s\n", " done.";

write, format="%s", "Comparing opacity array...";
if (anyof(opacity != pd2.copyopacity())) error, "CHECK FAILED";
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

end_section, "PatternDisk Astrobj";
