/*
    Copyright 2012-2014, 2019 Thibaut Paumard

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

begin_section, "PolishDoughnut Astrobj";

write, format="%s", "Creating default PolishDoughnut... ";
pd = gyoto_PolishDoughnut();
write, format="%s\n", "done.";

write, format="%s\n", "Printing object:";
write, format="%s\n", "-----------------------------";
pd;
write, format="%s\n", "-----------------------------";

// write, format="%s", "Setting lambda to 2... ";
// pd, lambda=2;
// write, format="%s\n", "done.";

// write, format = "%s", "Checking value of lambda... ";
// if (pd(get_lambda=1)==2) write, format="%s\n", "done."; \
//  else error, "CHECK FAILED!";

write, format= "%s" , "Instanciating Kerr metric... ";
gg = gyoto_KerrBL();
write, format= "%s\n", "done.";

write, format= "%s" , "Attaching metric to doughnut... ";
pd, metric=gg;
write, format= "%s\n", "done.";

write, format= "%s" , "Checking attached metric... ";
if (pd.metric()()==gg()) write, format="%s\n", "done."; \
 else error, "CHECK FAILED!";

write, format="%s\n", "Printing object:";
write, format="%s\n", "-----------------------------";
pd;
write, format="%s\n", "-----------------------------";

if (gyoto.haveUDUNITS() && gyoto.haveXerces()) {
  write, format="%s", "Creating PolishDoughnut from file... ";
  pd = gyoto_Scenery(GYOTO_EXAMPLES_DIR+"example-polish-doughnut.xml").astrobj();
  write, format="%s\n", "done.";

  write, format="%s\n", "Printing object:";
  write, format="%s\n", "-----------------------------";
  pd;
  write, format="%s\n", "-----------------------------";

  write, format="%s", "Reading PolishDoughnut scenery... ";
  sc = gyoto_Scenery(GYOTO_EXAMPLES_DIR+"example-polish-doughnut.xml") ;
  write, format="%s\n", "done.";

  write, format="%s", "Setting spectro... ";
  noop, sc.screen.spectro(
                  gyoto_SpectroUniform(kind="freqlog",
                                       nsamples=2,
                                       band=[10, 14],
                                       unit="Hz")
                  );
  noop, sc.astrobj.opticallythin(1);
  write, format="%s\n", "done.";

  write, format="%s", "Ray-tracing scenery... ";
  img = sc(,,"Spectrum")(,,avg);
  write, format="%s\n", "done.";

  write, format="%s", "Displaying image... ";
  fma; pli, img;
  write, format="%s\n", "done.";

  pause, 1000;

  noop, sc.screen.spectro(
                  gyoto_SpectroUniform(kind="freqlog",
                                       nsamples=10,
                                       band=[10, 14],
                                       unit="Hz")
                  );

  write, format="%s", "Integrating one spectrum with radiative transfer...\n";
  s1 = sc(10, 15, "Spectrum[J.m-2.s-1.sr-1.Hz-1]");
  write, format="%s\n", "done.";
  midpoints = sc.screen.spectro.midpoints();
  widths = sc.screen.spectro.widths();
  fma;
  logxy, 1, 1;
  plg, (s1*widths), midpoints, color="red";
  xytitles, "Frequency [Hz]";

  write, format="%s", "Integrating one bin spectrum with radiative transfer...\n";
  s2 = sc(10, 15, "BinSpectrum[J.m-2.s-1.sr-1]");
  write, format="%s\n", "done.";
  channels = sc.screen.spectro.channels();
  widths = sc.screen.spectro.widths();
  s22=array(double,numberof(s2)*2);
  s22(::2)=s2; s22(2::2)=s2;
  chan2=array(double,numberof(s2)*2);
  chan2(::2)=channels(1,); chan2(2::2)=channels(2,);
  plg, s22, chan2;
 }
if (batch()) {
  pause, 1000;
  winkill,0;
 }

end_section, "PolishDoughnut Astrobj";
