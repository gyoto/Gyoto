/*
    Copyright 2014 Thibaut Paumard

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
restore, gyoto;

// Put screen on Earth in KerrBL metric
met=KerrBL(mass=4e6, unit="sunmass");
scr=Screen(metric=met, distance=8, unit="kpc");
scr, inclination=pi/2;
scr, time=scr(distance=, unit="kpc"), unit="kpc";

alpha=50.*pi/180/3600/1e6; // 50µas
R0=scr(distance=, unit="geometrical");
Rs=1000.;
tinf=scr(time=, unit="geometrical_time");

//tol=10.^-(indgen(18)+3);
tol=spanl(1e-10, 1e-21, 18);
dmor=spanl(.5, 1e-4, numberof(tol));

// Choose integrator

integrators=["Legacy", "Legacy",
             "runge_kutta_cash_karp54", "runge_kutta_dopri5",
             "runge_kutta_fehlberg78"];
colors=["red", "red", "blue", "magenta", "black"];
types=["solid", "dashdot", "solid", "solid", "solid", "solid"];
dfl=array(double, numberof(tol), numberof(integrators));
ctime=dfl;
niter=long(dfl);

for(k=1; k<=numberof(integrators); ++k) {
  write, format="Trying integrator %s\n", integrators(k); 
  ph=Photon(metric=met, initcoord=scr, alpha, 0.);
  ph, integrator=integrators(k);
  ph, delta=tinf/2;
  ph, maxiter=long(1e10);

  // Two times Legacy: try various numerical parameters
  if (k==1) { legacy_type="SpecificIntegrator";
  } else if (k==2) { legacy_type="GenericIntegrator";
  } else legacy_type=[];
  
  if (k<=2) {
    noop, ph(metric=)(setparameter=legacy_type);
    noop, ph(metric=)(difftol=1e-6);
  }
  
  for (i=1; i<=numberof(tol); ++i) {
    if (k<=2) {
      write, format="Integrator: %s (%s), dm/R: %e\n",
        integrators(k), legacy_type, dmor(i); 
      noop, ph(metric=)(deltamaxoverr=dmor(i));
    } else write, format="Integrator: %s, tol: %e\n", integrators(k), tol(i); 
    verbose, 0; //norm *will* drift
    ph, abstol=tol(i), reltol=tol(i);
    ph, delta=tinf/2.;
    ph, reset=;
    tic;
    vel=ph(get_cartesian=[1, -1]*scr(time=, unit="geometrical_time"))(,4:6);
    ctime(i,k)=tac();
    norm=sqrt((vel^2)(,sum));
    cosine=(vel(1,)*vel(2,))(sum)/(norm(1,)*norm(2,));
    dfl(i,k)=acos(cosine);
    niter(i,k)=dimsof(ph(get_coord=))(2);
  }


 }                                        


mima=[1e-12, 10];
axis1toaxis2=Rs/(Rs+R0)*3600*1e6;

winkill, 0;
window, 0, style="boxed.gs";
// set the viewport to something nice
get_style, land, sys, leg, cleg;
sys.viewport=[0.1757,0.6143,0.509465,0.780535];
set_style, land, sys, leg, cleg;
//
fma;
for (k=3; k<=numberof(integrators); ++k)
  plg, abs(dfl(dif,k))/pi*180*3600*1e6, tol(:-1), color=colors(k), type=types(k), marks=0;
for (k=1; k<=2; ++k)
  plg, abs(dfl(dif,k))/pi*180*3600*1e6, (dmor(:-1))^3/1e8, color=colors(k), type=types(k), marks=0;
logxy, 1, 1;
limits, tol(1), tol(-1), 1e-3, 1e9;
xytitles, "AbsTol and RelTol or DeltaMaxOverR^3^/1e8", "Error on deflection angle [!mas]";
pdf, "doc-fig-dfl-tol.pdf";

winkill, 1;
window, 1, style="boxed.gs";
// set the viewport to something nice
get_style, land, sys, leg, cleg;
sys.viewport=[0.1757,0.6143,0.509465,0.780535];
set_style, land, sys, leg, cleg;
//
fma;
for (k=3; k<=numberof(integrators); ++k)
  plg, ctime(:-1,k), tol(:-1), color=colors(k), type=types(k), marks=0;
for (k=1; k<=2; ++k)
  plg, ctime(:-1,k), (dmor(:-1))^3/1e8, color=colors(k), type=types(k), marks=0, marks=0;
xytitles, "AbsTol and RelTol or DeltaMaxOverR^3^/1e8", "Computing time [s]";
logxy, 1, 1;
limits, tol(1), tol(-1), 0.1, 10;
pdf, "doc-fig-ctime-tol.pdf";

winkill, 2;
window, 2, style="boxed.gs";
// set the viewport to something nice
get_style, land, sys, leg, cleg;
sys.viewport=[0.1757,0.6143,0.509465,0.780535];
set_style, land, sys, leg, cleg;
//
fma;
plsys, 1;
for (k=1; k<=numberof(integrators); ++k) {
  i0=min(where(dfl(dif,k)!=0));
  plg, ctime(i0:-1,k), abs(dfl(dif,k)(i0:))/pi*180*3600e6,
    color=colors(k), type=types(k), marks=0;
 }
xytitles, "Error on deflection angle [!mas]", "Computing time [s]";
logxy, 1, 1;
limits, 1e4, 1e-4, 0.2, 4;
pdf, "doc-fig-dfl-ctime.pdf";

          
bdfx=array(long, numberof(colors));
bdfl=array(double, numberof(colors));
for (k=1; k<=numberof(colors); ++k) {
  a=abs(dfl(dif,k));
  if (anyof (a==0) ) a(where(a==0))=a(max);
  bdfx(k)=a(mnx);
  bdfl(k)=dfl(bdfx(k),k);
 }

dfl_ref=bdfl([1, 3, 4, 5])(avg);

(bdfl-dfl_ref)/pi*180;
