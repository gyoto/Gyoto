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

// This file is a port from tools/Orbite.C rev. 66
#include "gyoto.i"
width=15;
prec=12;
mm=1; // BH mass, still unsused
aa=0.995;
radius=1;

gg=gyoto_KerrBL(spin=aa,mass=mm);

//gg_gen=gyoto_AnalyticMetric_new(mm);

//To be changed :
gen=0;// 1=use generic metric, 0=use Kerr
name="StarOrbitXYZ.dat";// name of the saving file
tstop=-2000.;//integration stops when t=tstop

// initial conditions :

//***CI Kerr Equatorial***
ri=10.790954422036;//at apo ; cf Eric Eq. 5.57 et CondInitKerrEquat.nb
r2=ri*ri;
a2=aa*aa;
thetai=pi/2.;
sinu=sin(thetai);
cosi=cos(thetai);
sinu2=sinu*sinu;
cosi2=cosi*cosi;
sigma=r2+a2*cosi2;
phii=0.;
tobs=1000;//observation time (chosen)
robs=100;//observer distance
ttrav=robs;//approximate photon travel time
ti=tobs-4.*ttrav;

//cf Levin & Perez-Giz 2008 Ph Rv D 77, 1003005, Fig. 15 bottom right
EE=.921103;
LL=2.;
QQ=0.;//equatorial plane

//cf LevinPerez Eq. A7
lambda=1.-2*mm*ri/sigma;
xi=2.*mm*aa*ri*sinu2/sigma;
gamma = sinu2*(r2+a2+2.*mm*a2*ri*sinu2/sigma);
fact=1./(gamma*lambda+xi*xi);
phipointi=lambda*fact*LL+xi*fact*EE;
tpointi=-xi*fact*LL+gamma*fact*EE;

write, format= "phipointi %e tpointi %e\n", phipointi, tpointi;

//cf MTW 33.31d :
thetapoint2=(QQ-LL*LL/(tan(thetai)*tan(thetai))+a2*(EE*EE-1.)*cosi2)/(sigma*sigma);
//cout << "thetapoint2= " << thetapoint2 << endl;
if (-thetapoint2 < 1e-10) thetapoint2=0;// if thetapoint2 slightly < 0...
if (thetapoint2 < 0.) error , "thetapoint2 < 0 !!! aborting...";
thetapointi=sqrt(thetapoint2);
//cout << "thetapoint " << thetapointi << endl;
pthetai=sigma*thetapointi;
//cout << "pthetai= " << pthetai << endl;

posi=[ti,ri,thetai,phii];
//cout << "ri thi phi " << setprecision(prec) << ri << " " << thetai << " " << phii << endl;
if (!gen) {
  
  gtt   = gg(posi, 1, 1);
  grr   = gg(posi, 2, 2);
  gthth = gg(posi, 3, 3);
  gphph = gg(posi, 4, 4);
  gtph  = gg(posi, 1, 4);
  
 }else{
  gtt   =gg_gen(posi, 1, 1);
  grr   =gg_gen(posi, 2, 2);
  gthth =gg_gen(posi, 3, 3);
  gphph =gg_gen(posi, 4, 4);
  gtph  =gg_gen(posi, 1, 4);
 }

rpoint2=(-1.-gtt*tpointi*tpointi-2.*gtph*tpointi*phipointi-gthth*thetapointi*thetapointi-gphph*phipointi*phipointi)/grr;
if (-rpoint2 < 1e-10) rpoint2=0;// if rpoint2 slightly < 0...
//cout << "rpoint2 " << rpoint2 << endl;
if (rpoint2 < 0.) error, "rpoint2 < 0 !!! aborting...";

rpointi=sqrt(rpoint2);
Bigdelta=r2-2.*mm*ri+a2;
//cout << "Bigdelta" << Bigdelta << endl;
pri=sigma/Bigdelta*rpointi;
//cout << "pri= " << pri << endl;
  
cst=[1.,EE,LL,QQ];//4 Kerr cst of motion Âµ, E, L, Q  
pti=-EE; pphii=LL;
yinit=[ti,ri,thetai,phii,pti,pri,pthetai,pphii];


if (gen) init=[ti, ri, thetai, phii, tpointi, rpointi, thetapointi, phipointi];\
 else init = gg (get_coord = yinit, cst);

write,"\n";
write, format="phidot= %e tdot= %e\n", init(8), init(5);
write,"\n";

write,"dans Orbit.i yinit= "+pr1(yinit);
write,"dans Orbit.i  init= "+pr1(init);

t0=init(1);
r0=init(2);
theta0=init(3);
phi0=init(4);
tdot0=init(5);
rdot0=init(6);
thetadot0=init(7);
phidot0=init(8);

v= [rdot0/tdot0, thetadot0/tdot0, phidot0/tdot0];

if (gen) st=gyoto_Star(metric=gg_gen,radius=radius,initcoord=init,v) ;\
  else st=gyoto_Star(metric=gg, radius=radius, initcoord=init, v) ;
gyoto_Star_xFill,st,tstop;

gyoto_Star_get_coord,st,t,r;
gyoto_Star_get_xyz, st, x, y, z;
f=open(name, "w");
write, f, format="%15.12g  %15.12g  %15.12g  %15.12g\n", r, x, y, z;
close, f;
 
