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

// Yorick port of ../tools/RayTracing.C rev. 66
#include "gyoto.i"

/*
  This program aims at computing the null geodesics of photons from an
  observation screen to an astrophysical object (star orbit, fixed
  star or disk).  The final result is a list of illuminated pixels.
  
*/

//****CHOOSE COORDINATES***
ChooseCoord="BL"; //"BL", "KS" or "Kerr"
//*************************

//********************************************
// Definition of Metric and Observer position*
//********************************************

//Metric definition :
mm=1;  //BH mass
aa=0.7;  //BH spin 
write, "\n";
write, format="BH spin= %g\n", aa;

//Kerr metric construction
if (ChooseCoord=="BL") gg=gyoto_KerrBL(spin=aa,mass=mm) ; \
 else if (ChooseCoord=="KS") gg=gyoto_KerrKS(spin=aa,mass=mm) ;
  
//Observer position :

// Observer Boyer Lindquist (BL) coord :
r0=100.;//100;
/* 100 = star orbit ; 45 = disk*/
theta0=1.22;//0.785;//1.22;  //0.05;//1.309;//75 deg //pi/4.;//pi/2.;//x axis
//75 deg = 1.309 rad ; 70 deg = 1.22 rad
//45 deg = 0.785 rad
/* 0.05 = star orbit ; 1.309 = 75 deg = disk */
phi0=0;//0.;
/* 0. = star orbit ; pi/4. = disk */
tobs0=1000.;
delta_t=0.;//31.;//change delta_t if you wan't to make a movie...
write, format="Observing time delta_t= %g\n", delta_t;
tobs=tobs0+delta_t;//observation time
t0=tobs;//for time-reversed integration,  initial time for photon (not for star if ray-tracing a star...)

prec=8;width=15;fmts="%"+pr1(width)+"."+pr1(prec)+"g"; //format string

pos=array(double,8);
  
// Corresponding observer Kerr-Schild space and time coord (cf Hameury, Marck, Pelat 93):
psi0=phi0-0.5*aa/sqrt(mm*mm-aa*aa)*log(abs((1+(r0-mm)/(sqrt(mm*mm-aa*aa)))/(1-(r0-mm)/(sqrt(mm*mm-aa*aa)))));
x0=(r0*cos(psi0)-aa*sin(psi0))*sin(theta0);
y0=(r0*sin(psi0)+aa*cos(psi0))*sin(theta0);
z0=r0*cos(theta0);
T0=1000.;
  
if (ChooseCoord=="BL" || ChooseCoord=="Kerr") pos(1:4)=[t0, r0, theta0, phi0];\
else if (ChooseCoord=="KS") pos(1:4)=[T0, x0, y0, z0];

write, format="%s\n", "setting pos: "+pr1(pos);
screen=gyoto_Screen(metric=gg, observerpos=pos);

//******************
// Object creation *
//******************


{ // variables defined below are local in the C++ code,
  // THIS DOES NOT WORK IN YORICK !!!
  // everything is global
        
  //OBJECT == STAR ORBIT    
  //CI = ISCO in equatorial plane

  Z1=1.+(1.-aa*aa/(mm*mm))^(1./3.)*((1.+aa/mm)^(1./3.)+(1.-aa/mm)^(1./3.));
  Z2=sqrt(3.*aa*aa/(mm*mm)+Z1*Z1);//cf Bardeen 1972 2.21
  risco=mm*(3.+Z2-sqrt((3.-Z1)*(3.+Z1+2.*Z2)));//ISCO
  thetai=1.5708;
  phii=0.;
  rhor=mm+sqrt(mm*mm-aa*aa);//event horizon

  radius=0.5;//0.5;//0.4;//1.2;//NB : cf Trippe et al 2007 ("IR polarized flare") , R_hotspot < 0.3 * R_S = 0.6
  ri=risco+radius;// so that the point of the star closest to the BH is on the ISCO : the entire star must not be closer than ISCO to be able to orbit.
    //pay attention to the fact that changing the radius of the star thus changes its trajectory

  write, format="risco= "+fmts+"\n", risco;
  write, format="rhor= "+fmts+"\n", rhor;
  write, format="star radius= "+fmts+"\n", radius;
  write, format="ri +/- rstar= "+fmts+"  "+fmts+"\n", ri-radius, ri+radius;

  ttravel=r0;//approximate photon travel time (r0=c*tobs=tobs)
  ti=tobs0-4.*ttravel;

  pri=0.;//BL canonical momentum p_r
  pthetai=0.;//BL canonical momentum p_theta
  
  yinit=[ti,ri,thetai,phii,ti,pri,pthetai];
    
  //cf Bardeen 1972 2.12-2.15:
  //these equations assume that the particule mass is 1.
  if (aa==1.) {
    E=(ri+sqrt(mm*ri)-mm)/(ri^(3./4.)*sqrt(ri^(1./2.)+2.*sqrt(mm)));
    L=(mm*(ri^(3./2.)+sqrt(mm)*ri+mm*sqrt(ri)-mm^(3./2.)))/(ri^(3./4.)*sqrt(ri^(1./2.)+2.*sqrt(mm)));
  }else{
    E=(ri^(3./2.)-2*mm*sqrt(ri)+aa*sqrt(mm))/(ri^(3./4.)*sqrt(ri^(3./2.)-3.*mm*sqrt(ri)+2.*aa*sqrt(mm)));
    L=(sqrt(mm)*(ri*ri-2.*aa*sqrt(mm)*sqrt(ri)+aa*aa))/(ri^(3./4.)*sqrt(ri^(3./2.)-3.*mm*sqrt(ri)+2.*aa*sqrt(mm)));
  }
  Q=0.;
  cst=[1.,E,L,Q];//1. = star mass, not used in MakeCoord
    
    
  write,format="dans main E,L= "+fmts+"  "+fmts+"\n", cst(2), cst(3);
    
    
  //Proper time step, number of proper time steps
  //double delta=0.1;
  niter=230;//chosen so that the orbit is computed for [ti;tobs]
    
  init=gg(get_coord = yinit, cst);  
  t0=init(1);
  r0=init(2);
  theta0=init(3);
  //cout << "theta0 dans rt.C : " << theta0 << endl;
  phi0=init(4);
  tdot0=init(5);
  rdot0=init(6);
  thetadot0=init(7);
  phidot0=init(8);
  //cout << "coord init dans rt.C " << t0 << " " << r0 << " " << theta0 << " " << phi0 << " " << tdot0 << " " << rdot0 << " " << thetadot0 << " " << phidot0 << endl;

  //Kerr coordinates : (cf article null hypersurface Eric)
  tdot0_Kerr=tdot0+rdot0/((r0*r0+aa*aa)/(2.*mm*r0)-1.);
  phidot0_Kerr=phidot0+aa*rdot0/(r0*r0-2.*mm*r0+aa*aa);

  v= [rdot0/tdot0, thetadot0/tdot0, phidot0/tdot0];
  pos=[t0, r0,theta0,phi0];

  "gg:";
  gyoto_Metric_g(gg, pos, ChooseCoord);
  "radius: "+pr1(radius);
  "pos: "+pr1(pos);
  "v: "+ pr1(v);
  "ChooseCoord: "+ ChooseCoord;
  
  if (ChooseCoord=="BL")
    orbit = gyoto_Star(metric=gg, radius=radius, initcoord=pos, v) ; \
  else {
    v_Kerr= [rdot0/tdot0_Kerr, thetadot0/tdot0_Kerr, phidot0_Kerr/tdot0_Kerr];
    pos_Kerr=[t0,r0,theta0,phi0];
    //cout << endl;
    //cout << "pos init star= " << t0 << " " << r0 << " " << theta0 << " " << phi0 << endl;
    //cout << endl;

    orbit = gyoto_Star(metric=gg, radius=radius,
                       initcoord=pos_Kerr, v_Kerr) ;
  }
  
  //OBJECT == FIXED STAR
  //NB : theta_ein=sqrt(4.*M*dLS/(dL*dS))
  //   if (ChooseCoord=="BL") {
  //     StarPos=[10.,pi/2.,pi];  //BL star (fixed) space coord
  //   } else if (ChooseCoord=="KS") {
  //     StarPos=[-10.,0.,0.];//KS star (fixed) space coord
  //   }
  //FixedS = gyoto_FixedStar(position=StarPos,radius=radius,metric=gyotoKerrBL());
    
  //OBJECT == DISK
  // ThinDisk = gyoto_ThinInfiniteDisk_new(gg,ChooseCoord);
    
}

//***************************
// Null geodesic integration*
//***************************

//Proper time step, number of proper time steps
deltatau=.01;//1.;

pixels = open("pix.dat", "w");//illuminated pixels on screen

N=1001;//501;//1023;//501;//3001;//number of pixels in each direction
//angular pixel step :
//double delta=pi/(double(N));//all sky view : fov = pi
//warning : vel[3] in Metric::getRayCoord is not defined everywhere.
//double delta=pi/(4.*double(N));//limited view : fov = pi/4. for DISK
delta=pi/(10.*double(N));//limited view : fov = pi/10. for ORBIT
 
ph=gyoto_Photon_new();
coordout=ChooseCoord;

for (i=368;i<=368;i++) {
  //i=251 with a=0.9 integration doesn't progress
    
  write, "\n";
  write, "******************\n";
  write, format= "i = %d\n", i;
  write, "******************\n";
  write, "\n";

  xscr=delta*(i-(N+1)/2.);
  
  for (j=474;j<=474;j++) {
    //j=315 with a=0.9 integration doesn't progress

    yscr=delta*(j-(N+1)/2.);

    coord=screen(raycoord=[-xscr, yscr]);

    //1. Init cond :
    gyoto_Photon_setInitialCondition,ph, gg, orbit, coord;

    //2. delta :
    gyoto_Photon_setDelta,ph,deltatau;

    //3. hit :
    if (gyoto_Photon_hit(ph, 0.))
      write, pixels, format="%15i   %15d   %15.12g\n",
        i, j, gyoto_Star_get_redshift(orbit)^3;

    //4. save :
    //gyoto_Photon_save_xyz(ph,"xyz.dat");
  }
 }

close, pixels;
