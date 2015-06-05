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

#include "check-helpers.i"

func check_jacobian(gg, pos) {
  if (is_void(tolerance)) tolerance=1e-6;
  d=dimsof(pos);
  if (d(1)==0) {
    pos=[[0., 10., 1., 2.]];
    d=dimsof(pos);
  } else if (d(1)==1) {
    pos=[ pos ];
    d=dimsof(pos);
  }

  doing, "checking Jacobian", 0;

  for (n=1; n<=d(3); ++n) {
    j1=gg(jacobian=pos(,n));
    j2=array(double, 4, 4, 4);
    for (a=1; a<=4; ++a)
      for (i=1; i<=4; ++i)
        for (j=1; j<=4; ++j)
          j2(j,i,a)=gg(jacobian=pos(,n), j, i, a);
    j3=jacobian(gg, pos(,n));
    j4=jacfunc(gg, pos(,n));
    if (anyof(j1!=j2))
      error, "The two forms of the jacobian method don't yield the same result";
    if (max(abs(j1-j4))>tolerance)
      error, "The Jacobian is wrong";
    if (max(abs(j1-j3))>tolerance)
      error, "The Jacobian doesn't agree with its numerical estimate";
    dot;
  }
 done;
}

func ffunc(gg, pos) {
  spin_=gg.spin(); a2_=spin_^2;
  x=pos(2); y=pos(3); z=pos(4);
  x2=x*x; y2=y*y; z2=z*z; a2z2=a2_*z2;
  tau=x2+y2+z2-a2_;
  rho2=tau*tau+4.*a2z2; rho=sqrt(rho2);
  r2=0.5*(tau+rho);
  r=sqrt(r2); r3=r2*r; r4=r2*r2; r2_a2=r2+a2_;
  rx_ay=r*x+spin_*y; ry_ax=r*y-spin_*x;
  f=2.*r3/(r4+a2_*z2); fr2=f*r2;
  return    f;
}

func dffunc(gg, pos) {
  spin_=gg.spin(); a2_=spin_^2; a4=a2_*a2_;
  x=pos(2); y=pos(3); z=pos(4);
  x2=x*x; y2=y*y; z2=z*z; a2z2=a2_*z2;
  tau=x2+y2+z2-a2_;
  rho2=tau*tau+4.*a2z2; rho=sqrt(rho2);
  r2=0.5*(tau+rho);
  r=sqrt(r2); r3=r2*r; r4=r2*r2; r2_a2=r2+a2_;
  rx_ay=r*x+spin_*y; ry_ax=r*y-spin_*x;
  f=2.*r3/(r4+a2_*z2); fr2=f*r2;

  r4_a2z2=r4+a2z2;
  temp=-(2.*r3*(r4-3.*a2z2))/(r4_a2z2*r4_a2z2*rho);
  x2_y2_z2=x2+y2+z2;
  temp2=(a4+2.*r2*x2_y2_z2 - a2_* (x2_y2_z2 - 4.* z2 + rho));
  return    [
             0.,
	x*temp,	
	y*temp,


-((4.*r*z*(2.* a4*a2_ + (a2_ + 2.*r2)*x2_y2_z2*x2_y2_z2 + 
    a4*(-3.*x2 - 3.*y2 + z2 - 2.*rho) + 
    a2_*(x2 + y2 - z2)*rho))/(rho*temp2*temp2))
             

             ];
}

func dfnum(gg, pos, eps=) {
  grad=array(double, 4);
  if (is_void(eps)) eps=1e-10;
  
  for (i=1; i<=4; ++i) {
    delta=array(0., 4);
    delta(i)=eps;
    grad(i)=(ffunc(gg, pos+delta)-ffunc(gg, pos-delta))/(2.*eps);
  }

  return grad;
}

func kfunc(gg, pos) {
  spin_=gg.spin(); a2_=spin_^2;
  x=pos(2); y=pos(3); z=pos(4);
  x2=x*x; y2=y*y; z2=z*z; a2z2=a2_*z2;
  tau=x2+y2+z2-a2_;
  rho2=tau*tau+4.*a2z2; rho=sqrt(rho2);
  r2=0.5*(tau+rho);
  r=sqrt(r2); r3=r2*r; r4=r2*r2; r2_a2=r2+a2_;
  rx_ay=r*x+spin_*y; ry_ax=r*y-spin_*x;
  f=2.*r3/(r4+a2_*z2); fr2=f*r2;
  return    [
     1.,
     rx_ay/r2_a2,
     ry_ax/r2_a2,
     z/r
     ];
}

func dkfunc(gg, pos) {
  spin_=gg.spin(); a2_=spin_^2;
  x=pos(2); y=pos(3); z=pos(4);
  x2=x*x; y2=y*y; z2=z*z; a2z2=a2_*z2;
  tau=x2+y2+z2-a2_;
  rho2=tau*tau+4.*a2z2; rho=sqrt(rho2);
  r2=0.5*(tau+rho);
  r=sqrt(r2); r3=r2*r; r4=r2*r2; r2_a2=r2+a2_;
  rx_ay=r*x+spin_*y; ry_ax=r*y-spin_*x;
  f=2.*r3/(r4+a2_*z2); fr2=f*r2;

  frac1=1./(r2_a2*r2_a2*rho);
  frac2=z/(r2_a2*r*rho);
  frac3=-z/(r*rho);
  
  return [
	// d/dt
	[0., 0., 0., 0.],
	// d/dx
	[
	  0.,
	  (r3*(x2+rho)-rx_ay*x*(x2+y2+z2+rho)+a2_*(rx_ay*x+r*(x2+rho)))*frac1,
	  (x*(r3*y+a2_*(ry_ax+r*y)-ry_ax*(x2+y2+z2))-(spin_*r2_a2+ry_ax*x)*rho)*frac1,
	  x*frac3
	],
	// d/dy
	[
	  0.,
	  (a2_*(rx_ay+r*x)*y+r2_a2*spin_*rho-y*(-r3*x+rx_ay*(x2+y2+z2+rho)))*frac1,
	  (r3*(y2+rho)-ry_ax*y*(x2+y2+z2+rho)+a2_*(ry_ax*y+r*(y2+rho)))*frac1,
	  y*frac3

	],
	// d/dz
	[
	  0.,
	  ((a2_-r2)*x-2*spin_*r*y)*frac2,
	  ((a2_-r2)*y+2*spin_*r*x)*frac2,
	  (2.*r2- (z2*(a2_ + x2 + y2 + z2 + rho))/rho)/(2.*r3)
	]
      ]
}

func dknum(gg, pos, eps=) {
  grad=array(double, 4, 4);
  if (is_void(eps)) eps=1e-10;
  
  for (i=1; i<=4; ++i) {
    delta=array(0., 4);
    delta(i)=eps;
    grad(,i)=(kfunc(gg, pos+delta)-kfunc(gg, pos-delta))/(2.*eps);
  }

  return grad;
}

func jacfunc(gg, pos) {
  jac=array(double, 4, 4, 4);
  df=dffunc(gg, pos);
  dk=dkfunc(gg, pos);
  k=kfunc(gg, pos);
  f=ffunc(gg, pos);
  for(a=1; a<=4; ++a)
    for (mu=1; mu<=4; ++mu)
      for (nu=1; nu<=mu;++nu)
        jac(nu, mu, a)=jac(mu, nu, a)=df(a)*k(mu)*k(nu)+f*dk(mu, a)*k(nu)+f*k(mu)*dk(nu,a);
  return jac;
}

func christofunc(gg, pos) {
  dst=array(double, 4, 4, 4);
  gup=gg(gmunu_up=pos);
  jac=jacfunc(gg, pos);
  for (a=1; a<=4; ++a) {
    for (mu=1; mu<=4; ++mu) {
      for (nu=1; nu<=4; ++nu) {
	dst(nu, mu, a)=0.;
        for (i=1; i<=4; ++i) {
	  dst(nu, mu, a)+=0.5*gup(a,i)*
	    (jac(nu, i, mu)+jac(i, mu, nu)-jac(nu, mu, i));
	}
      }
    }
  }
  return dst;
}

restore, gyoto;

begin_section, "KerrKS metric", "Kerr in Kerr-Schild coordinates";

doing, "creating KerrKS metric";
gg=KerrKS(spin=0.);
done;

positions=[[0, 10., 12., 5.],
           [0, 5., 2., 7.],
           [0, -10., 0., 50.],
           [0, 0., 0., 10000.]];

check_gmunu, gg, positions;
check_gmunu_up, gg, positions;

doing, "checking dk";
pos=[0, 10., 12., 5.];
dk1=dkfunc(gg, pos);
dk3=dknum(gg, pos);
if (max(abs(dk1-dk3))>1e-5)
  error, "dk is wrong";
done;

doing, "checking df";
pos=[0, 10., 12., 5.];
df1=dffunc(gg, pos);
df3=dfnum(gg, pos);
if (max(abs(df1-df3))>1e-6)
  error, "df is wrong";
done;

check_jacobian, gg, positions;
check_christoffels, gg, positions;

end_section, "KerrKS metric", "Kerr in Kerr-Schild coordinates";
