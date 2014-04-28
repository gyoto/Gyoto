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
#include "gyoto_std.i"

// There is no specific Yorick interface for the Minkowski metric.
// The only limitation is that is is not possible to retrieve the kind
// of coordinate system in use.

func doing(msg) { write, format="%s... ", msg; }
func done {write, format="%s\n", "done.";}
func output(msg) {write, format="%s\n", msg;}

func christoffel(met, pos, eps=) {
  // works only for diagonal metric
  res=array(double, 4, 4, 4);
  grad=array(double, 4, 4, 4);
  if (is_void(eps)) eps=1e-10;
  
  g0=met(pos);

  for (i=1; i<=4; ++i) g0(i, i)=1./g0(i, i);
  
  for (i=1; i<=4; ++i) {
    delta=array(0., 4);
    delta(i)=eps;
    grad(,,i)=(met(pos+delta)-met(pos-delta))/(2.*eps);
  }

  for (a=1; a<=4; ++a) {
    for (mu=1; mu<=4; ++mu) {
      for (nu=1; nu<=4; ++nu) {
        for (i=1; i<=4; ++i) {
          res(nu, mu, a) += 0.5*g0(a, i)*(grad(nu,i,mu)+grad(i,mu,nu)-grad(nu, mu, i));
        }
      }
    }
  }
  return res;
  
}



restore, gyoto;

output, "**************************************************";
output, "CHECKING MINKOWSKI METRIC IN CARTESIAN COORDINATES";
output, "**************************************************";
  

doing, "creating Minkowski metric";
gg=Metric("Minkowski");
done;

doing, "creating Star in this metric";
st = Star(metric=gg, initcoord=[0., 0., 0., 0.], [0., 0., 0.]);
done;

doing, "integrating geodesic";
st, xfill=1000.;
done;

doing, "checking results";
data=st(get_coord=);
d=dimsof(data);
write, format="(integration produced %d rows)...", d(2);
if (!(allof(data(,2:)==([0., 0, 0, 1., 0, 0, 0](-:1:d(2), )))))
  error, "integration produced wrong results";
done;

doing, "checking another star";
pos=[0., 4., 2., 6.];
vel=[0.5, 0.2, 0.4];
st = Star(metric=gg, initcoord=pos, vel);
st, xfill=1000.;
data=st(get_coord=);
maxerr=max(abs(minmax(data(, 2:4)-data(,1)(,-)*vel(-,)-pos(-,2:4))));
if (maxerr>1e-10) error, "integration produced wrong results";
done;

doing, "checking a photon";
pos=[0., 4., 2., 6.];
vel=[0.5, 0.2, 0.4];
st = Photon(metric=gg, initcoord=pos, vel);
st, xfill=1000.;
data=st(get_coord=);

nrows=dimsof(data)(2);
norm=array(double, nrows);
for(i=1; i<=nrows; ++i)
  norm(i)=gg(scalarprod=data(i, 1:4), data(i, 5:8), data(i, 5:8));
if (max(norm) > 1e-10) "error: norm was not conserved";

maxerr=max(abs(minmax(data(, 2:4)-
                      data(,1)(,-)*vel(-,)/data(,5)(,-)-pos(-,2:4))));
if (maxerr>1e-10) error, "integration produced wrong results";
done;

output, "**************************************************";
output, "CHECKING MINKOWSKI METRIC IN SPHERICAL COORDINATES";
output, "**************************************************";

doing, "changing coordinate system";
gg, setparameter="Spherical";
done;

doing, "checking christoffels";
pos=[0, 10., pi/2, 0.];
Gamma1=gg(christoffel=pos);
Gamma2=array(double, 4, 4, 4);
for (a=1; a<=4; ++a)
  for (i=1; i<=4; ++i)
    for (j=1; j<=4; ++j)
      Gamma2(j,i,a)=gg(christoffel=pos, j, i, a);
Gamma3=christoffel(gg, pos);
if (anyof(Gamma1!=Gamma2))
  error, "The two forms of the christoffel method don't yield the same result";
if (max(abs(Gamma1-Gamma2))>1e-6)
  error, "The cristoffels don't agree with their numerical estimate";
done;

doing, "checking motionless star";
pos=[0., 4., 2., 6.];
vel=[0., 0., 0.];
st = Star(metric=gg, initcoord=pos, vel);
st, xfill=1000.;
data=st(get_coord=);
d=dimsof(data);
write, format="(integration produced %d rows)...", d(2);
if (!(allof(data(,2:)==(_(pos(2:),[1.],vel)(-:1:d(2), )))))
  error, "integration produced wrong results";
done;

doing, "checking moving star";
pos=[0., 10.791, pi/2., 0];
vel=[0., 0., 0.016664];
st = Star(metric=gg, initcoord=pos, vel);
st, adaptive=1;
st, delta=0.1;
st, xfill=1000.;
data=st(get_coord=);

dates=data(,1);
txyz=st(get_cartesian=dates);

d=dimsof(data);
write, format="(integration produced %d rows)...", d(2);
//if (!(allof(data(,2:)==(_(pos(2:),[1.],vel)(-:1:d(2), )))))
//  error, "integration produced wrong results";
done;



write, format="\n%s\n", "ALL TESTS PASSED";
