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

// There is no specific Yorick interface for the Minkowski metric.
// The only limitation is that is is not possible to retrieve the kind
// of coordinate system in use.

restore, gyoto;

begin_section, "Minkowski metric", "in cartesian coordinates";

positions=[[0, 10., 12., 5.],
           [0, 5., 2., 7.],
           [0, -10., 0., 50.],
           [0, 0., 0., 10000.]];

doing, "creating Minkowski metric";
gg=Metric("Minkowski");
done;

check_gmunu, gg, positions;
check_christoffels, gg, positions;

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

begin_section, "Minkowski metric", "in spherical coordinates";

doing, "changing coordinate system";
gg, setparameter="Spherical";
done;

positions=[[0, 10., pi/2., 0.],
           [100., 10., pi/4., 2.],
           [1000., 100., 1., 1.]];

check_gmunu, gg, positions;
check_christoffels, gg, positions;

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

end_section, "Minkowski metric";
