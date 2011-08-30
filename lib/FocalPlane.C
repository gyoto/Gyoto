/*
    Copyright 2011 Frederic Vincent, Thibaut Paumard

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

#include <GyotoFocalPlane.h>
using namespace Gyoto;

FocalPlane::FocalPlane(SmartPointer<Metric> gg, SmartPointer<Astrobj> obj,
		       double xmin, double ymin,
		       size_t nx, size_t ny, double dx, double dy):
  gg_(gg), obj_(obj),
  nx_(nx), ny_(ny), dx_(dx), dy_(dy),
  xmin_(xmin), ymin_(ymin)
{

}

FocalPlane::FocalPlane(SmartPointer<Metric> gg, SmartPointer<Astrobj> obj,
		       double xmin, double ymin, double xmax, double ymax,
		       size_t nx, size_t ny):
  gg_(gg), obj_(obj),
  nx_(nx), ny_(ny), dx_((xmax-xmin)/nx), dy_((ymax-ymin)/ny),
  xmin_(xmin), ymin_(ymin)
{

}

FocalPlane::~FocalPlane() {}

size_t FocalPlane::getNx() const { return nx_; }
size_t FocalPlane::getNy() const { return ny_; }
double FocalPlane::getDx() const { return dx_; }
double FocalPlane::getDy() const { return dy_; }
double FocalPlane::getXmin() const { return xmin_; }
double FocalPlane::getYmin() const { return ymin_; }
double FocalPlane::getXmax() const { return xmin_+(nx_-1)*dx_; }
double FocalPlane::getYmax() const { return ymin_+(ny_-1)*dy_; }
void FocalPlane::getX(double * dest) const {
  size_t i; for (i=0; i<nx_; ++i) dest[i]=xmin_+i*dx_;
}
void FocalPlane::getY(double * dest) const {
  size_t i; for (i=0; i<ny_; ++i) dest[i]=ymin_+i*dy_;
}

// void FocalPlane::hitMap(double tobs, double *dest, double deltatau) {
//   size_t i, j;
//   double x, y;
//   double coord[8];
//   int sys;
//   Photon ph;


//   // en substance
//   double mindate=gg_ -> getDistance() * N * obj_->getRmax();

//   for (j=0; j<ny; ++j) {
//     y=ymin_+j*dy_;
//     for (i=0; i<nx; ++i) {
//       x=xmin_+j*dx_;
//       // get coord and sys
//       gg->getObserverCoord(tobs, x, y, coord, sys);
//       ph.setInitialCondition(gg, obj, coord, deltatau, sys);
//       dest[j*nx+i]=ph.hit(obj->getTmin());
//     }
//   }
// }
