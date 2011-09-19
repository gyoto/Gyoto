/**
 * \file GyotoFocalPlane.h
 * \brief Observed image
 *
 *  A bunch of Gyoto::Photon instances
 */
/*
    Copyright 2011 Thibaut Paumard, Frederic Vincent

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

#ifndef __GyotoFocalPlane_H_ 
#define __GyotoFocalPlane_H_ 

namespace Gyoto {
  class FocalPlane;
}

#include <GyotoMetric.h>
#include <GyotoPhoton.h>
#include <GyotoAstrobj.h>

/**
 * \classe Gyoto::FocalPlane
 * \brief Observed image
 *
 * A grid of Gyoto::Photon instances
 */
class Gyoto::FocalPlane {
 protected:
  SmartPointer<Metric::Generic> gg_;
  SmartPointer<Astrobj::Generic> obj_;
  size_t nx_;
  size_t ny_;
  double dx_; ///< radians
  double dy_; ///< radians
  double xmin_; ///< radians
  double ymin_; ///< radians
 public:
  FocalPlane(SmartPointer<Metric::Generic> gg, SmartPointer<Astrobj::Generic> obj,
	     double xmin, double ymin,
	     size_t nx, size_t ny, double dx, double dy);
  FocalPlane(SmartPointer<Metric::Generic> gg, SmartPointer<Astrobj::Generic> obj,
	     double xmin, double ymin, double xmax, double ymax,
	     size_t nx, size_t ny);
  ~FocalPlane();

  /**
   * \param double dest[nx][ny]
   */
  void hitMap(double *dest);
  size_t getNx() const;
  size_t getNy() const;
  double getDx() const;
  double getDy() const;
  double getXmin() const;
  double getYmin() const;
  double getXmax() const;
  double getYmax() const;

  /**
   * \param double dest[] : an nx-length array of doubles which will
   * hold the x values upon completion.
   */
  void getX(double * dest) const; ///< Get X axis

  /**
   * \param double dest[] : an ny-length array of doubles which will
   * hold the x values upon completion.
   */
  void getY(double * dest) const; ///< Get Y axis

  //  void hitMap(double tobs, double *dest, double deltatau=0.1) ;

};

#endif
