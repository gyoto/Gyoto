/**
 * \file GyotoStarTrace.h
 * \brief Like a Star that would be on all points of its orbit at all time
 *
 */

/*
    Copyright 2013 Thibaut Paumard

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


#ifndef __GyotoStarTrace_H_ 
#define __GyotoStarTrace_H_ 

namespace Gyoto{
  namespace Astrobj { class StarTrace; }
}

#include <GyotoStar.h>

/**
 * \class Gyoto::Astrobj::StarTrace
 * \brief Like a Star that would be on all points of its orbit at all time
 *
 */
class Gyoto::Astrobj::StarTrace :
  public Gyoto::Astrobj::Star {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::StarTrace>;
  
  // Data : 
  // -----
  protected:
  double tmin_;
  double tmax_;
  double * x_; ///< Cartesian x
  double * y_; ///< Cartesian y
  double * z_; ///< Cartesian z

  // Constructors - Destructor
  // -------------------------
 public:
 /**
  * Create Star object and set initial condition.
  * \param gg Gyoto::SmartPointer to the Gyoto::Metric in this part of the Universe
  * \param radius star radius
  * \param pos initial 4-position
  * \param v   initial 3-velocity
  */
  StarTrace(SmartPointer<Metric::Generic> gg, double radius,
       double pos[4], double v[3]) ;                        ///< Standard constructor

 /**
  * Create Star object with undefined initial conditions. One needs to
  * set the coordinate system, the metric, and the initial position
  * and velocity before integrating the orbit. setInititialCondition()
  * can be used for that.
  */
  StarTrace(); ///< Default constructor
  
  StarTrace(const StarTrace& orig); ///< Copy constructor
  virtual StarTrace * clone() const ;

  virtual ~StarTrace() ;                        ///< Destructor
  
  using Star::xAllocate;
  void xAllocate(size_t);
  void xAllocateXYZ(); ///< Allocate x_, y_, z_
  using Star::xExpand;
  size_t xExpand(int);

  void computeXYZ(size_t i);
  void computeXYZ();

  using Star::setInitCoord;
  virtual void setInitCoord(const double coord[8], int dir);

  virtual void setMetric(SmartPointer<Metric::Generic> gg);

  virtual void xStore(size_t ind, double coord[8]) ;

  // Accessors
  // ---------
 public:
  virtual std::string className() const ; ///< "StarTrace"
  virtual std::string className_l() const ; ///< "startrace"

  double TMin();
  void TMin(double);
  double TMax();
  void TMax(double);

  virtual void setInitialCondition(double coord[8]); ///< Same as Worldline::setInitialCondition(gg, coord, sys,1)

  virtual int setParameter(std::string name,
			   std::string content,
			   std::string unit);

#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp) const ; ///< called from Factory
#endif

  virtual double operator()(double const coord[4]) ;

};


#endif
