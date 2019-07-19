/**
 * \file GyotoStarTrace.h
 * \brief Like a Star that would be on all points of its orbit at all time
 *
 * A StarTrace is a Star that is considerred to be simultaneously on
 * all the points of its orbit at all time. The purpose is to
 * precompute quickly an integrated image that can later be used as a
 * mask to efficiently compute many images of the underlying Star at
 * varying observing dates.
 */

/*
    Copyright 2013-2015 Thibaut Paumard

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
 * StarTrace inherits all the members and methods from Star. It has
 * two additional members, tmin_ and tmax_, which specify the time
 * interval of the Star's orbit that is to be considerred illuminated.
 *
 * A StarTrace is not (necessarily) continuous: the Star is
 * considerred to be present at all the locations computed by xFill(),
 * meaning that if the integration step is large compared to radius_,
 * the object will be a collection of discrete blobs. To ensure
 * continuity, one should use a non-adaptive step and specify a
 * reasonable step. Computation is also faster in optically thick
 * mode.
 *
 * \code
 * <Astrobj kind="StarTrace">
 *   ...
 *   <TMin> 600 </TMin>
 *   <TMax> 600 </TMax>
 *   <NonAdaptive/>
 *   <Radius> 1 </Radius>
 *   <Delta> 1 </Delta>
 * </Astrobj>
 * \endcode
 *
 */
class Gyoto::Astrobj::StarTrace :
  public Gyoto::Astrobj::Star {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::StarTrace>;
  
  // Data : 
  // -----
  protected:
  double tmin_; ///< Minimum date to consider on the underlying Star orbit
  double tmax_; ///< Maximum date to consider on the underlying Star orbit
  double * x_; ///< Cartesian x
  double * y_; ///< Cartesian y
  double * z_; ///< Cartesian z

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;

 /**
  * \brief Create Star object and set initial condition.
  *
  * \param gg Gyoto::SmartPointer to the Gyoto::Metric in this part of the Universe
  * \param radius star radius
  * \param pos initial 4-position
  * \param v   initial 3-velocity
  */
  StarTrace(SmartPointer<Metric::Generic> gg, double radius,
       double const pos[4], double const v[3]) ;

 /**
  * Create Star object with undefined initial conditions. One needs to
  * set the coordinate system, the metric, and the initial position
  * and velocity before integrating the orbit. setInititialCondition()
  * can be used for that.
  */
  StarTrace(); ///< Default constructor
  
  StarTrace(const StarTrace& orig); ///< Copy constructor


  /// Build StarTrace from Star
  StarTrace(const Star& o, double tmin, double tmax);

  virtual StarTrace * clone() const ;

  virtual ~StarTrace() ;                        ///< Destructor
  
  using Star::xAllocate;
  void xAllocate(size_t);
  void xAllocateXYZ(); ///< Allocate x_, y_, z_
  using Star::xExpand;
  size_t xExpand(int);

  void computeXYZ(size_t i); ///< Compute (and cache) x_, y_ and z_ for one date
  void computeXYZ(); ///< Compute (and cache) x_, y_ and z_

  using Star::setInitCoord;
  virtual void setInitCoord(const double coord[8], int dir = 0);

  using Generic::metric;
  virtual void metric(SmartPointer<Metric::Generic> gg);

  virtual void xStore(size_t ind, state_t const &coord, double tau) ;

  // Accessors
  // ---------
 public:
  virtual std::string className() const ; ///< "StarTrace"
  virtual std::string className_l() const ; ///< "startrace"

  double TMin()const; ///< Get tmin_
  void TMin(double); ///< Set tmin_
  double TMax()const; ///< Get tmax_
  void TMax(double); ///< Set tmax_

  using Star::setInitialCondition;
  virtual void setInitialCondition(double const coord[8]); ///< Same as Worldline::setInitialCondition(gg, coord, sys,1)

  virtual double operator()(double const coord[4]) ;

};


#endif
