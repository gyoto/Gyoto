/**
 * \file GyotoThinDisk.h
 * \brief Geometrically thin disks and rings
 *
 *  Gyoto::Astrobj::ThinDisk is a class to represent geometrically
 *  thin, optically thick or thin disks or rings in the equatorial
 *  plane of the object. It therefore assumes the metric has an
 *  equatorial plane, which orresponds to z==0 in a Cartesian
 *  coordinate system or to theta==M_PI/2 in a sperical coordinate
 *  system.
 *
 *  This calls is not abstract and can be used as is (it keeps the
 *  very simplistic Generic::emission() and Generci::transmission()),
 *  but it is also a base class to develop classes with more complex
 *  emission laws.
 *
 */

/*
    Copyright 2011-2015 Thibaut Paumard, Frederic Vincent

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


#ifndef __GyotoThinDisk_H_ 
#define __GyotoThinDisk_H_ 

namespace Gyoto{
  namespace Astrobj { class ThinDisk; }
}

#include <GyotoMetric.h>
#include <GyotoAstrobj.h>
#include <GyotoSpectrum.h>
#include <GyotoFunctors.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

/**
 * \class Gyoto::Astrobj::ThinDisk
 * \brief Geometrically thin disks and rings
\code
   <InnerRadius> rin_ </InnerRadius>
   <OuterRadius> rout_ </OuterRadius>
   <CounterRotating/>
\endcode
 * ThinDisk::setParameter() also takes care of calling
 * Generic::setParameter().
 */
class Gyoto::Astrobj::ThinDisk :
  public Gyoto::Astrobj::Generic,
  public Gyoto::Functor::Double_constDoubleArray
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::ThinDisk>;
  
  // Data : 
  // -----
 protected:
  double rin_ ; ///< disk inner radius in geometrical units
  double rout_ ; ///< disk outer radius in geometrical units

  /**
   * Geometrical thickness in geometrical units. Used only in the
   * optically thin regime (flag_radtransf_==1). Should be <<
   * rin_. Default: 1e-3.
   */
  double thickness_; ///< disk thickness
  int dir_; ///< 1 for corotating (default), -1 for counterrotating.
  unsigned int velocitykind_; ///< tag for VelocityKind

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;

  /**
   * Create direct ThinDisk object. When initializing a derived class,
   * always set kind.
   */
  ThinDisk(std::string kind="ThinDisk"); ///< Default constructor
  
  ThinDisk(const ThinDisk& orig); ///< Copy constructor
  virtual ThinDisk* clone () const ;

  virtual ~ThinDisk() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  virtual double innerRadius() const ; ///< Get rin_
  virtual double innerRadius(std::string const &) const ; ///< Get rin_
  virtual void   innerRadius(double); ///< Set rin_
  virtual void   innerRadius(double, std::string const &); ///< Set rin_
  virtual double outerRadius() const ; ///< Get rout_
  virtual double outerRadius(std::string const &) const ; ///< Get rout_
  virtual void   outerRadius(double); ///< Set rout_
  virtual void   outerRadius(double, std::string const &); ///< Set rout_
  virtual double thickness() const ; ///< Get thickness_
  virtual double thickness(std::string const &) const ; ///< Get thickness_
  virtual void   thickness(double); ///< Set thickness_
  virtual void   thickness(double, std::string const&); ///< Set thickness_
  virtual int    dir() const ; ///< Get dir_
  virtual void   dir(int); ///< Set dir_
  virtual bool   corotating() const; /// Get dir_==1
  virtual void   corotating(bool t); /// Set dir_=t?1:-1
  virtual std::string velocityKind() const ; ///< Get VelocityKind
  virtual void   velocityKind(std::string const&); ///< Set VelocityKind

  /**
   * A function which changes sign on the equatorial plane.
   */
  virtual double operator()(double const coord[]) ; ///< theta-pi/2 or z

  virtual double projectedRadius(double const coord[]) const ;
      ///< Projected radius of position coord on the equatorial plane

  virtual double sphericalPhi(double const coord[]) const;
      ///< Longitude

  /// Get fluid 4-velocity at point.
  /**
   * Fill vel with the 4-vector velocity of the fluid at 4-position
   * pos. getVelocity() should work at some distance from the
   * equatorial plane. The default implementation calls
   * Metric::Generic::circularVelocity().
   *
   * \param[in] pos 4-position at which to compute velocity;
   * \param[out] vel 4-velocity at pos.
   */
  virtual void getVelocity(double const pos[4], double vel[4])  ;

  public:
  virtual int Impact(Gyoto::Photon* ph, size_t index,
		     Astrobj::Properties *data=NULL) ;

};


#endif
