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

  // Constructors - Destructor
  // -------------------------
 public:

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
  virtual double getInnerRadius() const ; ///< Get rin_
  virtual double getInnerRadius(std::string unit) const ; ///< Get rin_
  virtual void   setInnerRadius(double); ///< Set rin_
  virtual void   setInnerRadius(double, std::string unit); ///< Set rin_
  virtual double getOuterRadius() const ; ///< Get rout_
  virtual double getOuterRadius(std::string unit) const ; ///< Get rout_
  virtual void   setOuterRadius(double); ///< Set rout_
  virtual void   setOuterRadius(double, std::string unit); ///< Set rout_
  virtual double getThickness() const ; ///< Get thickness_
  virtual double getThickness(std::string unit) const ; ///< Get thickness_
  virtual void   setThickness(double); ///< Set thickness_
  virtual void   setThickness(double, std::string unit); ///< Set thickness_
  virtual int    getDir() const ; ///< Get dir_
  virtual void   setDir(int); ///< Set dir_

  /**
   * A function which changes sign on the equatorial plane.
   */
  virtual double operator()(double const coord[]) ; ///< theta-pi/2 or z

  virtual double projectedRadius(double const coord[]) const ;
      ///< Projected radius of position coord on the equatorial plane

  virtual double sphericalPhi(double const coord[]) const;
      ///< Longitude

  /**
   * Used by Standard::Impact().
   *
   * Fill vel with the 4-vector velocity of the fluid at 4-position
   * pos. getVelocity() should work at some distance from the
   * equatorial plane. The default implementation calls
   * Metric::Generic::circularVelocity().
   *
   * \param pos input, 4-position at which to compute velocity;
   * \param vel output, 4-velocity at pos.
   */
  virtual void getVelocity(double const pos[4], double vel[4])  ;

 public:
  virtual int setParameter(std::string name,
			   std::string content,
			   std::string unit) ;

#ifdef GYOTO_USE_XERCES
  /**
   * The sub-classes implementations of the
   * Astrobj::Generic::fillElement() method should call
   * Astrobj::ThinDisk::fillElement() to fill the common bits.
   */
  virtual void fillElement(FactoryMessenger *fmp) const ;
  ///< Fill the generic XML bits
#endif

  public:
  virtual int Impact(Gyoto::Photon* ph, size_t index,
		     Astrobj::Properties *data=NULL) ;

};


#endif
