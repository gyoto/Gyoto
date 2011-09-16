/**
 * \file GyotoFixedStar.h
 * \brief Fixed (i.e. non-moving) star
 *
 *  The target of ray-traced Gyoto::Photon
 */

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


#ifndef __GyotoFixedStar_H_ 
#define __GyotoFixedStar_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  class FixedStar;
}

#include <GyotoAstrobj.h>
#include <GyotoMetric.h>

/**
 * \class Gyoto::FixedStar. 
 * \brief Fixed (i.e. non-moving) star
 *
 *  The target of ray-traced Gyoto::Photon
 */
class Gyoto::FixedStar : public Astrobj {
  friend class Gyoto::SmartPointer<Gyoto::FixedStar>;

 // Data : 
  // -----

 protected:
  
  double pos_[3];///< x, y, z or r, theta, phi
  double radius_;///< Star radius
  int use_generic_impact_; ///<Use Astrobj::Impact() or FixedStar::Impact_()
  
  // Constructors - Destructor
  // -------------------------
 public:
  
  /**
   * Everything is undefined, call setCoordSys(), setPos() and
   * setRadius().
   */
  FixedStar();///< Default constructor

  FixedStar(const FixedStar& orig);///< Copy constructor
  virtual FixedStar* clone() const;

  FixedStar(SmartPointer<Gyoto::Metric> gg, double StPsn[3], double radius);
                   ///< Standard constructor
  
  virtual ~FixedStar() ;                        ///< Destructor
  
 public:
  // Accessors
  // ---------
 public:
  double getRadius() const ; ///< Get radius_
  double const * getPos() const; ///< Get const pointer to pos_
  void getPos(double* dst) const; ///< Get a copy of the pos_ array
  //  const int getCoordSys() const; ///< Get coordinate system
  
  virtual void setMetric(SmartPointer<Metric> metric) ;
  void setRadius(double); ///< Set radius
  void setPos(const double[3]); ///< Set pos_ array
  //  void setCoordSys(int); ///< set coordinate system
  void useGenericImpact(int);
  
 public:
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(factoryMessenger *fmp) const ;
  static Astrobj::Subcontractor_t Subcontractor;
  static void Init();
  ///< called from Factory
#endif

  // Outputs
  // -------
 public:
  virtual double operator()(double const coord[4]) ;

  virtual int Impact(Photon *ph, size_t index, AstrobjProperties *data=NULL);


 protected:
  int Impact_(Photon *ph, size_t index, AstrobjProperties *data=NULL);
  virtual void getVelocity(double const pos[4], double vel[4]) ;

  double emission(double nu_em, double dsem, double cp[8], double co[8]=NULL)
    const;


};


#endif
