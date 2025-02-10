/**
 * \file GyotoSimThinDisk.h
 *
 */

/*
    Copyright 2024 Irene Urso


    Gyoto is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    Gyoto is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with Gyoto.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef __GyotoSimThinDisk_H_ 
#define __GyotoSimThinDisk_H_ 

namespace Gyoto{
  namespace Astrobj { class SimThinDisk; }
}

#include <GyotoSimBridge.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

class Gyoto::Astrobj::SimThinDisk :
  public Gyoto::Astrobj::SimBridge {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::SimThinDisk>;

  public:
  GYOTO_OBJECT; // This object has a (non-inherited) Property list

  SimThinDisk(); ///< Default constructor
  
  SimThinDisk(const SimThinDisk& orig); ///< Copy constructor
  SimThinDisk * clone() const ;
  
  protected:
  double rin_ ; ///< disk inner radius in geometrical units
  double rout_ ; ///< disk outer radius in geometrical units

  /**
   * Geometrical thickness in geometrical units. Used only in the
   * optically thin regime (flag_radtransf_==1). Should be <<
   * rin_. Default: 1e-3.
   */
  double thickness_; ///< disk thickness
 
  public:
  ~SimThinDisk() ;                        ///< Destructor

  virtual std::string className() const ; ///< "SimThinDisk"
  virtual std::string className_l() const ; ///< "simThindisk"
  
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

  virtual double operator()(double const coord[4]);
  
  virtual double projectedRadius(double const coord[]) const ; ///< Projected radius of position coord on the equatorial plane
  
  virtual int Impact(Gyoto::Photon* ph, size_t index, Astrobj::Properties *data=NULL) ;

};
#endif
