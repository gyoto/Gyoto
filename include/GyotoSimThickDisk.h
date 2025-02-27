/**
 * \file GyotoSimThickDisk.h
 *
 */

/*
    Copyright 2024 Nicolas Aimar


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


#ifndef __GyotoSimThickDisk_H_ 
#define __GyotoSimThickDisk_H_ 

namespace Gyoto{
  namespace Astrobj { class SimThickDisk; }
}

#include <GyotoSimBridge.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

class Gyoto::Astrobj::SimThickDisk :
  public Gyoto::Astrobj::SimBridge {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::SimThickDisk>;

  private:
  double HoverR_; ///< parameter which define the height of the disk depending on radius : h = HoverR * radius.

  public:
  GYOTO_OBJECT; // This object has a (non-inherited) Property list

  SimThickDisk(); ///< Default constructor
  
  SimThickDisk(const SimThickDisk& orig); ///< Copy constructor
  SimThickDisk * clone() const ;

  ~SimThickDisk() ;                        ///< Destructor

  virtual std::string className() const ; ///< "SimThickDisk"
  virtual std::string className_l() const ; ///< "SimThickDisk"

  void HoverR(double hh);
  double HoverR() const;

  virtual double operator()(double const coord[4]);

  void filename(std::string const &d); ///< Overload of the SimBridge function to set the default HoverR from FITS files


};
#endif