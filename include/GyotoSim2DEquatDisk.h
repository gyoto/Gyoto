/**
 * \file GyotoSim2DEquatDisk.h
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


#ifndef __GyotoSim2DEquatDisk_H_ 
#define __GyotoSim2DEquatDisk_H_ 

namespace Gyoto{
  namespace Astrobj { class Sim2DEquatDisk; }
}

#include <GyotoSimBridge.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

class Gyoto::Astrobj::Sim2DEquatDisk :
  public Gyoto::Astrobj::SimBridge {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Sim2DEquatDisk>;

  private:
  double HoverR_; ///< parameter which define the height of the disk depending on radius : h = HoverR * radius.

  public:
  GYOTO_OBJECT; // This object has a (non-inherited) Property list

  Sim2DEquatDisk(); ///< Default constructor
  
  Sim2DEquatDisk(const Sim2DEquatDisk& orig); ///< Copy constructor
  Sim2DEquatDisk * clone() const ;

  ~Sim2DEquatDisk() ;                        ///< Destructor

  virtual std::string className() const ; ///< "Sim2DEquatDisk"
  virtual std::string className_l() const ; ///< "sim2dequatdisk"

  void HoverR(double hh);
  double HoverR() const;

  virtual double operator()(double const coord[4]);


};
#endif