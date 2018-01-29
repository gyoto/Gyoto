/**
 * \file GyotoDynamicalDiskBolometric.h
 * \brief A geometrically thin, optically thick disk, evolving dynamically
 *
 *  The disk is described by a set of FITS files for a set of different times
 */

/*
  Copyright 2013, 2018 Frederic Vincent, Thibaut Paumard
  
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

#ifndef __GyotoDynamicalDiskBolometric_H_ 
#define __GyotoDynamicalDiskBolometric_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>

namespace Gyoto{
  namespace Astrobj { class DynamicalDiskBolometric; }
}

//#include <GyotoMetric.h>
#include <GyotoDynamicalDisk.h>

/**
 * \class Gyoto::Astrobj::DynamicalDiskBolometric
 * \brief Geometrically thin disk read from a set of FITS files
 * 
 *   This class describes a PatternDiskBB that evolves dynamically. 
 *   It is described by a set of FITS files.
 *
 */
class Gyoto::Astrobj::DynamicalDiskBolometric 
: public Astrobj::DynamicalDisk {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::DynamicalDiskBolometric>;
 private:
  
  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;

  DynamicalDiskBolometric(); ///< Standard constructor
  
  DynamicalDiskBolometric(const DynamicalDiskBolometric& ) ;///< Copy constructor
  virtual DynamicalDiskBolometric* clone () const; ///< Cloner
  
  virtual ~DynamicalDiskBolometric() ;                        ///< Destructor
  
  double emission(double nu_em, double dsem,
		  state_t const &,
		  double const coord_obj[8]) const;
    
  double bolometricEmission(double dsem, state_t const & cph, double const coord_obj[8]) const;

  void processHitQuantities(Photon* ph, 
			    state_t const &coord_ph_hit,
			    double const *coord_obj_hit, double dt,
			    Properties* data) const;
};

#endif
