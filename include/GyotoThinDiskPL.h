/**
 * \file GyotoThinDiskPL.h
 * \brief A subclass of ThinDisk emitting according to a powerlaw
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
    Copyright 2012-2018 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoThinDiskPL_H_ 
#define __GyotoThinDiskPL_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class ThinDiskPL; }
}

//#include <GyotoMetric.h>
#include <GyotoThinDisk.h>
#include <GyotoBlackBodySpectrum.h>

/**
 * \class Gyoto::Astrobj::ThinDiskPL
 * \brief Geometrically thin disk with black-body emission
 * 
 * Temperature varies with a power-law from the inner radius outwards:
 *
 * &T; = &T;<SUB>0</SUB> (r<SUB>cur</SUB>/r<SUB>inner</SUB>)<SUP>&alpha;</SUP>
 *
 */
class Gyoto::Astrobj::ThinDiskPL : public Astrobj::ThinDisk {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::ThinDiskPL>;
 private:
  double slope_; ///< Power law index
  double Tinner_; ///< Reference temperature assumed at inner radius
 protected:
  SmartPointer<Spectrum::BlackBody> spectrumBB_; ///< disk black body

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;
  
  ThinDiskPL(); ///< Standard constructor
  
  ThinDiskPL(const ThinDiskPL& ) ;///< Copy constructor
  virtual ThinDiskPL* clone () const; ///< Cloner
  
  virtual ~ThinDiskPL() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  void Slope(double);
  double Slope()const;
  void Tinner(double);
  double Tinner()const;

 public:
  using ThinDisk::emission;
  virtual double emission(double nu_em, double dsem,
			  state_t const &c_ph,double const c_obj[8]=NULL) const;
};

#endif
