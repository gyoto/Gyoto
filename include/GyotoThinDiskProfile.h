/**
 * \file GyotoThinDiskProfile.h
 * \brief A subclass of ThinDisk emitting according to some specified profile
 * that should be hardcoded in emission()
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
    Copyright 2020 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoThinDiskProfile_H_ 
#define __GyotoThinDiskProfile_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class ThinDiskProfile; }
}

//#include <GyotoMetric.h>
#include <GyotoThinDisk.h>

/**
 * \class Gyoto::Astrobj::ThinDiskProfile
 * \brief A subclass of ThinDisk emitting according to some specified profile
 * that should be hardcoded in emission()
 * 
 */
class Gyoto::Astrobj::ThinDiskProfile : public Astrobj::ThinDisk {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::ThinDiskProfile>;
 private:
  double* model_param_; ///< A vector containing an arbitrary number of parameters necessary to compute the disk image
  bool circular_motion_; ///< True if motion is circular, else radial fall
 protected:

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;
  
  ThinDiskProfile(); ///< Standard constructor
  
  ThinDiskProfile(const ThinDiskProfile& ) ;///< Copy constructor
  virtual ThinDiskProfile* clone () const; ///< Cloner
  
  virtual ~ThinDiskProfile() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  bool circularMotion() const;
  void circularMotion(bool circ);

  void model_param(std::vector<double> const &v);
  std::vector<double> model_param() const;

 public:
  virtual double emission(double nu_em, double dsem,
			  state_t const &c_ph,double const c_obj[8]=NULL) const;

  virtual void getVelocity(double const pos[4], double vel[4]);

  virtual void processHitQuantities(Photon* ph,
				    state_t const &coord_ph_hit,
				    double const *coord_obj_hit,
				    double dt,
				    Properties* data) const;

};

#endif
