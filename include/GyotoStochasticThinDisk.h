/**
 * \file GyotoStochasticThinDisk.h
 * \brief A subclass of ThinDisk emitting according to some stochastic profile
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
    Copyright 2025 Irene Urso

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

#ifndef __GyotoStochasticThinDisk_H_ 
#define __GyotoStochasticThinDisk_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class StochasticThinDisk; }
}

//#include <GyotoMetric.h>
#include <GyotoThinDisk.h>

/**
 * \class Gyoto::Astrobj::StochasticThinDisk
 * \brief A subclass of ThinDisk emitting according to some stochastic profile
 * 
 */
class Gyoto::Astrobj::StochasticThinDisk : public Astrobj::ThinDisk {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::StochasticThinDisk>;
 private:
  double* model_param_; ///< A vector containing an arbitrary number of parameters necessary to compute the disk image
  // Precomputed modal arrays (single 1D vectors)
  std::vector<double> Cmn_;     // Amplitudes
  std::vector<double> Nmn_;     // Amplitudes normalisation
  std::vector<double> Phimn_;   // Phi phases
  std::vector<double> Psimn_;   // Psi phases
  std::vector<double> Kmn_;   // Bessel zeros
 protected:
  unsigned int equationkind_; ///< tag for EquationKind
  unsigned int motionkind_; ///< tag for MotionKind
  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;
  
  StochasticThinDisk(); ///< Standard constructor
  
  StochasticThinDisk(const StochasticThinDisk& ) ;///< Copy constructor
  virtual StochasticThinDisk* clone () const; ///< Cloner
  
  virtual ~StochasticThinDisk() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  virtual std::string equationKind() const ; ///< Get EquationKind
  virtual void equationKind(std::string const&); ///< Set EquationKind
  virtual std::string motionKind() const ; ///< Get MotionKind
  virtual void motionKind(std::string const&); ///< Set MotionKind

  void model_param(std::vector<double> const &v);
  std::vector<double> model_param() const;
  
  virtual void modalQuantities();

 public:
  virtual double spectrum(double const alpha_r, double const alpha_theta, int const m, int const n) const;
  virtual double solution(double const c_obj[8]=NULL) const;
  virtual double envelope(double nu_em, state_t const &c_ph,double const c_obj[8]=NULL) const;
  
  virtual double emission(double nu_em, state_t const &c_ph,double const c_obj[8]=NULL) const;
  virtual void getVelocity(double const pos[4], double vel[4]);

  virtual void processHitQuantities(Photon* ph,
				    state_t const &coord_ph_hit,
				    double const *coord_obj_hit,
				    double dt,
				    Properties* data) const;

};

#endif
