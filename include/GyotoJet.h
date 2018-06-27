/**
 * \file GyotoJet.h
 * \brief Simple jet model from Zdziarski, Stawarz & Sikora (MNRAS,2017)
 *
 * This class implements model III of jets of the above paper.
 * It assumes angle-averaged synchrotron radiation 
 * emission by power-law electrons.
 */

/*
    Copyright 2017-2018 Frederic Vincent

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

#ifndef __GyotoJet_H_ 
#define __GyotoJet_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>

namespace Gyoto{
  namespace Astrobj { class Jet; }
}

//#include <GyotoMetric.h>
#include <GyotoStandardAstrobj.h>
#include <GyotoPowerLawSynchrotronSpectrum.h>


/**
 * \class Gyoto::Astrobj::Jet
 * \brief Simple jet model 
 *
 * This jet assumes angle-averaged synchrotron radiation 
 * emission by power-law electrons.
 */

class Gyoto::Astrobj::Jet
: public Astrobj::Standard,
  public Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Jet>;
 private:
  SmartPointer<Spectrum::PowerLawSynchrotron> spectrumPLSynch_;
  double jetOuterOpeningAngle_; ///< Jet outer opening angle
  double jetInnerOpeningAngle_; ///< Jet inner opening angle
  double jetBaseHeight_; ///< Height of the base of the jet (z value)
  double gammaJet_; ///< Constant Lorentz factor in jet
  double baseNumberDensity_; ///< electron nb density at jet base (cgs)
  double magneticParticlesEquipartitionRatio_; ///< Ratio of magnetic to particles energy density

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;
  
  Jet(); ///< Standard constructor
  
  Jet(const Jet& ) ;///< Copy constructor
  virtual Jet* clone () const; ///< Cloner
  
  virtual ~Jet() ;                        ///< Destructor
  
  // Accessors
  // ---------
 public:
  void jetOuterOpeningAngle(double ang);
  double jetOuterOpeningAngle() const;
  void jetInnerOpeningAngle(double ang);
  double jetInnerOpeningAngle() const;
  void jetBaseHeight(double hh);
  double jetBaseHeight() const;
  void gammaJet(double gam);
  double gammaJet() const;
  void baseNumberDensity(double ne);
  double baseNumberDensity()const;
  void magneticParticlesEquipartitionRatio(double rr);
  double magneticParticlesEquipartitionRatio()const;
  void expoPL(double index);
  double expoPL()const;

 public:
  using Generic::metric;
  virtual void metric(SmartPointer<Metric::Generic>);
    
  virtual double emission(double nu_em, double dsem,
			  double c_ph[8],double c_obj[8]=NULL) const;
  virtual double operator()(double const coord[4]) ;
  virtual void radiativeQ(double Inu[], double Taunu[], 
			  double nu_em[], size_t nbnu,
			  double dsem, double coord_ph[8],
			  double coord_obj[8]=NULL) const ;
  virtual void getVelocity(double const pos[4], double vel[4]) ;
  //virtual void updateSpin() ;
  //virtual void tell(Gyoto::Hook::Teller *msg);

};

#endif
