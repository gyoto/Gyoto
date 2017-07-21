/**
 * \file GyotoJet.h
 * \brief Simple jet model from Zdziarski, Stawarz & Sikora (MNRAS,2017)
 *
 * This class implements model III of jets of the above paper.
 * It assumes synchrotron radiation emission by power-law electrons.
 */

/*
    Copyright 2017 Frederic Vincent

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

/**
 * \class Gyoto::Astrobj::Jet
 * \brief Simple jet model from Zdziarski, Stawarz & Sikora (MNRAS,2017)
 *
 * This class implements model III of jets of the above paper.
 * It assumes synchrotron radiation emission by power-law electrons.
 */

class Gyoto::Astrobj::Jet
: public Astrobj::Standard,
  public Hook::Listener
{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Jet>;
 private:
  double aa_; ///< Kerr spin
  double baseJetHeight_; ///< Height of the base of the jet (z value)
  double baseJetRadiusOverHeight_; ///< Ratio base cylindrical radius / base height
  double gammaMax_; ///< Max Lorentz factor in jet
  double mdotJet_; ///< Mass-flow rate
  double alfvenRadiusCoef_; ///< Ratio of Alfven radius to light-cylinder radius
  double expoPL_; ///< Power-law index, ne(gamma) \propto gamma^{-p}

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
  void baseJetHeight(double hh);
  double baseJetHeight()const;
  void baseJetRadiusOverHeight(double par);
  double baseJetRadiusOverHeight()const;
  void gammaMax(double gam);
  double gammaMax()const;
  void mdotJet(double mdot);
  double mdotJet()const;
  void alfvenRadiusCoef(double coef);
  double alfvenRadiusCoef()const;
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
  virtual void updateSpin() ;
  virtual void tell(Gyoto::Hook::Teller *msg);
  void JetQuantitiesFromZ(const double zz, double qty[3]) const;
  void JetQuantitiesFromR(const double rr, double qty[2]) const;

  double emissionSynchro_PL_direction(double number_density_PL,
				      double nuem, double nuc,
				      double theta_mag) const;
  double emissionSynchro_PL_averaged(double number_density_PL,
				     double nuem, double nuc) const;
  double absorptionSynchro_PL_direction(double number_density_PL,
					double nuem, double nuc,
					double theta_mag) const ;
  double absorptionSynchro_PL_averaged(double number_density_PL,
				       double nuem, double nuc) const;
};

#endif
