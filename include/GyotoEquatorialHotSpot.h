/**
 * \file GyotoEquatorialHotSpot.h
 * \brief Equatorial hot spot
 *
 */

/*
    Copyright 2013, 2018 Frederic Vincent & Thibaut Paumard

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

#ifndef __GyotoEquatorialHotSpot_h
#define __GyotoEquatorialHotSpot_h

#include <GyotoAstrobj.h>
#include <GyotoThinDisk.h>
#include <GyotoWorldline.h>
#include <GyotoThermalSynchrotronSpectrum.h>
#include <string>


namespace Gyoto {
  namespace Astrobj {
    class EquatorialHotSpot;
  };
};

class Gyoto::Astrobj::EquatorialHotSpot
: public Gyoto::Astrobj::ThinDisk,
  public Gyoto::Worldline {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::EquatorialHotSpot>;

 private:
  double sizespot_;
  enum beaming_t {IsotropicBeaming=0, NormalBeaming=1,
		  RadialBeaming=2, IsotropicConstant=3};
  beaming_t beaming_;
  double beamangle_;
  SmartPointer<Spectrum::ThermalSynchrotron> spectrumThermalSynch_; // Thermal distribution synchrotron spectrum
  std::string magneticConfig_; ///< Specify the magnetic field configuration for polarisation

 public:
  GYOTO_OBJECT;
  GYOTO_WORLDLINE;
  using Gyoto::Worldline::deltaMax;
  using Gyoto::Astrobj::Generic::deltaMax;
  EquatorialHotSpot();
  EquatorialHotSpot(const EquatorialHotSpot &o);
  virtual ~EquatorialHotSpot();
  virtual EquatorialHotSpot * clone() const ;
  
  // Accessors for the Property list
  void spotRadSize(double t);
  double spotRadSize() const;

  void beaming(std::string const &b);
  std::string beaming() const;

  void beamAngle(double t);
  double beamAngle() const;

  void magneticConfiguration(std::string config);
  std::string magneticConfiguration() const;

  //

  double getMass() const;
  using Generic::metric;
  void metric(SmartPointer<Metric::Generic> gg);
  void setInitialCondition(double coord[8]);

  void getVelocity(double const pos[4], double vel[4]);

  double emission(double nu_em, double dsem,
		  state_t const &,
		  double const coord_obj[8]) const;

  void radiativeQ(double *Inu, double *Qnu, double *Unu,
		  double *Vnu,
		  Eigen::Matrix4d *Onu,
		  double const *nuem , size_t nbnu,
		  double dsem,
		  state_t const &cph,
		  double const *co) const;

  void radiativeQ(double Inu[], // output
		  double Taunu[], // output
		  double const nu_ems[], size_t nbnu, // input
		  double dsem,
		  state_t const &coord_ph,
		  double const coord_obj[8]) const;
  
  // needed for legacy XML files
  virtual int setParameter(std::string name,
		       std::string content,
		       std::string unit);
#ifdef GYOTO_USE_XERCES
  // needed for wait_pos_
  void setParameters(FactoryMessenger* fmp);
  virtual void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const ;
#endif

#endif
};
