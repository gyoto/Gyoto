/**
 * \file GyotoPlasmoid.h
 * \brief Plasmoid sphere formed by magnetic reconnection following a Star orbit, emitting synchrotron,
 * with two distributions of electrons:
 * one thermal at "low" temperature and one kappa at "high" temperature
 *
 */

/*
    Copyright 2019 Frederic Vincent, Thibaut Paumard, Nicolas Aimar

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


#ifndef __GyotoPlasmoid_H_ 
#define __GyotoPlasmoid_H_ 

namespace Gyoto{
  namespace Astrobj { class Plasmoid; }
}

#include <iostream>
#include <fstream>
#include <iomanip>
#include <GyotoMetric.h>
#include <GyotoUniformSphere.h>
#include <GyotoFitsRW.h>
#include <GyotoKappaDistributionSynchrotronSpectrum.h>
//#include <GyotoThermalSynchrotronSpectrum.h>
#ifdef GYOTO_USE_CFITSIO
#include <fitsio.h>
#endif

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

#include <string>

/**
 * \class Gyoto::Astrobj::Plasmoid
 * \brief Plasmoid Shere of plasma emitting synchrotron, following 
 * a trajectory specified in getVelocity (non-geodesic a priori)
 *
 */
class Gyoto::Astrobj::Plasmoid :
  public FitsRW,
  public Gyoto::Astrobj::UniformSphere{
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::Plasmoid>;
  
  // Data : 
  // -----
 private:
  double* posIni_; // 4-position of the plasmoid in spherical coordinates
  double* fourveldt_; // 4-velocity of the plasmoid in spherical coordinates (dxi/dt, not dtau) 
  std::string flag_; // type of motion "helical" or "equatorial"
  bool posSet_;
  double t_inj_;
  double radiusMax_; // Maximun radius of the Plasmoid in geometrical units
  std::string varyRadius_;
  // FITS FILE Quantities
  std::string filename_;
  double* freq_array_;
  double* jnu_array_;
  double* anu_array_;

  // Constructors - Destructor
  // -------------------------
 public:
  GYOTO_OBJECT; // This object has a (non-inherited) Property list

 /**
  * Create Plasmoid object with undefined initial conditions. One needs to
  * set the coordinate system, the metric, the type of motion, and the initial position
  * and velocity before integrating the orbit. initCoord()
  * can be used for that.
  */
  Plasmoid(); ///< Default constructor
  
  Plasmoid(const Plasmoid& orig); ///< Copy constructor
  virtual Plasmoid * clone() const ;

  virtual ~Plasmoid() ;                        ///< Destructor
  
 public:
  virtual std::string className() const ; ///< "Plasmoid"
  virtual std::string className_l() const ; ///< "inflate_star"

 public:
  void motionType(std::string const type);
  SmartPointer<Metric::Generic> metric() const;
  void metric(SmartPointer<Metric::Generic> gg);
  void initPosition(std::vector<double> const &v);
  std::vector<double> initPosition() const;
  void initVelocity(std::vector<double> const &v);
  std::vector<double> initVelocity() const;
  void initCoord(std::vector<double> const &v);
  std::vector<double> initCoord() const;
  void radiusMax(double rr);
  double radiusMax() const;
  void Radius(std::string vary);
  
  virtual void radiativeQ(double Inu[], double Taunu[], 
			  double const nu_em[], size_t nbnu,
			  double dsem, state_t const &coord_ph,
			  double const coord_obj[8]=NULL) const;

  void getCartesian(double const * const dates, size_t const n_dates,
          double * const x, double * const y,
          double * const z, double * const xprime=NULL,
          double * const yprime=NULL,
          double * const zprime=NULL);

  void getVelocity(double const pos[4], double vel[4]);

  int Impact(Gyoto::Photon* ph, size_t index,
         Astrobj::Properties *data=NULL);

  void file(std::string const &f);

  std::vector<size_t> fitsRead(std::string filename);

};


#endif
