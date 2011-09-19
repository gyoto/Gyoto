/**
 * \file GyotoMetric.h
 * \brief Base class for metric description
 * 
 * Classes which represent a metric (e.g. Gyoto::Kerr) should inherit
 * from Gyoto::Metric and implement all of the virtual methods. 
 *
 */

/*
    Copyright 2011 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoMetric_H_
#define __GyotoMetric_H_ 

#include <iostream>
#include <fstream>
#include <string>

namespace Gyoto {
  //class Photon;
  class Metric;
}

#include <GyotoWorldline.h>
#include <GyotoSmartPointer.h>
#include <GyotoAstrobj.h>
#include <GyotoRegister.h>

/**
 * \class Gyoto::Metric
 * \brief Base class for metrics
 *
 * Example: class Gyoto::Kerr 
 */
class Gyoto::Metric : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Metric>;

 private:
  std::string kind_;
  double mass_;     ///< Mass yielding geometrical unit (in kg).
  int coordkind_; ///< Kind of coordinates (cartesian-like, spherical-like, unspecified)

 public:
  const std::string getKind() const;
  void setKind(const std::string);
  int getRefCount();
  
  // Constructors - Destructor
  // -------------------------
  //Metric(const Metric& ) ;                ///< Copy constructor
  Metric();
  Metric(const int coordkind);
  Metric(const double mass, const int coordkind);

  virtual ~Metric() ;                        ///< Destructor
  
  // Mutators / assignment
  // ---------------------
  virtual Metric * clone() const ; ///< Virtual copy constructor

  void setMass(const double);        ///< Set mass used in unitLength()
  void setMass(const double, std::string unit);        ///< Set mass used in unitLength()

  // Accessors

  int getCoordKind() const; ///< Get coordinate kind
  void setCoordKind(int coordkind); ///< Set coordinate kind

  double getMass() const;        ///< Get mass used in unitLength()			      
  double unitLength() const ; ///< M * G / c^2, M is in kg, unitLength in meters.

 /**
   * \param coord[4] 4-position (geometrical units);
   * \param v[3]     3-velocity dx1/dx0, dx2/dx0, dx3/dx0;
   * \return tdot = dx0/dtau.
   */
  virtual double SysPrimeToTdot(const double coord[4], const double v[3]) const;
  ///<Compute tdot as a function of dr/dt, dtheta/dt and dphi/dt. Everything is in geometrical units.

  /**
   * \brief Yield circular valocity at a given position
   *
   * \param pos input: position,
   * \param vel output: velocity,
   * \param dir 1 for corotating, -1 for counterrotating.
   */
  virtual void circularVelocity(double const pos[4], double vel[4],
				double dir=1.) const ;

  /**
   * Set coord[4] so that the 4-velocity coord[4:7] is lightlike,
   * i.e. of norm 0. There may be up to two solutions. coord[4] is set
   * to the hightest. The lowest can be retrieved in tdot2. Everything
   * is expressed in geometrical units.
   *
   * \param coord[4] 8-position, coord[4] will be set according to the other elements;
   */
  virtual void nullifyCoord(double coord[8]) const;
  ///< Set tdot (coord[4]) such that coord is light-like. Everything is in geometrical units.

  /**
   * Set coord[4] so that the 4-velocity coord[4:7] is lightlike,
   * i.e. of norm 0. There may be up to two solutions. coord[4] is set
   * to the hightest. The lowest can be retrieved in tdot2. Everything
   * is expressed in geometrical units.
   *
   * \param coord[4] 8-position, coord[4] will be set according to the other elements;
   * \param tdot2    will be set to the smallest solution
   */
  virtual void nullifyCoord(double coord[8], double& tdot2) const;
  ///< Set tdot (coord[4]) such that coord is light-like and return other possible tdot


  /**
   * Compute the scalarproduct of the two quadrivectors u1 and u2 in
   * this Metric, at point pos expressed in coordinate system sys.
   * \param pos[4] 4-position;
   * \param u1[4] 1st quadrivector;
   * \param u2[4] 2nd quadrivector;
   * \return u1*u2
   */
  virtual double ScalarProd(const double pos[4],
		    const double u1[4], const double u2[4]) const; ///< Scalar product

  virtual double Norm3D(double* pos) const; ///< not clear
 
  /**
   * Metrics implementations are free to express lengths and distances
   * in whatever unit they see fit (presumably most often geometrical
   * units). This function returns this unit in SI (meters).
   */

  // Outputs
#ifdef GYOTO_USE_XERCES
  /**
   * Metrics implementations should impement fillElement to save their
   * parameters to XML and call the Metric::fillElement(fmp) for the
   * shared properties
   */

  virtual void fillElement(FactoryMessenger *fmp) ; /// < called from Factory
  void processGenericParameters(Gyoto::FactoryMessenger *fmp) ;
#endif

  /// Display
  //  friend std::ostream& operator<<(std::ostream& , const Metric& ) ;
  //  std::ostream& print(std::ostream&) const ;
 
  /**
   * 
   * \param x[4]     4-position at which to compute the coefficient;
   * \param 0<=mu<=3 1st index of coefficient;
   * \param 0<=nu<=3 2nd index of coefficient;
   * \param sys      coordinate systemp in which x is expressed (see
   *                  GyotoMetric.h)
   * \return Metric coefficient $g_{mu, nu}$ at point x 
   */
  virtual double gmunu(const double * x,
		       int mu, int nu) const
    = 0 ; ///< Metric coefficients

  /**
   * Value of Christoffel symbol $\Gamma^{\alpha}_{\mu\nu}$ at point 
   * $(x_{1},x_{2},x_{3})$
   */  
  virtual double christoffel(const double coord[8],
			     const int alpha, const int mu, const int nu) const = 0;

  virtual int myrk4(Worldline * line, const double coord[8], double h, double res[8]) const;
  
  virtual int myrk4_adaptive(Gyoto::Worldline* line, const double coord[8],
			     double lastnorm, double normref,
			     double coordnew[8], double h0, double& h1) const;

  /**
   * The integrating loop will ask this the Metric through this method
   * whether or not it is happy to conitnue the integration.
   * Typically, the Metric should answer 0 when everything is fine, 1
   * when too close to the event horizon, inside the BH...
   *
   * \param coord[8] coordinates to check.
   */
  virtual int isStopCondition(double const * const coord) const;

  /**
   * F function such as dy/dtau=F(y,cst)
   */
  virtual int diff(const double y[8], double res[8]) const ;

  /**
   * Set Metric-specific constants of motion. Used e.g. in KerrBL.
   */
  virtual void setParticleProperties(Gyoto::Worldline* line,
				     const double * coord) const;
  
  //// REGISTER STUFF ////
 public:
  typedef SmartPointer<Metric> Subcontractor_t(FactoryMessenger*);
  static Register::Entry * Register_;
  static Subcontractor_t* getSubcontractor(std::string);
  static void Register(std::string, Subcontractor_t*);
  static void initRegister();


};

#endif
