/**
 * \file GyotoMetric.h
 * \brief Base class for metric description
 * 
 * Classes which represent a metric (e.g. Gyoto::Kerr) should inherit
 * from Gyoto::Metric::Generic and implement all of the virtual
 * methods plus at least one of the gmunu methods and one of the
 * christoffel methods.
 *
 */

/*
    Copyright 2011-2016, 2018-2020 Frederic Vincent, Thibaut Paumard

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
#include <vector>

#include <GyotoSmartPointer.h>
#include <GyotoObject.h>
#include <GyotoAstrobj.h>
#include <GyotoRegister.h>
#include <GyotoHooks.h>
#include <GyotoDefs.h>

namespace Gyoto {
  namespace Metric {

    class Generic;

    /// A function to build instances of a specific Metric::Generic sub-class
    /**
     * This is a more specific version of the
     * SmartPointee::Subcontractor_t type. A Metric::Subcontrator_t is
     * called by the Gyoto::Factory to build an instance of the kind
     * of metric specified in an XML file (see Register()). The
     * Factory and Subcontractor_t function communicate through a
     * Gyoto::FactoryMessenger.
     */
    typedef SmartPointer<Metric::Generic> Subcontractor_t(FactoryMessenger*, std::vector<std::string> const &);


    /** 
     * \brief Subcontractor template
     *
     * Instead of reimplementing the wheel, your subcontractor can simply be
     * Gyoto::Metric::Subcontractor<MyKind>
     *
     * \tparam T Sub-class of Metric::Generic 
     */
    template<typename T> SmartPointer<Metric::Generic> Subcontractor
      (FactoryMessenger* fmp, std::vector<std::string> const &plugins) {
      SmartPointer<T> gg = new T();
      gg -> plugins(plugins);
#ifdef GYOTO_USE_XERCES
      if (fmp) gg -> setParameters(fmp);
#endif
      return gg;
    }

    /// Query the Metric register
    /**
     * Query the Metric register to get the Metric::Subcontractor_t
     * correspondig to a given kind name. This function is normally
     * called only from the Factory. If plugin is specified, only a
     * subcontractor matching both name and plugin will be returned,
     * loading the plug-in if necessary. If plugin is the empty
     * string, then the first subcontractor matching name will be
     * returned, and the name of the plug-in it belongs to will be
     * returned in plugin upon output.
     *
     * \param[in] name e.g. "KerrBL"
     * \param[inout] plugin e.g. "stdplug".
     * \param[in] errmode int=0. If errmode==0, failure to find a
     *        registered Metric by that name is an error. Else, simply
     *        return NULL pointer in that case.
     * \return pointer to the corresponding subcontractor.
     */
    Gyoto::Metric::Subcontractor_t* getSubcontractor(std::string name,
						     std::vector<std::string> &plugin,
						     int errmode=0);

    /// The Metric register
    /**
     * Use the Metric::initRegister() once in your program to
     * initiliaze it, the Metric::Register() function to fill it, and
     * the Metric::getSubcontractor() function to query it.
     */
    extern Register::Entry * Register_;

    /// Make a Metric kind known to the Factory
    /**
     * Register a new Metric::Generic sub-class so that the
     * Gyoto::Factory knows it.
     *
     * \param kind The kind name which identifies this object type in
     * an XML file, as in &lt;Metric kind="name"&gt;
     *
     * \param scp A pointer to the subcontractor, which will
     * communicate with the Gyoto::Factory to build an instance of
     * the class from its XML description
     */
     void Register(std::string kind, Gyoto::Metric::Subcontractor_t* scp);

     /// Empty the Metric register.
     /**
      *  This must be called once. It is called by
      *  Gyoto::Register::init().
      */
     void initRegister();

  }

  /* Documented elswhere */
  class Worldline;
}

/**
 * \namespace Gyoto::Metric
 * \brief Access to metrics
 * 
 * Objects which describe space-time geometry must inherit from the
 * Gyoto::Metric::Generic class.
 *
 * To be usable, a Metric::Generic sub-class should register a
 * Metric::Subcontractor_t function using the Metric::Register()
 * function. See also \ref writing_plugins_page .
 */
/**
 * \class Gyoto::Metric::Generic
 * \brief Base class for metrics
 *
 * Example: class Gyoto::Metric::KerrBL
 *
 * See Gyoto::Metric for an introduction.
 *
 */
class Gyoto::Metric::Generic
: public Gyoto::SmartPointee,
  public Gyoto::Object,
  public Gyoto::Hook::Teller
{
  friend class Gyoto::SmartPointer<Gyoto::Metric::Generic>;

 private:
  double mass_;     ///< Mass yielding geometrical unit (in kg).
  int coordkind_; ///< Kind of coordinates (cartesian-like, spherical-like, unspecified)

  /**
   * \brief Whether some virtual methods are implemented
   *
   * The default implementations of some methods call
   * one-another. This member is used internally to avoid infinite
   * recursion.
   */
  int __defaultfeatures;

 protected:
  double delta_min_; ///< Minimum integration step for the adaptive integrator
  double delta_max_; ///< Maximum integration step for the adaptive integrator

  /**
   * \brief Numerical tuning parameter
   *
   * Ensure that delta (the numerical integration step) is never
   * larger than a fraction of the distance between the current
   * location and the center of the coordinate system.
   *
   * For investigations close to the event horizon, 0.5 is usually
   * fine. If high accuracy is needed long after deflection (weak
   * lensing), then this must be smaller. A good test is to look at a
   * MinDistance map for a FixedStar: it must be smooth.
   */
  double delta_max_over_r_;

  bool keplerian_; ///< 1 if circularVelocity should return the Newtonian Keplerian velocity, in r^-3/2

 protected:
  /**
   * \brief Set coordkind_
   *
   * coordkind(int coordkind) is protected because, for most Metrics,
   * it should not be changed in runtime.
   */
  void coordKind(int coordkind); ///< Set coordinate kind


 public:
  GYOTO_OBJECT;

  int getRefCount();
  
  // Constructors - Destructor
  // -------------------------
  Generic();
  Generic(const int coordkind, const std::string &name); ///< Constructor setting Generic::coordkind_ and kind_
  Generic(Generic const &o); ///< Copy constructor
  virtual ~Generic() ;                        ///< Destructor
  
  // Mutators / assignment
  // ---------------------
  virtual Generic * clone() const ; ///< Virtual copy constructor

  virtual void mass(const double);        ///< Set mass used in unitLength()
  virtual void mass(const double, const std::string &unit);        ///< Set mass used in unitLength()

  // Accessors

  int coordKind() const; ///< Get coordinate kind

  double mass() const;        ///< Get mass used in unitLength()
  double mass(const std::string &unit) const; ///< Get mass used in unitLength()

  /**
   * Metrics implementations are free to express lengths and distances
   * in whatever unit they see fit (presumably most often geometrical
   * units). This function returns this unit in SI (meters).
   */
  double unitLength() const ; ///< M * G / c^2, M is in kg, unitLength in meters
  double unitLength(const std::string &unit) const ; ///< unitLength expressed in specified unit

  /**
   * Returns the marginally bound radius
   * Should be implemented in derived classes if useful
   * If called on the base class, returns an error
   */
  virtual double getRmb() const;

  /**
   * Returns the marginally stable (ISCO) radius 
   * Should be implemented in derived classes if useful
   * If called on the base class, returns an error
   */
  virtual double getRms() const;

  /**
   * Returns the specific angular momentum l=-u_phi/u_t
   * Should be implemented in derived classes if useful
   * If called on the base class, returns an error
   */
  virtual double getSpecificAngularMomentum(double rr) const;

  /**
   * Returns potential W=-ln(|u_t|) for a cst specific angular momentum l_cst
   * Should be implemented in derived classes if useful
   * If called on the base class, returns an error
   */
  virtual double getPotential(double const pos[4], double l_cst) const;

  /**
   * Get delta_min_
   */
  double deltaMin() const;

  /**
   * Set delta_min_
   */
  void deltaMin(double h1);

  /**
   * Get delta_max_
   */
  double deltaMax() const;

  /**
   * Get delta max at a given position
   *
   * \param pos 4-position
   * \param[optional] delta_max_external external constraint on delta_max
   * \return the smallest value between delta_max_,
   * delta_max_external, and R*delta_max_over_r_ where R is pos[1] in
   * spherical coordinates and max(x1, x2, x3) in Cartesian
   * coordinates.
   */
  virtual double deltaMax(double const pos[8], double delta_max_external) const;

  /**
   * Set delta_max_
   */
  void deltaMax(double h1);

  double deltaMaxOverR() const; ///< Get delta_max_over_r_
  void deltaMaxOverR(double t); ///< Set delta_max_over_r_

  bool keplerian() const; ///< Get keplerian_
  void keplerian(bool); ///< Set keplerian_

  virtual void cartesianVelocity(double const coord[8], double vel[3]);
  ///< Compute xprime, yprime and zprime from 8-coordinates

  /**
   * \param coord 4-position (geometrical units);
   * \param v     3-velocity dx1/dx0, dx2/dx0, dx3/dx0;
   * \return tdot = dx0/dtau.
   */
  virtual double SysPrimeToTdot(const double coord[4], const double v[3]) const;
  ///<Compute tdot as a function of dr/dt, dtheta/dt and dphi/dt. Everything is in geometrical units.

  /**
   * \brief Yield circular velocity at a given position
   * 
   * Give the velocity of a massive particle in circular orbit at the
   * given position projected onto the equatorial plane. Such a
   * velocity may not exist everywhere (or anywhere) for a given
   * metric. This method is intended to be used by Astrobj classes
   * such as Torus or ThinDisk.
   *
   * If keplerian_ is set to true, this method should return the
   * Keplerian velcity instead (derived classes should ensure this,
   * see KerrBL::circularVelocity() for instance).
   *
   * The default implementation throws an error if keplerian_ is set
   * to false.
   *
   * \param pos input: position,
   * \param vel output: velocity,
   * \param dir 1 for corotating, -1 for counterrotating.
   */
  virtual void circularVelocity(double const pos[4], double vel[4],
				double dir=1.) const ;

  /**
   * \brief Yield ZAMO velocity at a given position
   *
   * Give the velocity of a zero angular momentul observer (whatever
   * is closest to "at rest"). The default implementation simply
   * projects (1, 0, 0, 0) othogonally along ephi and normalizes it,
   * thus ensuring that vel is orthogonal to ephi.
   *
   * \param pos input: position,
   * \param vel output: velocity,
   */
  virtual void zamoVelocity(double const pos[4], double vel[4]) const ;

  /**
   * Set coord[4] so that the 4-velocity coord[4:7] is lightlike,
   * i.e. of norm 0. There may be up to two solutions. coord[4] is set
   * to the hightest. The lowest can be retrieved using
   * nullifyCoord(double coord[8], double& tdot2) const. Everything is
   * expressed in geometrical units.
   *
   * \param[in,out] coord 8-position, coord[4] will be set according
   * to the other elements;
   */
  virtual void nullifyCoord(double coord[8]) const;
  ///< Set tdot (coord[4]) such that coord is light-like. Everything is in geometrical units.

  /**
   * Set coord[4] so that the 4-velocity coord[4:7] is lightlike,
   * i.e. of norm 0. There may be up to two solutions. coord[4] is set
   * to the hightest. The lowest can be retrieved in tdot2. Everything
   * is expressed in geometrical units.
   *
   * \param[in,out] coord 8-position, coord[4] will be set according
   * to the other elements;
   * \param[out] tdot2    will be set to the smallest solution
   */
  virtual void nullifyCoord(double coord[8], double& tdot2) const;
  ///< Set tdot (coord[4]) such that coord is light-like and return other possible tdot

  /**
   * \brief Normalize fourvelvel to -1
   *
   * First computes threevel as xiprime=xidot/x0dot for i in {1, 2,
   * 3}, then computes x0dot using SyPrimeToTdot, then computes again
   * xidot as xidot=xiprime*x0dot.
   *
   * \param[in,out] coord 8-position, coord[4-7] will be set according
   * to the other elements;
   */
  virtual void normalizeFourVel(double coord[8]) const;

  /**
   * \brief Normalize fourvelvel to -1
   *
   * First computes threevel as xiprime=xidot/x0dot for i in {1, 2,
   * 3}, then computes x0dot using SyPrimeToTdot, then computes again
   * xidot as xidot=xiprime*x0dot.
   *
   * \param[in] pos 4-position;
   * \param[in,out] fourvel 4-velocity, will be renormalized.
   */
  virtual void normalizeFourVel(double const pos[4], double fourvel[4]) const;


  /**
   * Compute the scalarproduct of the two quadrivectors u1 and u2 in
   * this Metric, at point pos expressed in coordinate system sys.
   * \param pos 4-position;
   * \param u1 1st quadrivector;
   * \param u2 2nd quadrivector;
   * \return u1*u2
   */
  virtual double ScalarProd(const double pos[4],
		    const double u1[4], const double u2[4]) const; ///< Scalar product

  /**
   * Compute the norm of the quadrivector u1 in
   * this Metric, at point pos expressed in coordinate system sys.
   * \param[in] pos 4-position;
   * \param[in] u1 quadrivector;
   * \return ||u1||
   */
  double norm(const double pos[4],
		      const double u1[4]) const; ///< Scalar product

  ///\brief multiply vector by scalar
  void multiplyFourVect(double vect[4], double a) const;

  ///\brief add second vector to first one
  void addFourVect(double u1[4], double const u2[4]) const;

  ///\brief project u1 orthogonally to u2 at pos
  void projectFourVect(double const pos[4], double u1[4], double const u2[4]) const;

  /**
   * \brief Computes dual 1-form
   * Compute the dual 1-form of 4-vector.
   * \param IN_ARRAY1_1 4-position;
   * \param IN_ARRAY1_2 quadrivector;
   * \param[out] ARGOUT_ARRAY1 output 1-form
   */
  void dualOneForm(double const IN_ARRAY1_1[4], double const IN_ARRAY1_2[4], double ARGOUT_ARRAY1[4]) const ;

  /**
   * \brief Computes the orthonormal local tetrad of the observer
   *
   * \param obskind  input: kind of observer (eg: "ZAMO","KeplerianObserver"...)
   * \param pos      input: position,
   * \param fourvel output: observer 4-velocity (norm -1)
   * \param screen1 output: first vector in the screen plane
   * \param screen2 output: second vector in the screen plane
   * \param screen3 output: vector normal to the screen
   */
  virtual void observerTetrad(obskind_t obskind,
			      double const pos[4], double fourvel[4],
			      double screen1[4], double screen2[4],
			      double screen3[4]) const ;

  /**
   * \brief Computes the orthonormal local tetrad of the observer
   *
   * \param[in]  pos     position,
   * \param[in]  fourvel observer 4-velocity (norm -1)
   * \param[out] screen1 first vector in the screen plane
   * \param[out] screen2 second vector in the screen plane
   * \param[out] screen3 vector normal to the screen
   */
  virtual void observerTetrad(double const pos[4], double fourvel[4],
			      double screen1[4], double screen2[4],
			      double screen3[4]) const ;

  /**
   * \brief Apply Gram-Schmidt orthonormalization to a basis
   *
   * On input, u0 to u3 must be four non-zero norm, independent
   * 4-vectors. On output, they will form an orthonormal basis.
   *
   * \param[in]     pos position,
   * \param[in,out] u0 basis vector
   * \param[in,out] u1 basis vector
   * \param[in,out] u2 basis vector
   * \param[in,out] u3 basis vector
   */
  void GramSchmidt(double const pos[4], double u0[4],
		   double u1[4], double u2[4], double u3[4]) const;
  // Outputs

  /**
   * \brief Metric coefficients
   *
   * The default implementation calls Metric:: gmunu(double g[4][4], const double * pos) const
   * 
   * \param x  4-position at which to compute the coefficient;
   * \param mu 1st index of coefficient, 0&le;&mu;&le;3;
   * \param nu 2nd index of coefficient, 0&le;&nu;&le;3;
   * \return Metric coefficient g<SUB>&mu;,&nu;</SUB> at point x 
   */
  virtual double gmunu(double const x[4], int mu, int nu) const;

  /**
   * \brief Metric coefficients
   *
   * The default implementation calls double gmunu(const double * x, int mu, int nu) const.
   *
   * \param[out] ARGOUT_ARRAY2 (g) 4x4 array to store the coeefficients
   * \param[in] IN_ARRAY1  (x) 4-position at which to compute the coefficients;
   * \return Metric coefficient g<SUB>&mu;,&nu;</SUB> at point x 
   */
  // Keep argument names for swig!
  virtual void gmunu(double ARGOUT_ARRAY2[4][4], double const IN_ARRAY1[4]) const;

  /**
   * \brief Metric contravariant coefficients
   *
   * The default implementation calls Metric:: gmunu_up(double g[4][4], const double * pos) const
   *
   * \param x  4-position at which to compute the coefficient;
   * \param mu 1st index of coefficient, 0&le;&mu;&le;3;
   * \param nu 2nd index of coefficient, 0&le;&nu;&le;3;
   * \return Metric coefficient g<SUP>&mu;,&nu;</SUP> at point x
   */
  virtual double gmunu_up(double const x[4], int mu, int nu) const;

  /**
   * \brief Metric contravariant coefficients
   *
   * The default implementation inverts the covariant coefficients matrix.
   **/
  // Keep argument names for swig!
  virtual void gmunu_up(double ARGOUT_ARRAY2[4][4], const double IN_ARRAY1[4]) const;

  /**
   * \brief Derivatives of the metric covariant coefficients
   *
   * The default implementation evaluates them numerically. The gmunu
   * matrix is assumed to be symmetrical but no other assumptions are
   * made at the moment.
   **/
  // Keep argument names for swig!
  virtual void jacobian(double ARGOUT_ARRAY3[4][4][4], const double IN_ARRAY1[4]) const;

  /**
   * \brief gmunu_up() and jacobian() in one go
   **/
  // Keep argument names for swig!
  virtual void gmunu_up_and_jacobian(double ARGOUT_ARRAY2[4][4], double ARGOUT_ARRAY3[4][4][4], const double IN_ARRAY1[4]) const;

  /**
   * \brief Computes lapse scalar and shift vector at coord.
   *
   **/
  virtual void computeNBeta(const double coord[4],double &NN,double beta[3]) const;

  /**
   * \brief Chistoffel symbol
   *
   * Value of Christoffel symbol
   * &Gamma;<SUP>&alpha;</SUP><SUB>&mu;&nu;</SUB> at point
   * (x<SUB>1</SUB>, x<SUB>2</SUB>, x<SUB>3</SUB>).
   */  
  virtual double christoffel(const double coord[4],
			     const int alpha, const int mu, const int nu) const;

  /**
   * \brief Chistoffel symbol
   *
   * Value of Christoffel symbol
   * &Gamma;<SUP>&alpha;</SUP><SUB>&mu;&nu;</SUB> at point
   * (x<SUB>1</SUB>, x<SUB>2</SUB>, x<SUB>3</SUB>).
   *
   * \return 1 on error, 0 otherwise
   */  
  virtual int christoffel(double dst[4][4][4], const double coord[4]) const ;



  /**
   * \brief RK4 integrator
   */
  virtual int myrk4(Worldline * line, state_t const &coord, double h, state_t &res) const;
  /// Obsolete, update your code
  virtual int myrk4(Worldline * line, const double coord[8], double h, double res[8]) const = delete;
  
  /**
   * \brief RK4 integrator with adaptive step
   */
  virtual int myrk4_adaptive(Gyoto::Worldline* line, state_t const &coord,
			     double lastnorm, double normref,
			     state_t &coordnew, double h0, double& h1,
			     double deltamax=GYOTO_DEFAULT_DELTA_MAX) const;
  /// Obsolete, update your code
  virtual int myrk4_adaptive(Gyoto::Worldline* line, const double coord[8],
                             double lastnorm, double normref,
			     double coordnew[8], double h0, double& h1,
                             double deltamax=GYOTO_DEFAULT_DELTA_MAX) const = delete;


  /**
   * \brief Check whether integration should stop
   *
   * The integrating loop will ask this the Metric through this method
   * whether or not it is happy to continue the integration.
   * Typically, the Metric should answer 0 when everything is fine, 1
   * when too close to the event horizon, inside the BH...
   *
   * \param coord 8-coordinate vector to check.
   */
  virtual int isStopCondition(double const coord[8]) const;

  /**
   * \brief F function such as dx/dt=F(x,cst)
   */
  virtual int diff(state_t const &x, state_t &dxdt, double mass) const ;
  /// Obsolete, update your code
  virtual int diff(state_t const &x, state_t &dxdt)  const = delete;
  virtual int diff(const double y[8], double res[8]) const = delete ;

  /**
   * \brief F function such as dx/dt=F(x,cst) for 3+1 case
   */
  virtual int diff31(state_t const &x, state_t &dxdt, double mass) const ;

  /**
   * \brief Set Metric-specific constants of motion. Used e.g. in KerrBL.
   */
  virtual void setParticleProperties(Gyoto::Worldline* line,
				     double const coord[8]) const;
  

};

#endif
