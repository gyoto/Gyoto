/**
 * \file GyotoNumericalMetricLorene.h
 * \brief Base class for 3+1 numerical metrics computed by LORENE
 * 
 *
 */

/*
 *   Copyright (c) 2014-2018, 2020 Frederic Vincent, Thibaut Paumard
 *
 *
 */

#ifndef __GyotoNumericalMetricLoreneMetric_H_
#define __GyotoNumericalMetricLoreneMetric_H_ 

#include <iostream>
#include <fstream>

namespace Gyoto {
  namespace Metric { class NumericalMetricLorene; }
  class FactoryMessenger;
}

// Forward declarations of Lorene classes
namespace Lorene {
  class Scalar;
  class Vector;
  class Sym_tensor;
  class Valeur;
}

#include <GyotoMetric.h>
#include <GyotoWorldline.h>
#include <GyotoSmartPointer.h>

#ifdef GYOTO_USE_XERCES
#include <GyotoRegister.h>
#endif

/**
 * \class Gyoto::NumericalMetricLorene
 * \brief Class for 3+1 numerical metrics computed by LORENE.
 *        This class can handle (so far) any kind of LORENE metric,
 *        stars, collapsing stars, Kerr, boson star e.g.
 */
class Gyoto::Metric::NumericalMetricLorene
: public Gyoto::Metric::Generic
{
  friend class Gyoto::SmartPointer<Gyoto::Metric::NumericalMetricLorene>;

 private:
  char* filename_; ///< Lorene .d data file(s) path
  bool mapet_; ///< Kind of Lorene mapping: 'false' for Map_af, 'true' for Map_et
  bool axisymCirc_; ///< True if sacetime is axisymmetric and circular
  bool bosonstarcircular_; ///< 1 to implement the circular velocity of a boson star
  int has_surface_; ///< 1 if the metric source has a surface
  int has_acceleration_vector_; ///< 1 if the metric source provides an
                                ///< acceleration vector
  int specify_marginalorbits_; ///< 1 if marginal orbits are specified in file
  double horizon_; ///< Value of horizon (or any innermost limit)
  double initial_time_; ///< Time at which (first) metric is given
  Lorene::Scalar** lapse_tab_;
  Lorene::Vector** shift_tab_;
  Lorene::Sym_tensor** gamcov_tab_;
  Lorene::Sym_tensor** gamcon_tab_;
  Lorene::Sym_tensor** kij_tab_;
  double* times_; ///< Coordinate times at which metrics are given
  int nb_times_; ///< Nb of time slices
  Lorene::Valeur** nssurf_tab_; ///< Metric source (e.g. star) surface (if any)
  Lorene::Vector** vsurf_tab_; ///< 4-velocity at surface (if any)
  Lorene::Vector** accel_tab_; ///< 4-acceleration at surface (if any)
  Lorene::Scalar** lorentz_tab_; ///< Lorentz factor at surface (if any)
  Lorene::Valeur** hor_tab_; ///< Apparent horizon (if any)
  double risco_; ///< ISCO coordinate radius
  double rico_; ///< Innermost circular orbit coordinate radius
  double rmb_; ///< Marginally bound orbit coordinate radius

  void free(); ///< deallocate memory

 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;
  NumericalMetricLorene(); ///< Constructor
  NumericalMetricLorene(const NumericalMetricLorene&); ///< Copy constructor
  virtual NumericalMetricLorene* clone() const ;
  virtual ~NumericalMetricLorene() ;        ///< Destructor

  /**
   * Access functions to get or set private attributes
   */
  virtual void setMetricSource();

  void        directory(std::string const &dir) ;
  std::string directory() const ;
  double initialTime() const ;
  void   initialTime(double t0);
  double horizon() const ;
  void   horizon(double t0);
  double rico() const ;
  void   rico(double r0);
  bool hasSurface() const;
  void hasSurface(bool s);
  bool hasAccelerationVector() const;
  void hasAccelerationVector(bool aa);
  bool bosonstarcircular() const;
  void bosonstarcircular(bool);
  bool specifyMarginalOrbits() const;
  void specifyMarginalOrbits(bool s);
  bool mapEt() const;
  void mapEt(bool s);
  bool axisymCirc() const;
  void axisymCirc(bool s);

  Lorene::Vector** getShift_tab() const;
  Lorene::Scalar** getLapse_tab() const;
  Lorene::Sym_tensor** getGamcon_tab() const;
  Lorene::Sym_tensor** getGamcov_tab() const;
  double* getTimes() const;
  int getNbtimes() const;
  Lorene::Valeur** getNssurf_tab() const;
  Lorene::Vector** getVsurf_tab() const;
  Lorene::Vector** getAccel_tab() const;
  Lorene::Scalar** getLorentz_tab() const;
  Lorene::Valeur** getHor_tab() const;
  double getRms() const;
  double getRmb() const;
  void setLapse_tab(Lorene::Scalar* lapse, int ii);
  void setShift_tab(Lorene::Vector* shift, int ii);
  void setGamcov_tab(Lorene::Sym_tensor* gamcov, int ii);
  void setGamcon_tab(Lorene::Sym_tensor* gamcon, int ii);
  void setKij_tab(Lorene::Sym_tensor* kij, int ii);
  void setTimes(double time,int ii);

  virtual double getSpecificAngularMomentum(double rr) const;
  virtual double getPotential(double const pos[4], double l_cst) const;
  
  /**
   * Compute lapse and shift at given coordinates
   */
  virtual void computeNBeta(const double coord[4],double &NN,double beta[3]) const;//Compute lapse and shift at coord

  /**
   * 4-Metric
   */
  using Generic::gmunu;
  double gmunu(const double x[4], int mu, int nu) const ;

  double gmunu(const double x[3], int indice_time, int mu, int nu) const ;

  virtual void gmunu_up(double ARGOUT_ARRAY2[4][4], const double IN_ARRAY1[4]) const ;

  void gmunu_up(double gup[4][4], const double x[4], int indice_time) const ;

  void gmunu_di(const double pos[4],
		double gmunudr[4][4],
		double gmunudth[4][4]) const ;
  
  void gmunu_di(const double pos[4],
		int indice_time,
		double gmunudr[4][4],
		double gmunudth[4][4]) const ;

  virtual void jacobian(double ARGOUT_ARRAY3[4][4][4], const double IN_ARRAY1[4]) const ;

  /**
   * \brief r derivative of contravariant 4-metric
   */
  double gmunu_up_dr(const double x[4], int mu, int nu) const ;

  double gmunu_up_dr(const double x[3], int indice_time, int mu, int nu) const ;
  
  double christoffel(const double coord[4], const int alpha, const int mu,
		     const int nu) const ;
  double christoffel(const double coord[4],
		     const int alpha, 
		     const int mu, const int nu,
		     const int indice_time) const;
  virtual int christoffel(double dst[4][4][4],
			  const double coord[4]) const;
  int christoffel(double dst[4][4][4],
		  const double coord[4],
		  const int indice_time) const;


  /**
   * \brief 3rd order interpolation routine
   *
   * Interpolation at order 3 at point tt, the considered function
   * taking the values "values" at time indices "indices".
   *
   */
  double Interpol3rdOrder(double tt, int indice_time, double values[4]) const;

  /**
   * \brief Computation of horizon value
   */
  double computeHorizon(const double* pos) const;
  double computeHorizon(const double* pos, int indice) const;

  /**
   * F function such as d(coord)/d(tau)=F(coord)
   */
  //using Generic::diff;
  virtual int diff(state_t const &coord, state_t &res, double mass) const;
  int diff(double tt, const double y[7], double res[7]) const ;
  virtual int diff(const double y[7], double res[7], int indice_time) const ;
  virtual int diff31(state_t const &x, state_t &dxdt, double mass) const ;

  /**
   * \brief Yield circular velocity at a given position
   * 
   * Give the velocity of a massive particle in circular orbit at the
   * given position projected onto the equatorial plane. Such a
   * velocity may not exist everywhere (or anywhere) for a given
   * metric. This method is intended to be used by Astrobj classes
   * such as Torus or ThinDisk.
   *
   * This circular velocity should be implemented for all specific
   * numerical metric used.
   *
   * If bosonstarcircular_ is set to true, this method returns the
   * boson star circular velocity.
   *
   * \param coor input: position,
   * \param vel output: velocity,
   * \param dir 1 for corotating, -1 for counterrotating.
   */
  void circularVelocity(double const coor[4], double vel[3],
			double dir) const ;
  void circularVelocity(double const coor[4], double vel[3],
			double dir, int indice_time) const ;

};

#endif
