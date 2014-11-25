/**
 * \file GyotoNumericalMetricLorene.h
 * \brief Base class for 3+1 numerical metrics computed by LORENE
 * 
 *
 */

/*
 *   Copyright (c) 2012 Frederic Vincent
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
#include <GyotoWIP.h>

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
: public WIP, public Gyoto::Metric::Generic
{
  friend class Gyoto::SmartPointer<Gyoto::Metric::NumericalMetricLorene>;

 private:
  char* filename_; ///< Lorene .d data file(s) path
  char* mapkind_; ///< Kind of Lorene mapping Map_af or Map_et
  int has_surface_; ///< 1 if the metric source has a surface
  int specify_marginalorbits_; ///< 1 if marginal orbits are specified in file
  double horizon_; ///< Value of horizon (or any innermost limit)
  double r_refine_; ///< Refine integration below this r
  double h0_refine_; ///< Imposed integration step for refined integration
  int refine_; ///< 1 if refined integration needed
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
  Lorene::Scalar** lorentz_tab_; ///< Lorentz factor at surface (if any)
  Lorene::Valeur** hor_tab_; ///< Apparent horizon (if any)
  double risco_; ///< ISCO coordinate radius
  double rmb_; ///< Marginally bound orbit coordinate radius

  void free(); ///< deallocate memory

 public:
  NumericalMetricLorene(); ///< Constructor
  NumericalMetricLorene(const NumericalMetricLorene&); ///< Copy constructor
  virtual NumericalMetricLorene* clone() const ;
  virtual ~NumericalMetricLorene() ;        ///< Destructor

  /**
   * Access functions to get or set private attributes
   */
  virtual void setMetricSource();
  char const * getFileName() const;
  Lorene::Vector** getShift_tab() const;
  Lorene::Scalar** getLapse_tab() const;
  Lorene::Sym_tensor** getGamcon_tab() const;
  double* getTimes() const;
  int getNbtimes() const;
  Lorene::Valeur** getNssurf_tab() const;
  Lorene::Vector** getVsurf_tab() const;
  Lorene::Scalar** getLorentz_tab() const;
  Lorene::Valeur** getHor_tab() const;
  double getRisco() const;
  double getRmb() const;
  void setLapse_tab(Lorene::Scalar* lapse, int ii);
  void setShift_tab(Lorene::Vector* shift, int ii);
  void setGamcov_tab(Lorene::Sym_tensor* gamcov, int ii);
  void setGamcon_tab(Lorene::Sym_tensor* gamcon, int ii);
  void setKij_tab(Lorene::Sym_tensor* kij, int ii);
  void setTimes(double time,int ii);
  
  /**
   * Runge-Kutta integrator at order 4
   */
  //using Generic::myrk4; //--> why using this?
  virtual int myrk4(double tt, const double coord[7], double h, double res[7]) const;
  virtual int myrk4(Worldline* line, const double coord[8], 
	    double h, double res[8]) const;

  /**
   * Adaptive Runge-Kutta
   */
  int myrk4_adaptive(Gyoto::Worldline* line, const double coord[8], double lastnorm, double normref, double coordnew[8], double h0, double& h1, double h1max) const;

  int myrk4_adaptive(double tt, const double coor[7], double lastnorm, double normref, double coornew[7], const double cst[2], double& tdot_used, double h0, double& h1, double& hused, double h1max) const;
  ///< With energy integration also, coor=[E,r,th,ph,dE/dt,Vr,Vth,Vph]

  /**
   * Reverse spatial vector if going throough 0, without horizon
   */
  void reverseR(double tt, double coord[4]) const;

  /**
   * Compute lapse and shift at given coordinates
   */
  void computeNBeta(const double coord[4],double &NN,double beta[3]) const;//Compute lapse and shift at coord

  /**
   * 4-Metric
   */
  double gmunu(const double x[4], int mu, int nu) const ;

  double gmunu(const double x[3], int indice_time, int mu, int nu) const ;

  /**
   * r derivative of contravariant 4-metric
   */
  double gmunu_up_dr(const double x[4], int mu, int nu) const ;

  double gmunu_up_dr(const double x[3], int indice_time, int mu, int nu) const ;
  
  double christoffel(const double coord[8], const int alpha, const int mu, 
		     const int nu) const ;
  /**
   * 3-Christoffels
   */
  double christoffel3(const double coord[6], const int indice_time, const int ii, 
		      const int jj, const int kk) const ; //3D Christoffel

  void setParticleProperties(Worldline * line, const double* coord) const;

  /**
   * 3rd order interpolation routine
   */
  double Interpol3rdOrder(double tt, int indice_time, double values[4]) const;
  /*Interpolation at order 3 at point tt, the considered function 
    taking the values "values" at time indices "indices".*/

  /**
   * Computation of horizon value
   */
  double computeHorizon(const double* pos) const;
  double computeHorizon(const double* pos, int indice) const;

  /**
   * F function such as d(coord)/d(tau)=F(coord)
   */
  using Generic::diff;
  int diff(double tt, const double y[7], double res[7]) const ;
  virtual int diff(const double y[7], double res[7], int indice_time) const ;

  virtual int setParameter(std::string, std::string, std::string);
#ifdef GYOTO_USE_XERCES
  virtual void fillElement(FactoryMessenger *fmp); /// < called from Factory
  virtual void setParameters(FactoryMessenger *fmp);
#endif

};

#endif
