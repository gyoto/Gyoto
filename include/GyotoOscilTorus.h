/**
 * \file GyotoOscilTorus.h
 * \brief Slender torus subject to realistic Blaes 2006 oscillation modes.
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
 *   Copyright (c) 2013, 2018 Frederic Vincent & Thibaut Paumard
 *
 *
 */

#ifndef __OscilTorus_h
#define __OscilTorus_h

#include <GyotoAstrobj.h>
#include <GyotoStandardAstrobj.h>
#include <GyotoKerrBL.h>


namespace Gyoto {
  namespace Astrobj {
    class OscilTorus;
  };
  class FactoryMessenger;
  namespace Spectrum {
    class Generic;
  }
};

class Gyoto::Astrobj::OscilTorus
: public Gyoto::Astrobj::Standard,
  public Hook::Listener {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::OscilTorus>;

 private:
  // Members corresponding to properties:  
  /**
   * \brief Large Radius
   *
   * Distance from the center of the coordinate system to the center
   * of the torus tube. The (square of the) radius of a vertical
   * cross-section is stored in critical_value_.
   */
  double c_;
  unsigned long mode_;
  double polycst_; ///< Polytropic constant
  double polyindex_; ///< Polytropic index
  double central_density_; ///< Central density
  enum perturb_t {Radial=1, Vertical=2, X=3, Plus=4, Breathing=5};
  perturb_t perturb_kind_;
  std::string emitting_area_; ///< Only for mode=0, file containing time series of cross section area
  double perturb_intens_; ///< Perturbation intensity

  // Cached values:
  SmartPointer<Gyoto::Metric::KerrBL> kerrbl_;
  std::vector<double> tt_;
  std::vector<double> area_; // tt_ and area_ contain area of cross section at time tt
  size_t nbt_; ///< numberof tt_
  int with_cross_; ///< is 1 if cross section data are given

  double sigma_; ///< perturbation rescaled pulsation
  double alpha_; ///< perturbation normalization (oder-unity)
  double w1_; ///< factors appearing in perturbed surf func
  double w2_; 
  double omr2_; ///< epicyclic freq at torus center
  double omth2_;
  double Omegac_; ///< Omega and l at torus center
  double lc_;
  double g_rr_; ///< metric coef at torus center
  double g_thth_;
  bool hold_;
  
 public:
  GYOTO_OBJECT;
  OscilTorus();
  OscilTorus(const OscilTorus &o);
  virtual ~OscilTorus();
  virtual OscilTorus * clone() const ;

  GYOTO_OBJECT_ACCESSORS_UNIT(largeRadius);
  GYOTO_OBJECT_ACCESSORS(unsigned long, mode);
  GYOTO_OBJECT_ACCESSORS(double, polyCst);
  GYOTO_OBJECT_ACCESSORS(double, polyIndex);
  GYOTO_OBJECT_ACCESSORS(double, centralDensity);
  GYOTO_OBJECT_ACCESSORS_STRING(perturbKind);
  GYOTO_OBJECT_ACCESSORS(double, perturbIntens);
  GYOTO_OBJECT_ACCESSORS_STRING(emittingArea);
  using Generic::metric;
  virtual void metric(Gyoto::SmartPointer<Gyoto::Metric::Generic>);

  virtual double operator()(double const coord[4]) ;
  virtual void getVelocity(double const pos[4], double vel[4]) ;
  double emission(double nuem,double,state_t const &,double const *) const;

#ifdef GYOTO_USE_XERCES
  virtual void setParameters(Gyoto::FactoryMessenger *fmp) ;
#endif

  virtual void updateCachedValues();
  void computeXbYb(const double * pos, double & xb, double & yb) ;

  // Hook::Listener API //
 public:
  /**
   * \brief Update cached values
   *
   * Calls updateCachedValues().
   *
   * See Hook::Listener::tell()
   */
  virtual void tell(Gyoto::Hook::Teller *msg);

};

#endif
