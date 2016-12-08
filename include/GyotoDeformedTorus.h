/**
 * \file GyotoDeformedTorus.h
 * \brief Slender torus subject to simple time-periodic deformations
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
 *   Copyright (c) 2013 Frederic Vincent
 *
 *
 */

#ifndef __DeformedTorus_h
#define __DeformedTorus_h

#include <GyotoAstrobj.h>
#include <GyotoStandardAstrobj.h>
#include <GyotoKerrBL.h>


namespace Gyoto {
  namespace Astrobj {
    class DeformedTorus;
  };
  class FactoryMessenger;
  namespace Spectrum {
    class Generic;
  }
};

class Gyoto::Astrobj::DeformedTorus
: public Gyoto::Astrobj::Standard {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj::DeformedTorus>;

 private:
  SmartPointer<Gyoto::Metric::KerrBL> gg_;
  SmartPointer<Spectrum::Generic> spectrum_;
  double r_center_;
  int mode_;
  double param_beta_;
  double param_beta_st_;
  double param_eta_;
  int perturb_kind_;
 public:
  GYOTO_OBJECT;
  DeformedTorus();
  DeformedTorus(const DeformedTorus &o);
  virtual ~DeformedTorus();
  virtual DeformedTorus * clone() const ;

  // Standard accessors
  GYOTO_OBJECT_ACCESSORS(SmartPointer<Spectrum::Generic>, spectrum);
  GYOTO_OBJECT_ACCESSORS(double, Rcenter);
  GYOTO_OBJECT_ACCESSORS(double, ParamBeta);
  GYOTO_OBJECT_ACCESSORS(double, ParamBetaSt);
  GYOTO_OBJECT_ACCESSORS(double, ParamEta);
  GYOTO_OBJECT_ACCESSORS(long, Mode);
  GYOTO_OBJECT_ACCESSORS(long, PerturbKind);

  using Generic::metric;
  virtual void metric(Gyoto::SmartPointer<Gyoto::Metric::Generic>);
  virtual double operator()(double const coord[4]) ;
  ///< Called by Astrobj::Generic::Impact()
  virtual void getVelocity(double const pos[4], double vel[4]) ;
  /*virtual int Impact(Gyoto::Photon* ph, size_t index,
    Astrobj::Properties *data=NULL);*/
  double emission(double nuem,double,double*,double*) const;
#endif
};
