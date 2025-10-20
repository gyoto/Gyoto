/**
 *  \file GyotoKonoplyaRezzollaZhidenko.h
 *  \brief Axisymmetric parametrised metric of Konoplya\&Rezzolla\&Zhidenko 2016
 *         See the paper: PRD, 93, 064015
 *         Kerr metric retrieved for the functions given in Table 1 of Cárdenas-Avendaño & Held PRD109(2024)064052
 *         when all the horizon and asymptotics deformation parameters are set to zero
 *         For the time being we allow non-zero values for :
 *         Horizon deformations parameters of Ni et al. JCAP09(2016)014
 *         Asymptotics physical deformations parameters of eq. B26-B27 in Cárdenas-Avendaño & Held PRD109(2024)064052
 *         Only one deformation parameter at a time
 */

/*
    Copyright 2025 Irene Urso

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

#ifndef __GyotoKonoplyaRezzollaZhidenko_h
#define __GyotoKonoplyaRezzollaZhidenko_h

#include <GyotoMetric.h>

namespace Gyoto {
  namespace Metric {
    class KonoplyaRezzollaZhidenko;
  };
};

class Gyoto::Metric::KonoplyaRezzollaZhidenko
: public Gyoto::Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::KonoplyaRezzollaZhidenko>;
 private:
  double spin_;
  double spin2_;
  double spin3_;
  double spin4_;
  double rms_; ///< Provide marginally stable orbits if needed
  double* deltashorizon_; ///< The ẟ-parameter vector [δ1,ẟ2,ẟ3,ẟ4,ẟ5,ẟ6]=[a01,w01,w21,b01,b21,->a21&k21] used in Ni et al. JCAP09(2016)014
  double* deltasasymptotics_; ///< The parameter vector [δepsilon0,ẟw00,ẟa00,ẟb00] used in Cárdenas-Avendaño & Held PRD109(2024)064052
 public:
  GYOTO_OBJECT;
  KonoplyaRezzollaZhidenko();
  KonoplyaRezzollaZhidenko(const KonoplyaRezzollaZhidenko & orig);
  virtual ~KonoplyaRezzollaZhidenko();
  virtual KonoplyaRezzollaZhidenko * clone() const ;

  // accessors
  void spin(const double val); ///< Set spin
  double spin() const ; ///< Returns spin
  GYOTO_OBJECT_ACCESSORS(double, rms);
  void deltashorizon(std::vector<double> const &v);
  std::vector<double> deltashorizon() const;
  void deltasasymptotics(std::vector<double> const &v);
  std::vector<double> deltasasymptotics() const;


  using Generic::gmunu;
  double gmunu(double const x[4], int mu, int nu) const ;
  void gmunu_up(double ARGOUT_ARRAY2[4][4], const double IN_ARRAY1[4]) const ;
  double gmunu_up(double const x[4], int mu, int nu) const ;
  double Definer0() const;
  enum class AsymptoticParameter {epsilon0, k00, w00, a20, a00, b00, epsilon2, b20, w20, k20};
  double DefineAsymptoticParameters(AsymptoticParameter, const double r0) const;
  enum class HorizonParameter {a01, a21, k21, k22, k23, w01, w21, b01, b21};
  double DefineHorizonParameters(HorizonParameter, const double r0) const;
  double N2(const double rr, const double th) const;
  double B(const double rr,  const double th) const;
  double Sigma(const double rr,  const double th) const;
  double W(const double rr,  const double th) const;
  double K2(const double rr,  const double th) const;
  double drN2(const double rr, const double th) const;
  double dthN2(const double rr, const double th) const;
  double drB(const double rr, const double th) const;
  double dthB(const double rr, const double th) const;
  double drSigma(const double rr, const double th) const;
  double dthSigma(const double rr, const double th) const;
  double drW(const double rr, const double th) const;
  double dthW(const double rr, const double th) const;
  double drK2(const double rr, const double th) const;
  double dthK2(const double rr, const double th) const;
  double KeplerianSpecificAngularMomentum(const double rr, const double th) const;
  double KeplerianAngularVelocity(const double rr, const double th, const double ell) const;
  double KeplerianEnergy(const double rr, const double th, const double Omega) const;
  double KeplerianAngularMomentum(const double rr, const double th, const double Omega) const;
  using Generic::christoffel;
  int christoffel(double dst[4][4][4], double const pos[4]) const ;
  int isStopCondition(double const coord[8]) const;
  virtual double getRms() const;
  virtual void circularVelocity(double const pos[4], double vel [4],
				double dir=1.) const ;

 
#endif
};
