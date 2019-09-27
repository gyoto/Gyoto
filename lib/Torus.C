/*
    Copyright 2011, 2018 Thibaut Paumard

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

#include "GyotoTorus.h"
#include "GyotoStandardAstrobj.h"
#include "GyotoUtils.h"
#include "GyotoBlackBodySpectrum.h"
#include "GyotoPowerLawSpectrum.h"
#include "GyotoMetric.h"
#include "GyotoProperty.h"
#include "GyotoFactoryMessenger.h"

#include <float.h>
#include <cmath>
#include <iostream>
#include <cstdlib>

using namespace Gyoto;
using namespace Gyoto::Astrobj;
using namespace std;

GYOTO_PROPERTY_START(Torus,
		     "Geometrical Torus in circular rotation.")
GYOTO_PROPERTY_SPECTRUM(Torus, Spectrum, spectrum,
			"Emission law.")
GYOTO_PROPERTY_SPECTRUM(Torus, Opacity, opacity,
			"Absorption law.")
GYOTO_PROPERTY_DOUBLE(Torus, SmallRadius, smallRadius,
		      "Minor radius, radius of a meridian circle.")
GYOTO_PROPERTY_DOUBLE(Torus, LargeRadius, largeRadius,
		      "Major radius, distance from centre of tube to centre of torus. ")
GYOTO_PROPERTY_END(Torus, Standard::properties)

Torus::Torus() : Standard("Torus"),
	  c_(3.5)
{
  critical_value_ = 0.25; // 0.5*0.5
  safety_value_ = 0.3;
  spectrum_ = new Spectrum::BlackBody(1000000.);
  opacity_ = new Spectrum::PowerLaw(0., 1.);
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();

}

Torus::Torus(const Torus& o)
  : Standard(o),
    c_(o.c_),
    spectrum_(o.spectrum_()?o.spectrum_->clone():NULL),
    opacity_(o.opacity_()?o.opacity_->clone():NULL),
    spectrumThermalSynch_(o.spectrumThermalSynch_()?o.spectrumThermalSynch_->clone():NULL)
{}

Torus::~Torus() {}

Torus* Torus::clone() const { return new Torus(*this); }

double Torus::largeRadius() const { return c_; }
double Torus::largeRadius(string unit) const {
  return Units::FromGeometrical(largeRadius(), unit, gg_);
}
double Torus::smallRadius() const { return sqrt(critical_value_); }
double Torus::smallRadius(string unit) const {
  return Units::FromGeometrical(smallRadius(), unit, gg_);
}

void Torus::largeRadius(double c) { c_ = c; }
void Torus::largeRadius(double c, string unit) {
  largeRadius(Units::ToGeometrical(c, unit, gg_));
}

void Torus::smallRadius(double a) {
  critical_value_ = a*a;
  safety_value_ = critical_value_ * 1.1;
}
void Torus::smallRadius(double c, string unit) {
  smallRadius(Units::ToGeometrical(c, unit, gg_));
}

SmartPointer<Spectrum::Generic> Torus::spectrum() const { return spectrum_; }
void Torus::spectrum(SmartPointer<Spectrum::Generic> sp) {spectrum_=sp;}

SmartPointer<Spectrum::Generic> Torus::opacity() const { return opacity_; }
void Torus::opacity(SmartPointer<Spectrum::Generic> sp) {opacity_=sp;}

double Torus::rMax() {
  if (rmax_==DBL_MAX) {
    rmax_ = 3.*(c_+sqrt(critical_value_));
  }
  return rmax_ ;
}

double Torus::emission(double nu_em, double dsem, state_t const &, double const *) const {
  if (flag_radtransf_) return (*spectrum_)(nu_em, (*opacity_)(nu_em), dsem);
  return (*spectrum_)(nu_em);
}

double Torus::transmission(double nuem, double dsem, state_t const &, double const *) const {
  if (!flag_radtransf_) return 0.;
  double opac = (*opacity_)(nuem);
  if (debug())
    cerr << "DEBUG: Torus::transmission(nuem="<<nuem<<", dsem="<<dsem<<"), "
	 << "opacity=" << opac << "\n";
  if (!opac) return 1.;
  return exp(-opac*dsem);
}

double Torus::integrateEmission(double nu1, double nu2, double dsem,
			       state_t const &, double const *) const {
  if (flag_radtransf_)
    return spectrum_->integrate(nu1, nu2, opacity_(), dsem);
  return spectrum_->integrate(nu1, nu2);
}

double Torus::operator()(double const pos[4]) {
  double drproj, h;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    drproj = pos[1]*sin(pos[2])-c_;
    h = pos[1]*cos(pos[2]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    h = pos[3];
    drproj = sqrt(pos[1]*pos[1]+pos[2]*pos[2])-c_;
    break;
  default:
    GYOTO_ERROR("Torus::distance(): unknown coordinate system kind");
    h=0.,drproj=0.;
  }
  return drproj*drproj + h*h;
}

double Torus::deltaMax(double * coord) {
  double d2 = (*this)(coord);
  if (d2<critical_value_) d2 = critical_value_;
  return 0.1 * sqrt(d2);
}

void Torus::getVelocity(double const pos[4], double vel[4]) {
  double pos2[4] = {pos[0]};
  switch (gg_ -> coordKind()) {
  case GYOTO_COORDKIND_CARTESIAN:
    pos2[1] = pos[1];
    pos2[2] = pos[2];
    pos2[3] = 0.;
    break;
  case GYOTO_COORDKIND_SPHERICAL:
    pos2[1] = pos[1] * sin(pos[2]);
    pos2[2] = M_PI*0.5;
    pos2[3] = pos[3];
    break;
  default:
    GYOTO_ERROR("Torus::getVelocity(): unknown coordkind");
  }
  gg_ -> circularVelocity(pos2, vel);
}

void Torus::radiativeQ(double Inu[], // output
		       double Taunu[], // output
		       double const nu_ems[], size_t nbnu, // input
		       double dsem,
		       state_t const &coord_ph,
		       double const coord_obj[8]) const {
  
  double number_density = 4.8e5;

  double temperature = 9.5e10;

  double thetae = GYOTO_BOLTZMANN_CGS*temperature
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);

  double magnetizationParameter = 0.01;

  double BB = sqrt(4.*M_PI*magnetizationParameter
		   *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
		   *number_density);

  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  //cout << "jet stuff= " << coord_ph[1] << " " << coord_ph[2] << " " << zz << " " << rcyljetbase << " " << rcyl << " " << number_density << " " << thetae << " " << temperatureSlope_ << " " << nu0 << endl;
  //cout << "jet zz,rcyl,th,ph,ne,Te= " <<  zz << " " << rcyl << " " << coord_ph[2] << " " << coord_ph[3] << " " << number_density << " " << temperature << endl;
  // Use that line for Compton study:
  //cout <<  zz << " " << rcyl << " " << number_density << " " << temperature << endl;

  // Emission and absorption synchrotron coefs
  double jnu_synch[nbnu], anu_synch[nbnu];
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    jnu_synch[ii]=-1.;
    anu_synch[ii]=-1.;
  }

  // THERMAL SYNCHROTRON
  spectrumThermalSynch_->temperature(temperature);
  spectrumThermalSynch_->numberdensityCGS(number_density);
  spectrumThermalSynch_->angle_averaged(1); // impose angle-averaging
  spectrumThermalSynch_->angle_B_pem(0.);   // so we don't care about angle
  spectrumThermalSynch_->cyclotron_freq(nu0);
  double besselK2 = bessk(2, 1./thetae);
  spectrumThermalSynch_->besselK2(besselK2);
  
  spectrumThermalSynch_->radiativeQ(jnu_synch,anu_synch,
			nu_ems,nbnu);

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){

    double jnu_tot = jnu_synch[ii],
      anu_tot = anu_synch[ii];

    //cout << "in jet stuff: " << number_density << " " << nu0 << " " << thetae << " " << hypergeom << " " << jnu_tot << " " << anu_tot << " " << dsem << endl;

    //cout << "at r,th= " << coord_ph[1] << " " << coord_ph[2] << endl;
    //cout << "jet jnu anu kappa= " << jnu_tot << " " << anu_tot << endl; //x" " << jnu_tot/anu_tot << " " << dsem << endl;

    // expm1 is a precise implementation of exp(x)-1
    double em1=std::expm1(-anu_tot * dsem * gg_->unitLength());
    Taunu[ii] = em1+1.;
    Inu[ii] = anu_tot == 0. ? jnu_tot * dsem * gg_->unitLength() :
      -jnu_tot / anu_tot * em1;

    if (Inu[ii]<0.)
      GYOTO_ERROR("In Torus::radiativeQ: Inu<0");
    if (Inu[ii]!=Inu[ii] or Taunu[ii]!=Taunu[ii])
      GYOTO_ERROR("In Torus::radiativeQ: Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Taunu[ii]==Taunu[ii]+1.)
      GYOTO_ERROR("In Torus::radiativeQ: Inu or Taunu is infinite");

  }
}
