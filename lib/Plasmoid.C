/*
    Copyright 2019, 2020 Frederic Vincent, Thibaut Paumard, Nicolas Aimar

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

#include "GyotoUtils.h"
#include "GyotoPlasmoid.h"
#include "GyotoPhoton.h"
#include "GyotoFactoryMessenger.h"

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <float.h>
#include <sstream>
#include <string.h>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Plasmoid, "Synchrotron-emitting orbiting plasmoid heated by magnetic reconnection")
GYOTO_PROPERTY_DOUBLE_UNIT(Plasmoid, NumberDensity, numberDensity,
               "cgs number density, constant through plasmoid")
GYOTO_PROPERTY_DOUBLE_UNIT(Plasmoid, TimeRef, timeRef,
                "time of reconnection event")
GYOTO_PROPERTY_DOUBLE(Plasmoid, TemperatureInitial, temperatureInitial,
              "temperature initial, constant through initial population of e-")
GYOTO_PROPERTY_DOUBLE(Plasmoid, TemperatureReconnection, temperatureReconnection,
              "temperature reconnection, constant through reconnected population of e-")
GYOTO_PROPERTY_DOUBLE(Plasmoid, MagnetizationParameter,
              magnetizationParameter,
              "magnetization parameter")
GYOTO_PROPERTY_DOUBLE(Plasmoid, KappaIndex, kappaIndex,
		      "PL index of kappa-synchrotron")
GYOTO_PROPERTY_END(Plasmoid, Star::properties)

Plasmoid::Plasmoid() :
  Star(),
  numberDensity_cgs_(1.),
  timeRef_(-10.),
  temperatureInitial_(1.),
  temperatureReconnection_(1.),
  magnetizationParameter_(1.),
  kappaIndexMin_(3.5),
  spectrumKappaSynch_(NULL),
  spectrumThermalSynch_(NULL)
{
  kind_="Plasmoid";
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
  spectrumKappaSynch_ = new Spectrum::KappaDistributionSynchrotron();
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
}

Plasmoid::Plasmoid(const Plasmoid& orig) :
  Star(orig),
  numberDensity_cgs_(orig.numberDensity_cgs_),
  timeRef_(orig.timeRef_),
  temperatureInitial_(orig.temperatureInitial_),
  temperatureReconnection_(orig.temperatureReconnection_),
  magnetizationParameter_(orig.magnetizationParameter_),
  kappaIndexMin_(orig.kappaIndexMin_),
  spectrumKappaSynch_(NULL),
  spectrumThermalSynch_(NULL)
{
  if (orig.spectrumKappaSynch_()) spectrumKappaSynch_=orig.spectrumKappaSynch_->clone();
  if (orig.spectrumThermalSynch_()) spectrumThermalSynch_=orig.spectrumThermalSynch_->clone();
}


Plasmoid* Plasmoid::clone() const { return new Plasmoid(*this); }

Plasmoid::~Plasmoid() {
  if (debug()) cerr << "DEBUG: Plasmoid::~Plasmoid()\n";
}

string Plasmoid::className() const { return  string("Plasmoid"); }
string Plasmoid::className_l() const { return  string("Plasmoid"); }


void Plasmoid::radiativeQ(double Inu[], // output
              double Taunu[], // output
              double const nu_ems[], size_t nbnu, // input
              double dsem,
              state_t const &coord_ph,
              double const coord_obj[8]) const {

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  double tcur=coord_ph[0]; //*GMoc3/60.; // in min
  double Rad = radius("cm");
  double t_inj=Rad/(GYOTO_C_CGS*pow(magnetizationParameter_/(magnetizationParameter_+1),0.5))/60.; //in min; //injection time, i.e. time during e- are heated and accelerated due to the reconnection, see [D. Ball et al., 2018]

  double number_density_ini=numberDensity_cgs_; // number density of e- which follow the initial thermal distribution, before the reconnection
  double number_density_rec=0.; // umber density of e- which follow the kappa distribution after the reconnection

  double reconnection_rate = 0.1;
  double n_dot=reconnection_rate*numberDensity_cgs_*pow(magnetizationParameter_/(magnetizationParameter_+1),0.5)/radius(); //"Reconnection rate", see [D. Ball et al., 2020] (ie Psaltis paper)

  double thetae_ini = GYOTO_BOLTZMANN_CGS*temperatureInitial_
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS); //Dimentionless temperature of the initial thermal population of e-
  double thetae_rec = GYOTO_BOLTZMANN_CGS*temperatureReconnection_
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS); //Dimentionless temperature of the population of e- after the reconnection

  
  double BB = sqrt(8.*M_PI*magnetizationParameter_
           *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
           *numberDensity_cgs_); // Magnetic field
  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  double sigma_thomson=8.*M_PI*pow(GYOTO_ELECTRON_CLASSICAL_RADIUS_CGS,2.)/3.; // Thomson's cross section 
  double AA = (4./3.*sigma_thomson*GYOTO_C_CGS*pow(BB,2.)/8.*M_PI)/(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS); // Coefficient of integration from [D. Ball et al., 2020] for cooling
  cout << AA << endl;

  double gamma_max = DBL_MAX; // Initial value of gamma_max=gamma_min=3*thetae_rec
  double gamma_max_0 = DBL_MAX;
  double kappaIndex = 100.; // Initialisation of the power law indice; do not set over 100 -> pow(xx,kappaIndex)=inf
 
  // COMPUTE VALUES IN FUNCTION OF PHASE
  if (tcur<=timeRef_)
  {
    number_density_ini=numberDensity_cgs_;
    number_density_rec=0.;
  }
  else if (tcur<=timeRef_+t_inj) // HEATING TIME
  {
    number_density_rec=n_dot*(tcur-timeRef_);

    number_density_ini=numberDensity_cgs_ - number_density_rec;

    //talk with Andreas and Anton; evolution of the kappa index
    kappaIndex = kappaIndex - (kappaIndex - kappaIndexMin_ )*(tcur-timeRef_)/t_inj; //linear
    //kappaIndex = kappaIndex*exp(-(tcur-timeRef_)) + kappaIndexMin_; // decreasing exponential
    //kappaIndex = kappaIndexMin_; // No varying kappa
  }
  else // COOLING TIME
  {
    // evolution of the number densities
    number_density_rec=n_dot*t_inj;
    number_density_ini=numberDensity_cgs_ - number_density_rec;

    kappaIndex = kappaIndexMin_;
    

    thetae_rec=max(thetae_rec*pow(1+AA*3.*thetae_rec*(tcur-(t_inj+timeRef_)),-1),thetae_ini); //see [Radiative Process in Astrophysics, G. Rybicki and Lightman, 1986];theta=gamma/3

    gamma_max_0 = 1.e10; //pow(DBL_MIN*(4.*M_PI*(1.-pow(DBL_MAX,1.-(kappaIndex-1.))))/((kappaIndex-1.)-1.),-1./(kappaIndex-1.));
    gamma_max = max(gamma_max_0*pow(1+AA*gamma_max_0*(tcur-(t_inj+timeRef_)),-1.),3*thetae_rec);
    
    /*double nu = *nu_ems;
    double gamma0 = sqrt(nu/nu0)*pow(1. - AA*sqrt(nu/nu0)*(tcur-(t_inj+timeRef_)),-1.);
    *nu_ems = pow(gamma0,2.)*nu0;*/

    
  }

  double hypergeom = Gyoto::hypergeom(kappaIndex, thetae_rec);
  double bessel_K2 = Gyoto::bessk(2,thetae_ini);
  
  // Defining jnus, anus
  double jnu_synch_kappa[nbnu], anu_synch_kappa[nbnu];
  double jnu_synch_th[nbnu], anu_synch_th[nbnu];
  
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    // [ exp(-anu*ds) will explose ]
    jnu_synch_kappa[ii]=-1.;
    anu_synch_kappa[ii]=-1.;
    jnu_synch_th[ii]=-1.;
    anu_synch_th[ii]=-1.;
  }

  //KAPPA SYNCHRO
  spectrumKappaSynch_->kappaindex(kappaIndex);
  spectrumKappaSynch_->angle_averaged(1);
  spectrumKappaSynch_->angle_B_pem(0.); // avg so we don't care
  spectrumKappaSynch_->cyclotron_freq(nu0);
  spectrumKappaSynch_->numberdensityCGS(number_density_rec);
  spectrumKappaSynch_->thetae(thetae_rec);
  spectrumKappaSynch_->hypergeometric(hypergeom);
  spectrumKappaSynch_->gamma_max(gamma_max);
  
  spectrumKappaSynch_->radiativeQ(jnu_synch_kappa,anu_synch_kappa,
                  nu_ems,nbnu);

  //THERMAL SYNCHRO
  spectrumThermalSynch_->temperature(thetae_ini);
  spectrumThermalSynch_->angle_B_pem(0.); // avg so we don't care
  spectrumThermalSynch_->cyclotron_freq(nu0);
  spectrumThermalSynch_->angle_averaged(1);
  spectrumThermalSynch_->numberdensityCGS(number_density_ini);
  spectrumThermalSynch_->besselK2(bessel_K2);
  
  spectrumThermalSynch_->radiativeQ(jnu_synch_th,anu_synch_th,
                  nu_ems,nbnu);


  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){
    double jnu_tot = jnu_synch_th[ii] + jnu_synch_kappa[ii],
      anu_tot = anu_synch_th[ii] + anu_synch_kappa[ii];
    //cout << "At r,th= " << coord_ph[1] << " " << coord_ph[2] << endl;
    //cout << "in unif stuff: " << number_density << " " << nu0 << " " << thetae << " " << hypergeom << " " << jnu_tot << " " << anu_tot << " " << dsem << endl;

    // expm1 is a precise implementation of exp(x)-1
    double em1=std::expm1(-anu_tot * dsem * gg_->unitLength());
    Taunu[ii] = em1+1.;
    Inu[ii] = anu_tot == 0. ? jnu_tot * dsem * gg_->unitLength() :
      -jnu_tot / anu_tot * em1;
    
    if (Inu[ii]<0.)
      GYOTO_ERROR("In Plasmoid::radiativeQ: Inu<0");
    if (Inu[ii]!=Inu[ii] or Taunu[ii]!=Taunu[ii])
      GYOTO_ERROR("In Plasmoid::radiativeQ: Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Taunu[ii]==Taunu[ii]+1.)
      GYOTO_ERROR("In Plasmoid::radiativeQ: Inu or Taunu is infinite");
    
  }

}


double Plasmoid::numberDensity() const {
  // Converts internal cgs central enthalpy to SI
  double dens=numberDensity_cgs_;
# ifdef HAVE_UDUNITS
  dens = Units::Converter("cm-3", "m-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
        << endl ;
# endif
  return dens; }

double Plasmoid::numberDensity(string const &unit) const
{
  double dens = numberDensity();
  if (unit != "") {
# ifdef HAVE_UDUNITS
    dens = Units::Converter("m-3", unit)(dens);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
          << endl ;
# endif
  }
  return dens;
}

void Plasmoid::numberDensity(double dens) {
# ifdef HAVE_UDUNITS
  dens = Units::Converter("m-3", "cm-3")(dens);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
        << endl ;
# endif
  numberDensity_cgs_=dens;
}

void Plasmoid::numberDensity(double dens, string const &unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    dens = Units::Converter(unit, "m-3")(dens);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
          << endl ;
# endif
  }
  numberDensity(dens);
}

double Plasmoid::temperatureInitial() const {
    return temperatureInitial_; }

void Plasmoid::temperatureInitial(double tt) {
    temperatureInitial_ = tt; }

double Plasmoid::temperatureReconnection() const {
    return temperatureReconnection_; }

void Plasmoid::temperatureReconnection(double tt) {
    temperatureReconnection_ = tt; }

void Plasmoid::magnetizationParameter(double rr) {
    magnetizationParameter_=rr; }

double Plasmoid::magnetizationParameter()const{
    return magnetizationParameter_; }
  
void Plasmoid::kappaIndex(double kk) {
    kappaIndexMin_=kk; }
    
double Plasmoid::kappaIndex() const {
    return kappaIndexMin_; }

double Plasmoid::timeRef() const {
  // Converts internal M-unit time to SI
  double tt=timeRef_;
# ifdef HAVE_UDUNITS
  if (gg_)
    tt = Units::ToSeconds(tt,"geometrical_time",gg_);
  else
    GYOTO_SEVERE << "Cannot convert to seconds as metric is not set!" << endl;
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
    << endl ;
# endif
  return tt;
}

double Plasmoid::timeRef(string const &unit) const
{
  double tt = timeRef();
  if (unit != "") {
# ifdef HAVE_UDUNITS
    tt = Units::Converter("s", unit)(tt);
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
      << endl ;
# endif
  }
  return tt;
}

void Plasmoid::timeRef(double tt) {
# ifdef HAVE_UDUNITS
  tt = Units::ToGeometricalTime(tt, "s", gg_);
# else
  GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
    << endl ;
# endif
  timeRef_ = tt;
}

void Plasmoid::timeRef(double tt, string const &unit) {
  if (unit != "") {
# ifdef HAVE_UDUNITS
    if (gg_)
      tt = Units::ToSeconds(tt,unit,gg_);
  else
    GYOTO_SEVERE << "Cannot convert to seconds as metric is not set!" << endl;
# else
    GYOTO_WARNING << "Units ignored, please recompile Gyoto with --with-udunits"
      << endl ;
# endif
  }
  timeRef(tt);
}
