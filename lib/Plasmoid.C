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
GYOTO_PROPERTY_DOUBLE(Plasmoid, TemperatureInitial, temperatureInitial,
              "temperature initial, constant through initial population of e-")
GYOTO_PROPERTY_DOUBLE(Plasmoid, TemperatureReconnection, temperatureReconnection,
              "temperature reconnection, constant through reconnected population of e-")
GYOTO_PROPERTY_DOUBLE(Plasmoid, MagnetizationParameter,
              magnetizationParameter,
              "magnetization parameter")
GYOTO_PROPERTY_DOUBLE(Plasmoid, PowerIndex, powerIndex,
              "PL index")
GYOTO_PROPERTY_END(Plasmoid, Star::properties)

Plasmoid::Plasmoid() :
  Star(),
  numberDensity_cgs_(1.),
  temperature_ini_(1.),
  temperature_rec_(1.),
  magnetizationParameter(1.),
  powerIndex(1.),
  spectrumPowerLawSynch_(NULL)
  {
      kind_="Plasmoid";
  # ifdef GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "done." << endl;
  # endif
      spectrumPowerLawSynch_ = new Spectrum::PowerLawSynchrotronSpectrum();
  }
  spectrumLowThermalSynch_(NULL)
  {
      kind_="Plasmoid";
  # ifdef GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "done." << endl;
  # endif
      spectrumLowThermalSynch_ = new Spectrum::ThermalSynchrotronSpectrum();
  }
  spectrumHighThermalSynch_(NULL)
  {
      kind_="Plasmoid";
  # ifdef GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "done." << endl;
  # endif
      spectrumHighThermalSynch_ = new Spectrum::ThermalSynchrotronSpectrum();
  }

Plasmoid::Plasmoid(const Plasmoid& orig) :
  Star(orig),
  numberDensity_cgs_(orig.numberDensity_cgs_),
  temperature_ini_(orig.temperature_ini_),
  temperature_rec_(orig.temperature_rec_),
  powerIndex_(orig.powerIndex_),
  magnetizationParameter_(orig.magnetizationParameter_),
  spectrumPowerLawSynch_(NULL)
{
  if (orig.spectrumPowerLawSynch_()) spectrumPowerLawSynch_=orig.spectrumPowerLawSynch_->clone();
}
  spectrumLowThermalSynch_(NULL)
{
  if (orig.spectrumLowThermalSynch_()) spectrumLowThermalSynch_=orig.spectrumLowThermalSynch_->clone();
}
  spectrumHighThermalSynch_(NULL)
{
  if (orig.spectrumHighThermalSynch_()) spectrumHighThermalSynch_=orig.spectrumHighThermalSynch_->clone();
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
  double t_inj=radius()/pow(magnetizationParameter_/(magnetizationParameter_+1),0.5); //injection time, i.e. time during e- are heated and accelerated due to the reconnection, see [D. Ball et al., 2018]
  
  double number_density_ini=numberDensity_cgs_; // number density of e- which follow the initial thermal distribution, before the reconnection
  double number_density_rec=0.; // umber density of e- which follow the kappa distribution after the reconnection

  double n_dot=0.1*numberDensity_cgs_*pow(magnetizationParameter_/(magnetizationParameter_+1),0.5)/radius(); //Reconnection rate, see [D. Ball et al., 2020] (ie Psaltis paper)

  double thetae_ini = GYOTO_BOLTZMANN_CGS*temperature_ini_
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS); //Dimentionless temperature of the initial thermal population of e-
  double thetae_rec = GYOTO_BOLTZMANN_CGS*temperature_rec_
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS); //Dimentionless temperature of the population of e- after the reconnection

  double ELECTRON_RADIUS_CGS = 2.8179403227e-13 // Classical radius of the electron in cgs
  
  double BB = sqrt(8.*M_PI*magnetizationParameter_
           *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
           *number_density_rec); // Magnetic field
  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  // UPDATE NUMBER DENSITIES AND THETA VALUES
  if (tcur<0.)
  {
    number_density_ini=numberDensity_cgs_;
    number_density_rec=0.;
  }
  else if (tcur<=t_inj)
  {
    number_density_rec=n_dot*tcur;
    number_density_ini=numberDensity_cgs_ - number_density_rec;
  }
  else
  {
    number_density_rec=n_dot*tinj;
    number_density_ini=numberDensity_cgs_ - number_density_rec;

    double sigma_thomson=8.*M_PI*pow(ELECTRON_RADIUS_CGS,2)/3.
    double AA=(4./3.*sigma_thomson*GYOTO_C_CGS*pow(BB,2)/8.*M_PI)/(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS); // Coefficient of integration from [D. Ball et al., 2020]
    thetae_rec=max(thetae_rec*pow(1+AA*3.*thetae_rec*t,-1),thetae_ini); //see [Radiative Process in Astrophysics, G. Rybicki and Lightman, 1986];theta=gamma/3
  }
  
  double bessel_K2_ini = Gyoto::bessk(2,thetae_ini);
  double bessel_K2_rec = Gyoto::bessk(2,thetae_rec);
  
  // Defining jnus, anus
  double jnu_synch_pl[nbnu], anu_synch_pl[nbnu];
  double jnu_synch_th_ini[nbnu], anu_synch_th_ini[nbnu];
  double jnu_synch_th_rec[nbnu], anu_synch_th_rec[nbnu];
  
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    // [ exp(-anu*ds) will explose ]
    jnu_synch_pl[ii]=-1.;
    anu_synch_pl[ii]=-1.;
    jnu_synch_th_ini[ii]=-1.;
    anu_synch_th_ini[ii]=-1.;
    jnu_synch_th_rec[ii]=-1.;
    anu_synch_th_rec[ii]=-1.;
  }
  
  // POWER LAW SYNCHRO
  spectrumPowerLawSynch_->PLindex_(powerIndex_);
  spectrumPowerLawSynch_->numberdensityCGS(number_density_rec);
  spectrumPowerLawSynch_->angle_averaged(1);
  spectrumPowerLawSynch_->angle_B_pem(0.); // avg so we don't care
  spectrumPowerLawSynch_->cyclotron_freq(nu0);

  spectrumPowerLawSynch_->radiativeQ(jnu_synch_pl,anu_synch_pl, 3.*thetae_rec,
                  nu_ems,nbnu);

  // THERMAL INI SYNCHRO
  spectrumLowThermalSynch_->T_(thetae_ini);
  spectrumLowThermalSynch_->numberdensityCGS_(number_density_ini);
  spectrumLowThermalSynch_->angle_B_pem(0.); // avg so we don't care
  spectrumLowThermalSynch_->cyclotron_freq_(nu0);
  spectrumLowThermalSynch_->angle_averaged(1);
  spectrumLowThermalSynch_->bessel_K2(bessel_K2_ini);

  spectrumLowThermalSynch_->radiativeQ(jnu_synch_th_ini,anu_synch_th_ini,
                  nu_ems,nbnu);

  // THERMAL REC SYNCHRO
  spectrumHighThermalSynch_->T_(thetae_rec);
  spectrumHighThermalSynch_->numberdensityCGS_(number_density_rec);
  spectrumHighThermalSynch_->angle_B_pem(0.); // avg so we don't care
  spectrumHighThermalSynch_->cyclotron_freq_(nu0);
  spectrumHighThermalSynch_->angle_averaged(1);
  spectrumHighThermalSynch_->bessel_K2(bessel_K2_rec);

  spectrumHighThermalSynch_->radiativeQ(jnu_synch_th_rec,anu_synch_th_rec,
                  nu_ems,nbnu);


  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){
    double jnu_tot = jnu_synch_th_ini[ii] + jnu_synch_th_rec[ii]+ jnu_synch_pl[ii],
      anu_tot = anu_synch_th_ini[ii] + anu_synch_th_rec[ii] + anu_synch_kappa[ii];
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

double Plasmoid::temperature_ini() const {
    return temperature_; }

void Plasmoid::temperature_ini(double tt) {
    temperature_ini_ = tt; }

double Plasmoid::temperature_rec() const {
    return temperature_; }

void Plasmoid::temperature_rec(double tt) {
    temperature_rec_ = tt; }

void Plasmoid::magnetizationParameter(double rr) {
  magnetizationParameter_=rr;}

double Plasmoid::magnetizationParameter()const{
  return magnetizationParameter_;}

double Plasmoid::powerIndex() const {
    return powerIndex_; }

void Plasmoid::powerIndex(double ind) {
    powerIndex_ = ind; }
