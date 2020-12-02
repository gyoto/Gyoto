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
#include "GyotoWorldline.h"
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
GYOTO_PROPERTY_VECTOR_DOUBLE(Plasmoid, InitPosition, initPosition,
           "(t,r,theta,phi) initial position of plasmoid")
GYOTO_PROPERTY_VECTOR_DOUBLE(Plasmoid, InitVelocity, initVelocity,
           "(dr/dt,dtheta/dt,dphi/dt) initial 3-velocity "
           "of plasmoid")
GYOTO_PROPERTY_DOUBLE_UNIT(Plasmoid, NumberDensity, numberDensity,
               "cgs number density, constant through plasmoid")
GYOTO_PROPERTY_DOUBLE(Plasmoid, TemperatureReconnection, temperatureReconnection,
               "Temperature de reconnection")
GYOTO_PROPERTY_DOUBLE(Plasmoid, MagnetizationParameter,
              magnetizationParameter,
              "magnetization parameter")
GYOTO_PROPERTY_DOUBLE(Plasmoid, PLIndex, PLIndex,
		      "PL index of kappa-synchrotron")
GYOTO_PROPERTY_END(Plasmoid, UniformSphere::properties)

Plasmoid::Plasmoid() :
  UniformSphere("Plasmoid"),
  gg_(NULL),
  flag_("None"),
  numberDensity_cgs_(1.),
  temperatureReconnection_(1.),
  magnetizationParameter_(1.),
  PLIndex_(3.5),
  posSet(false),
  spectrumThermalSynch_(NULL)
  //spectrumPLSynch_(NULL)
{
  kind_="Plasmoid";
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
  //spectrumPLSynch_ = new Spectrum::PowerLawSynchrotron();
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
}

Plasmoid::Plasmoid(const Plasmoid& orig) :
  UniformSphere(orig),
  gg_(orig.gg_),
  flag_(orig.flag_),
  numberDensity_cgs_(orig.numberDensity_cgs_),
  temperatureReconnection_(orig.temperatureReconnection_),
  magnetizationParameter_(orig.magnetizationParameter_),
  PLIndex_(orig.PLIndex_),
  posSet(orig.posSet),
  //spectrumPLSynch_(NULL)
  spectrumThermalSynch_(NULL)
{
  //if (orig.spectrumPLSynch_()) spectrumPLSynch_=orig.spectrumPLSynch_->clone();
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
  double tcur=coord_ph[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C/60.; // in min
  double t0 = pos_[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C/60.;  // t0 in min

  double vrec_cgs = 0.1*GYOTO_C_CGS*pow(magnetizationParameter_/(magnetizationParameter_+1),0.5);
  double t_inj=radius("cm")/vrec_cgs/60.; //in min; //injection time, i.e. time during e- are heated and accelerated due to the reconnection, see [D. Ball et al., 2018]
  //cout << "tcur=" << tcur << ", t0=" << t0 << ", t_inj=" << t_inj << endl;

  //double number_density_ini=numberDensity_cgs_; // number density of e- which follow the initial thermal distribution, before the reconnection
  double number_density_rec=0.; // number density of e- which follow the kappa distribution after the reconnection

  double n_dot=numberDensity_cgs_*vrec_cgs/radius("cm"); //"Reconnection rate", see [D. Ball et al., 2020] (ie Psaltis paper)
  //cout << "n_dot=" << n_dot << endl;

  double tempRec=temperatureReconnection_;
  double thetae_rec = GYOTO_BOLTZMANN_CGS*temperatureReconnection_
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS); //Dimentionless temperature of the population of e- after the reconnection
  double gamma_min=3.*thetae_rec;

  
  double BB = sqrt(4.*M_PI*magnetizationParameter_
           *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
           *numberDensity_cgs_); // Magnetic field
  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  double sigma_thomson=8.*M_PI*pow(GYOTO_ELECTRON_CLASSICAL_RADIUS_CGS,2.)/3.; // Thomson's cross section 
  double AA = (4./3.*sigma_thomson*GYOTO_C_CGS*pow(BB,2.))/(8.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS); // Coefficient of integration from [D. Ball et al., 2020] for cooling
  //cout << "AA=" << AA << ", B=" << BB << endl;

  double gamma_max = DBL_MAX;
  double gamma_max_0 = DBL_MAX;
 
  // COMPUTE VALUES IN FUNCTION OF PHASE
  if (tcur<=t0)
  {
    number_density_rec=0.;
  }
  else if (tcur<=t0+t_inj) // HEATING TIME
  {
    number_density_rec=n_dot*(tcur-t0)*60.;
  }
  else // COOLING TIME
  {
    // evolution of the number densities
    number_density_rec=n_dot*t_inj*60.;
    
    thetae_rec=thetae_rec*pow(1+AA*thetae_rec*(tcur-(t_inj+t0))*60.,-1.);
    tempRec=thetae_rec*GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS/GYOTO_BOLTZMANN_CGS;

    /*
    gamma_max_0 = 1.e10; //pow(DBL_MIN*(4.*M_PI*(1.-pow(DBL_MAX,1.-(kappaIndex-1.))))/((kappaIndex-1.)-1.),-1./(kappaIndex-1.));
    gamma_max = max(gamma_max_0*pow(1+AA*gamma_max_0*(tcur-(t_inj+t0)),-1.),gamma_min);
    */
  }
  
  //cout << "n/n0=" << number_density_rec/numberDensity_cgs_ << endl;
  //cout << "thetae=" << thetae_rec <<endl;

  // Defining jnus, anus
  double jnu_synch[nbnu];
  double anu_synch[nbnu];
  
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    // [ exp(-anu*ds) will explose ]
    jnu_synch[ii]=-1.;
    anu_synch[ii]=-1.;
  }

  /*
  //PL SYNCHRO
  spectrumPLSynch_->PLindex(PLIndex_);
  spectrumPLSynch_->angle_averaged(1);
  spectrumPLSynch_->angle_B_pem(0.); // avg so we don't care
  spectrumPLSynch_->cyclotron_freq(nu0);
  spectrumPLSynch_->numberdensityCGS(number_density_rec);
  spectrumPLSynch_->gamma_min(gamma_min);
  spectrumPLSynch_->gamma_max(gamma_max);
  
  spectrumPLSynch_->radiativeQ(jnu_synch,anu_synch,
                  nu_ems,nbnu);
  */
  
  double besselK2 = Gyoto::bessk(2,thetae_rec);

  // THERMAL SYNCHRO
  spectrumThermalSynch_->temperature(tempRec);
  spectrumThermalSynch_->numberdensityCGS(number_density_rec);
  spectrumThermalSynch_->angle_B_pem(0.);
  spectrumThermalSynch_->cyclotron_freq(nu0);
  spectrumThermalSynch_->angle_averaged(1);
  spectrumThermalSynch_->besselK2(besselK2);

  spectrumThermalSynch_->radiativeQ(jnu_synch,anu_synch,
                        nu_ems,nbnu);

  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){
    double jnu_tot = jnu_synch[ii],
      anu_tot = anu_synch[ii];
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

void Plasmoid::motionType(std::string const type){
  if (type=="Helicoidal" || type=="Equatorial")
  {
    flag_=type;
  }
  else
    GYOTO_ERROR("In Plasmoid::motonType: motion not recognized, please enter a valid motion type (Helicoidal or Equatorial)");
}

SmartPointer<Metric::Generic> Plasmoid::metric() const { return gg_; }

void Plasmoid::metric(SmartPointer<Metric::Generic> gg) {
  UniformSphere::metric(gg);
  gg_=gg;
}

void Plasmoid::initPosition(std::vector<double> const &v) {
  pos_[0] = v[0];
  pos_[1] = v[1];
  pos_[2] = v[2];
  pos_[3] = v[3];
  posSet=true;
}

std::vector<double> Plasmoid::initPosition() const {
  std::vector<double> v (4, 0.);
  v[0] = pos_[0];
  v[1] = pos_[1];
  v[2] = pos_[2];
  v[3] = pos_[3];
  return v;
}

void Plasmoid::initVelocity(std::vector<double> const &v) {
  if (!posSet)
  	GYOTO_ERROR("In Plasmoid::initVelocity initial Position not defined");
  vel_[1] = v[0];
  vel_[2] = v[1];
  vel_[3] = v[2];
  vel_[0] = 1.;

  double sum = 0;
  double g[4][4];

  gg_->gmunu(g, pos_);

  for (int i=0;i<4;++i) {
    for (int j=0;j<4;++j) {
      sum+=g[i][j]*vel_[i]*vel_[j];
    }
  }
  if (sum>=0)
 	GYOTO_ERROR("In Plasmoid::initVelocity Initial Velocity over C");

}

std::vector<double> Plasmoid::initVelocity() const {
  std::vector<double> v (3, 0.);
  v[0] = vel_[1];
  v[1] = vel_[2];
  v[2] = vel_[3];
  return v;
}

void Plasmoid::initCoord(std::vector<double> const &v) {
  pos_[0] = v[0];
  pos_[1] = v[1];
  pos_[2] = v[2];
  pos_[3] = v[3];
  vel_[0] = v[4];
  vel_[1] = v[5];
  vel_[2] = v[6];
  vel_[3] = v[7];
}

std::vector<double> Plasmoid::initCoord() const {
  std::vector<double> v (8, 0.);
  v[0] = pos_[0];
  v[1] = pos_[1];
  v[2] = pos_[2];
  v[3] = pos_[3];
  v[4] = vel_[0];
  v[5] = vel_[1];
  v[6] = vel_[2];
  v[7] = vel_[3];
  return v;
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

double Plasmoid::temperatureReconnection() const {
    return temperatureReconnection_;}

void Plasmoid::temperatureReconnection(double tt) {
    temperatureReconnection_=tt;}
  
void Plasmoid::magnetizationParameter(double rr) {
    magnetizationParameter_=rr; }

double Plasmoid::magnetizationParameter()const{
    return magnetizationParameter_; }
  
void Plasmoid::PLIndex(double kk) {
    PLIndex_=kk; }
    
double Plasmoid::PLIndex() const {
    return PLIndex_; }

void Plasmoid::getCartesian(double const * const dates, size_t const n_dates,
          double * const x, double * const y, double * const z, 
          double * const xprime, double * const yprime, double * const zprime){
  // this yields the position of the center of the UnifSphere
  // at time t
  // vel_ is the initial 3-velocity dxi/dt
  // vel is the 4-velocity dxnu/dtau

  if (n_dates!=1)
    GYOTO_ERROR("In Plasmoid::getCartesian n_dates!=1");

  if (flag_=="None")
      GYOTO_ERROR("In Plasmoid::getCartesian Motion not defined; motionType('Helicoidal' or 'Equatorial'");

  double tt=dates[0];
  
  double r, theta, phi; // spherical coordinates
  double vel[4];
  
  if (flag_=="Helicoidal") // Helicoidal ejection
  {
    r = pos_[1]+vel_[1]*(tt-pos_[0]);
    theta = pos_[2];
    phi = pos_[3] + pos_[1]*pos_[1]*vel_[3]/vel_[1]*(pow(pos_[1],-1.)-pow(r,-1.)); // result of integrale of vphi over time
    //cout << phi << endl;

  }
  else // Equatorial motion (Keplerian orbit)
  {
    if (pos_[2]!=M_PI/2.)
      cout << "Warning input theta value incompatible with 'Equatorial' motion. Theta fixed to pi/2." << endl;
    getVelocity(pos_, vel);

    r = pos_[1];
    theta = M_PI/2.;
    phi = pos_[3] + vel[3]/vel[0]*(tt-pos_[0]);

  }
  // Convertion into cartesian coordinates
  x[0] = r*sin(theta)*cos(phi);
  y[0] = r*sin(theta)*sin(phi);
  z[0] = r*cos(theta);

  if (xprime!=NULL && yprime!=NULL && zprime!=NULL)
  {
    xprime[0] = r*sin(theta)*sin(phi)*vel[2];
    yprime[0] = -r*sin(theta)*cos(phi)*vel[2];
    zprime[0] = 0.;
  }
}

void Plasmoid::getVelocity(double const pos[4], double vel[4]){
  if (!gg_)
    GYOTO_ERROR("In Plasmoid::getVelocity Metric not set");
  if (flag_=="None")
    GYOTO_ERROR("In Plasmoid::getVelocity Motion not defined; motionType('Helicoidal' or 'Equatorial'");
  
  if (flag_=="Helicoidal") // Helicoidal case
  {
  	vel[0] = 1.;
	vel[1] = vel_[1];
	vel[2] = 0.;
	vel[3] = vel_[3]*pow(pos_[1]/pos[1],2.); // conservation of the Newtonian angular momentum [Ball et al. 2020]
	gg_->normalizeFourVel(pos, vel);

  }
  else // Equatorial case
  {
    gg_->circularVelocity(pos, vel);
  }
}