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
GYOTO_PROPERTY_DOUBLE(Plasmoid, RadiusMax, radiusMax,
		      "Maximun radius of the Plasmoid")
GYOTO_PROPERTY_END(Plasmoid, UniformSphere::properties)

Plasmoid::Plasmoid() :  
  UniformSphere("Plasmoid"),
  gg_(NULL),
  flag_("None"),
  numberDensity_cgs_(1.),
  temperatureReconnection_(1.),
  magnetizationParameter_(1.),
  PLIndex_(3.5),
  posSet_(false),
  posIni_(NULL),
  fourveldt_(NULL),
  radiusMax_(1.),
  varyRadius_("None"),
  spectrumPLSynchHigh_(NULL),
  spectrumPLSynchLow_(NULL),
  spectrumkappa_(NULL)
{
  kind_="Plasmoid";
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
  spectrumPLSynchHigh_ = new Spectrum::PowerLawSynchrotron();
  spectrumPLSynchLow_ = new Spectrum::PowerLawSynchrotron();
  spectrumkappa_ = new Spectrum::KappaDistributionSynchrotron();

  posIni_= new double[4];
  fourveldt_= new double[4];
}

Plasmoid::Plasmoid(const Plasmoid& orig) :
  UniformSphere(orig),
  gg_(orig.gg_),
  flag_(orig.flag_),
  numberDensity_cgs_(orig.numberDensity_cgs_),
  temperatureReconnection_(orig.temperatureReconnection_),
  magnetizationParameter_(orig.magnetizationParameter_),
  PLIndex_(orig.PLIndex_),
  posSet_(orig.posSet_),
  posIni_(NULL),
  fourveldt_(NULL),
  radiusMax_(orig.radiusMax_),
  varyRadius_(orig.varyRadius_),
  spectrumPLSynchHigh_(NULL),
  spectrumPLSynchLow_(NULL),
  spectrumkappa_(NULL)
{
  if (orig.spectrumPLSynchHigh_()) spectrumPLSynchHigh_=orig.spectrumPLSynchHigh_->clone();
  if (orig.spectrumPLSynchLow_()) spectrumPLSynchLow_=orig.spectrumPLSynchLow_->clone();
  if (orig.spectrumkappa_()) spectrumkappa_=orig.spectrumkappa_->clone();

  if(orig.posIni_){
	  posIni_= new double[4];
	  memcpy(posIni_,orig.posIni_, 4*sizeof(double));
  }

  if(orig.fourveldt_){
	  fourveldt_= new double[4];
	  memcpy(fourveldt_,orig.fourveldt_, 4*sizeof(double));
  }
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
  double t0 = posIni_[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C/60.;  // t0 in min

  double rmax_cgs = radiusMax_*GYOTO_G_OVER_C_SQUARE_CGS*gg_->mass()*1.e3;
  double vrec_cgs = 0.1*GYOTO_C_CGS*pow(magnetizationParameter_/(magnetizationParameter_+1),0.5);
  double t_inj=rmax_cgs/(vrec_cgs)/60.; //in min; //injection time, i.e. time during e- are heated and accelerated due to the reconnection, see [D. Ball et al., 2018]
  //cout << "tcur=" << tcur << ", t0=" << t0 << ", t_inj=" << t_inj << endl;

  double number_density_rec=0.; // number density of e- which follow the kappa distribution after the reconnection

  double n_dot=numberDensity_cgs_*vrec_cgs/rmax_cgs; //"Reconnection rate", see [D. Ball et al., 2020] (ie Psaltis paper)
  //cout << "n_dot=" << n_dot << endl;

  double thetae_rec = GYOTO_BOLTZMANN_CGS*temperatureReconnection_
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS); //Dimentionless temperature of the population of e- after the reconnection
  
  double BB = sqrt(4.*M_PI*magnetizationParameter_
           *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
           *numberDensity_cgs_); // Magnetic field
  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  double sigma_thomson=8.*M_PI*pow(GYOTO_ELECTRON_CLASSICAL_RADIUS_CGS,2.)/3.; // Thomson's cross section 
  double AA = (4./3.*sigma_thomson*GYOTO_C_CGS*pow(BB,2.))/(8.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS); // Coefficient of integration from [D. Ball et al., 2020] for cooling
  //cout << "AA=" << AA << ", B=" << BB << endl;

  double gamma_max = DBL_MAX;
  double gamma_change=4.*thetae_rec;
 
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
    
    gamma_max = gamma_max*pow(1.+AA*gamma_max*(tcur-(t_inj+t0)*60.),-1.);
    gamma_change = gamma_change*pow(1.+AA*gamma_change*(tcur-(t_inj+t0)*60.),-1.);

  }
  // Defining jnus, anus
  double jnu_synch_pl_low[nbnu], jnu_synch_pl_high[nbnu], jnu[nbnu];
  double anu_synch_pl_low[nbnu], anu_synch_pl_high[nbnu], anu[nbnu];
  
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    // [ exp(-anu*ds) will explose ]
    jnu_synch_pl_low[ii]=-1.;
    anu_synch_pl_low[ii]=-1.;
    jnu_synch_pl_high[ii]=-1.;
    anu_synch_pl_high[ii]=-1.;
    jnu[ii]=-1.;
    anu[ii]=-1.;
  }

  
  //PL SYNCHRO LOW
  spectrumPLSynchLow_->PLindex(-2.);
  spectrumPLSynchLow_->angle_averaged(1);
  spectrumPLSynchLow_->angle_B_pem(0.); // avg so we don't care
  spectrumPLSynchLow_->cyclotron_freq(nu0);
  spectrumPLSynchLow_->numberdensityCGS(number_density_rec);
  spectrumPLSynchLow_->gamma_min(1.);
  spectrumPLSynchLow_->gamma_max(gamma_change);
  
  spectrumPLSynchLow_->radiativeQ(jnu_synch_pl_low,anu_synch_pl_low,
                  nu_ems,nbnu);
  

  // PL SYNCHRO HIGH
  spectrumPLSynchHigh_->PLindex(PLIndex_);
  spectrumPLSynchHigh_->angle_averaged(1);
  spectrumPLSynchHigh_->angle_B_pem(0.); // avg so we don't care
  spectrumPLSynchHigh_->cyclotron_freq(nu0);
  spectrumPLSynchHigh_->numberdensityCGS(number_density_rec);
  spectrumPLSynchHigh_->gamma_min(gamma_change);
  spectrumPLSynchHigh_->gamma_max(gamma_max);
  
  spectrumPLSynchHigh_->radiativeQ(jnu_synch_pl_high,anu_synch_pl_high,
                  nu_ems,nbnu);
  
  /*
  // Kappa SYNCHRO HIGH
  double hypergeom = Gyoto::hypergeom(PLIndex_+1., thetae_rec);

  spectrumkappa_->kappaindex(PLIndex_+1.);
  spectrumkappa_->angle_averaged(1);
  spectrumkappa_->angle_B_pem(0.); // avg so we don't care
  spectrumkappa_->cyclotron_freq(nu0);
  spectrumkappa_->numberdensityCGS(number_density_rec);
  spectrumkappa_->thetae(thetae_rec);
  spectrumkappa_->hypergeometric(hypergeom);
  
  spectrumkappa_->radiativeQ(jnu,anu,
                  nu_ems,nbnu);
  */




  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){
    double jnu_tot = jnu_synch_pl_low[ii] + jnu_synch_pl_high[ii],
      anu_tot = anu_synch_pl_low[ii] + anu_synch_pl_high[ii];
    //double jnu_tot = jnu[ii], anu_tot = anu[ii];

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
  if (type=="Helical" || type=="Equatorial")
  {
    flag_=type;
  }
  else
    GYOTO_ERROR("In Plasmoid::motonType: motion not recognized, please enter a valid motion type (Helical or Equatorial)");
}

SmartPointer<Metric::Generic> Plasmoid::metric() const { return gg_; }

void Plasmoid::metric(SmartPointer<Metric::Generic> gg) {
  UniformSphere::metric(gg);
  gg_=gg;
}

void Plasmoid::initPosition(std::vector<double> const &v) {
  posIni_[0] = v[0];
  posIni_[1] = v[1];
  posIni_[2] = v[2];
  posIni_[3] = v[3];
  posSet_=true;
}

std::vector<double> Plasmoid::initPosition() const {
  std::vector<double> v (4, 0.);
  v[0] = posIni_[0];
  v[1] = posIni_[1];
  v[2] = posIni_[2];
  v[3] = posIni_[3];
  return v;
}

void Plasmoid::initVelocity(std::vector<double> const &v) {
  if (!posSet_)
  	GYOTO_ERROR("In Plasmoid::initVelocity initial Position not defined");
  fourveldt_[1] = v[0];
  fourveldt_[2] = v[1];
  fourveldt_[3] = v[2];
  fourveldt_[0] = 1.;

  double sum = 0;
  double g[4][4];

  gg_->gmunu(g, posIni_);

  for (int i=0;i<4;++i) {
    for (int j=0;j<4;++j) {
      sum+=g[i][j]*fourveldt_[i]*fourveldt_[j];
    }
  }
  if (sum>=0)
 	GYOTO_ERROR("In Plasmoid::initVelocity Initial Velocity over C");

}

std::vector<double> Plasmoid::initVelocity() const {
  std::vector<double> v (3, 0.);
  v[0] = fourveldt_[1];
  v[1] = fourveldt_[2];
  v[2] = fourveldt_[3];
  return v;
}

void Plasmoid::initCoord(std::vector<double> const &v) {
  posIni_[0] = v[0];
  posIni_[1] = v[1];
  posIni_[2] = v[2];
  posIni_[3] = v[3];
  fourveldt_[0] = v[4];
  fourveldt_[1] = v[5];
  fourveldt_[2] = v[6];
  fourveldt_[3] = v[7];
}

std::vector<double> Plasmoid::initCoord() const {
  std::vector<double> v (8, 0.);
  v[0] = posIni_[0];
  v[1] = posIni_[1];
  v[2] = posIni_[2];
  v[3] = posIni_[3];
  v[4] = fourveldt_[0];
  v[5] = fourveldt_[1];
  v[6] = fourveldt_[2];
  v[7] = fourveldt_[3];
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

void Plasmoid::radiusMax(double rr) {
	if (rr<0.2)
		GYOTO_ERROR("In Plasmoid::radiusMax radiusMax<0.2 (minimum value)");
	radiusMax_=rr;
}

double Plasmoid::radiusMax() const {
	return radiusMax_;
}

void Plasmoid::Radius(std::string vary) {
  if (vary=="Constant" || vary=="Varying") varyRadius_=vary;
  else
    GYOTO_ERROR("In Plasmoid::Radius operation on radius not recognized, please enter a valid operation (Constant or Varying)");
}

void Plasmoid::getCartesian(double const * const dates, size_t const n_dates,
          double * const x, double * const y, double * const z, 
          double * const xprime, double * const yprime, double * const zprime){
  // this yields the position of the center of the UnifSphere
  // at time t
  // fourveldt_ is the initial 3-velocity dxi/dt
  // vel is the 4-velocity dxnu/dtau

  if (n_dates!=1)
    GYOTO_ERROR("In Plasmoid::getCartesian n_dates!=1");

  if (flag_=="None")
      GYOTO_ERROR("In Plasmoid::getCartesian Motion not defined; motionType('Helical' or 'Equatorial'");

  double tt=dates[0];
  
  double r, theta, phi; // spherical coordinates
  double vel[4];
  
  if (flag_=="Helical") // Helical ejection
  {
    r = posIni_[1]+fourveldt_[1]*(tt-posIni_[0]);
    theta = posIni_[2];
    phi = posIni_[3] + posIni_[1]*posIni_[1]*fourveldt_[3]/fourveldt_[1]*(pow(posIni_[1],-1.)-pow(r,-1.)); // result of integrale of vphi over time
    //cout << phi << endl;

  }
  else // Equatorial motion (Keplerian orbit)
  {
    if (posIni_[2]!=M_PI/2.)
      cout << "Warning input theta value incompatible with 'Equatorial' motion. Theta fixed to pi/2." << endl;
    getVelocity(posIni_, vel);

    r = posIni_[1];
    theta = M_PI/2.;
    phi = posIni_[3] + vel[3]/vel[0]*(tt-posIni_[0]);

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
    GYOTO_ERROR("In Plasmoid::getVelocity Motion not defined; motionType('Helical' or 'Equatorial'");
  
  if (flag_=="Helical") // Helical case
  {
  	vel[0] = 1.;
	vel[1] = fourveldt_[1];
	vel[2] = 0.;
	vel[3] = fourveldt_[3]*pow(posIni_[1]/pos[1],2.); // conservation of the Newtonian angular momentum [Ball et al. 2020]
	gg_->normalizeFourVel(pos, vel);

  }
  else // Equatorial case
  {
    gg_->circularVelocity(pos, vel);
  }
}


int Plasmoid::Impact(Photon* ph, size_t index, Properties *data){
	// Overload function of StandardAstrobj::Impact
	// This function update the radius of the plasmoid 
	// which increase linearly during the injection phase
	// before calling the StandardAstrobj function

	double radiusMin = 0.2;
	double vrec_cgs = 0.1*GYOTO_C_CGS*pow(magnetizationParameter_/(magnetizationParameter_+1),0.5);
  double t_inj=radiusMax_*GYOTO_G_OVER_C_SQUARE_CGS*gg_->mass()*1.e3/(vrec_cgs)/60.; //in min;
  double t0 = posIni_[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C/60.;  // t0 in min

  size_t sz = ph -> parallelTransport()?16:8;
  state_t p1(sz);
  ph->getCoord(index, p1);
  double tcur = p1[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C/60.; //tcur in min

  
  if (varyRadius_== "Varying")
  {
    if (tcur<=t0) radius(radiusMin);
    else if (tcur<=t0+t_inj) radius(radiusMin+(radiusMax_-radiusMin)*(tcur-t0)/t_inj);
    else radius(radiusMax_);
  }
	else if (varyRadius_== "Constant") radius(radiusMax_);
  else{
    GYOTO_ERROR("In Plasmoid::Impact operation on radius not recognized. Use Radius('Constant' or 'Varying')");
  }

	return Standard::Impact(ph, index, data);
}















void Plasmoid::radiativeQ(double Inu[], double Qnu[], double Unu[], double Vnu[], // output
              double alphaInu[], double alphaQnu[], double alphaUnu[], double alphaVnu[], // output
              double rQnu[], double rUnu[], double rVnu[], // outut
              double const nu_ems[], size_t nbnu, double dsem,
              state_t const &coord_ph, double const coord_obj[8]) const {

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  double tcur=coord_ph[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C/60.; // in min
  double t0 = posIni_[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C/60.;  // t0 in min

  double rmax_cgs = radiusMax_*GYOTO_G_OVER_C_SQUARE_CGS*gg_->mass()*1.e3;
  double vrec_cgs = 0.1*GYOTO_C_CGS*pow(magnetizationParameter_/(magnetizationParameter_+1),0.5);
  double t_inj=rmax_cgs/(vrec_cgs)/60.; //in min; //injection time, i.e. time during e- are heated and accelerated due to the reconnection, see [D. Ball et al., 2018]
  //cout << "tcur=" << tcur << ", t0=" << t0 << ", t_inj=" << t_inj << endl;

  double number_density_rec=0.; // number density of e- which follow the kappa distribution after the reconnection

  double n_dot=numberDensity_cgs_*vrec_cgs/rmax_cgs; //"Reconnection rate", see [D. Ball et al., 2020] (ie Psaltis paper)
  //cout << "n_dot=" << n_dot << endl;

  double tempRec=temperatureReconnection_;
  double thetae_rec = GYOTO_BOLTZMANN_CGS*temperatureReconnection_
    /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS); //Dimentionless temperature of the population of e- after the reconnection
  if (thetae_rec<20.)
    GYOTO_ERROR("In radiativeQ : Theta too low, the model of the distribution is no valid, please increase temperatureReconnection(double)");

  
  double BB = sqrt(4.*M_PI*magnetizationParameter_
           *GYOTO_PROTON_MASS_CGS * GYOTO_C_CGS * GYOTO_C_CGS
           *numberDensity_cgs_); // Magnetic field
  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

  double sigma_thomson=8.*M_PI*pow(GYOTO_ELECTRON_CLASSICAL_RADIUS_CGS,2.)/3.; // Thomson's cross section 
  double AA = (4./3.*sigma_thomson*GYOTO_C_CGS*pow(BB,2.))/(8.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS); // Coefficient of integration from [D. Ball et al., 2020] for cooling
  //cout << "AA=" << AA << ", B=" << BB << endl;

  double gamma_max = DBL_MAX;
  double gamma_change=2.*thetae_rec;
 
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

    gamma_max = gamma_max*pow(1.+AA*gamma_max*(tcur-(t_inj+t0)*60.),-1.);
    gamma_change = gamma_change*pow(1.+AA*gamma_change*(tcur-(t_inj+t0)*60.),-1.);
  }
  //cout << radius() << endl;

  // Defining jnus, anus
  double jInu[nbnu], jQnu[nbnu], jUnu[nbnu], jVnu[nbnu];
  double aInu[nbnu], aQnu[nbnu], aUnu[nbnu], aVnu[nbnu];
  double rotQnu[nbnu], rotUnu[nbnu], rotVnu[nbnu];
  
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    // [ exp(-anu*ds) will explose ]
    jInu[ii]=-1.;
    jQnu[ii]=-1.;
    jUnu[ii]=-1.;
    jVnu[ii]=-1.;
    aInu[ii]=-1.;
    aQnu[ii]=-1.;
    aUnu[ii]=-1.;
    aVnu[ii]=-1.;
    rotQnu[ii]=-1.;
    rotUnu[ii]=-1.;
    rotVnu[ii]=-1.;
  }
  
  

  // PL SYNCHRO HIGH
  spectrumPLSynchHigh_->PLindex(PLIndex_);
  spectrumPLSynchHigh_->angle_averaged(1);
  spectrumPLSynchHigh_->angle_B_pem(0.); // avg so we don't care
  spectrumPLSynchHigh_->cyclotron_freq(nu0);
  spectrumPLSynchHigh_->numberdensityCGS(0.87*number_density_rec);
  spectrumPLSynchHigh_->gamma_min(gamma_change);
  spectrumPLSynchHigh_->gamma_max(gamma_max);

  //spectrumThermalSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu,
  //                                  aInu, aQnu, aUnu, aVnu,
  //                                  rotQnu, rotUnu, rotVnu, nu_ems,nbnu);


  double angle_B = spectrumPLSynchHigh_->angle_B_pem();
  double Xhi=0.; // Express Xhi in function of angle_B (angle_B = theta in the VadeMecum)

  for (size_t ii=0; ii<nbnu; ++ii){
    double jnu_tot[4], anu_tot[4], rnu_tot[3], Onu[4][4];
    jnu_tot[0] = jInu[ii];
    jnu_tot[1] = jQnu[ii];
    jnu_tot[2] = jUnu[ii];
    jnu_tot[3] = jVnu[ii];
    anu_tot[0] = aInu[ii];
    anu_tot[1] = aQnu[ii];
    anu_tot[2] = aUnu[ii];
    anu_tot[3] = aVnu[ii];
    rnu_tot[0] = rotQnu[ii];
    rnu_tot[1] = rotUnu[ii];
    rnu_tot[2] = rotVnu[ii];

    Omatrix(Onu, anu_tot, rnu_tot, Xhi, dsem);

    // matrix-vector product of Onu x jnu_tot
    double increment[4]={0};
    for (int ll=0;ll<4;ll++){
      for (int cc=0;cc<4;cc++){
        increment[ll] += Onu[ll][cc]*jnu_tot[cc];
      }
    }

    // WRINTING VALUES IN OUTPUTS ARRAYS
    Inu[ii] = increment[0]*dsem*gg_->unitLength();
    Qnu[ii] = increment[1]*dsem*gg_->unitLength();
    Unu[ii] = increment[2]*dsem*gg_->unitLength();
    Vnu[ii] = increment[3]*dsem*gg_->unitLength();

    alphaInu[ii] = anu_tot[0];
    alphaQnu[ii] = anu_tot[1];
    alphaUnu[ii] = anu_tot[2];
    alphaVnu[ii] = anu_tot[3];

    rQnu[ii] = rnu_tot[0];
    rUnu[ii] = rnu_tot[1];
    rVnu[ii] = rnu_tot[2];

    // Trouver une securite !!!!
  }
}

void Plasmoid::Omatrix(double Onu[4][4], double alphanu[4], double rnu[3], double Xhi, double dsem) const{
  double alphasqrt, rsqrt, lambda1, lambda2, Theta, sigma;
  
  double aI=alphanu[0], aV=alphanu[3];
  double aQ=alphanu[1]*cos(2*Xhi)+alphanu[2]*sin(2*Xhi);
  double aU=alphanu[2]*cos(2*Xhi)-alphanu[1]*sin(2*Xhi);
  double rQ=rnu[0]*cos(2*Xhi)+rnu[1]*sin(2*Xhi);
  double rU=rnu[1]*cos(2*Xhi)-rnu[0]*sin(2*Xhi);
  double rV=rnu[2];

  alphasqrt = aQ*aQ+aU*aU+aV*aV;
  rsqrt = rQ*rQ+rU*rU+rV*rV;
  lambda1 = pow(pow(pow(alphasqrt-rsqrt,2.)/4.+pow(aQ*rQ+aU*rU+aV*rV,2.),0.5)+pow(alphasqrt-rsqrt,2.)/2.,0.5);
  lambda2 = pow(pow(pow(alphasqrt-rsqrt,2.)/4.+pow(aQ*rQ+aU*rU+aV*rV,2.),0.5)-pow(alphasqrt-rsqrt,2.)/2.,0.5);
  Theta = pow(lambda1,2)+pow(lambda2,2);
  sigma = (aQ*rQ+aU*rU+aV*rV)/abs(aQ*rQ+aU*rU+aV*rV);

  double M1[4][4]={0}, M2[4][4]={0}, M3[4][4]={0}, M4[4][4]={0};

  // Fill of M1
  for (int ii=0;ii<4;ii++){
    M1[ii][ii]=1.;
  }

  // Fill of M2
  M2[0][1]=lambda2*aQ-sigma*lambda1*rQ;
  M2[0][2]=lambda2*aU-sigma*lambda1*rU;
  M2[0][3]=lambda2*aV-sigma*lambda1*rV;
  M2[1][0]=M2[0][1];
  M2[1][2]=sigma*lambda1*aV+lambda2*rV;
  M2[1][3]=-sigma*lambda1*aU-lambda2*rU;
  M2[2][0]=M2[0][2];
  M2[2][1]=-M2[1][2];
  M2[2][3]=sigma*lambda1*aQ+lambda2*rQ;
  M2[3][0]=M2[0][3];
  M2[3][1]=-M2[1][3];
  M2[3][2]=-M2[2][3];

  // Fill of M3
  M3[0][1]=lambda1*aQ+sigma*lambda2*rQ;
  M3[0][2]=lambda1*aU+sigma*lambda2*rU;
  M3[0][3]=lambda1*aV+sigma*lambda2*rV;
  M3[1][0]=M3[0][1];
  M3[1][2]=-sigma*lambda2*aV+lambda1*rV;
  M3[1][3]=sigma*lambda2*aU-lambda1*rU;
  M3[2][0]=M3[0][2];
  M3[2][1]=-M3[1][2];
  M3[2][3]=-sigma*lambda2*aQ+lambda1*rQ;
  M3[3][0]=M3[0][3];
  M3[3][1]=-M3[1][3];
  M3[3][2]=-M3[2][3];

  // Fill of M4
  M4[0][0]= (alphasqrt+rsqrt)/2.;
  M4[0][1]=aV*rU-aU*rV;
  M4[0][2]=aQ*rV-aV*rQ;
  M4[0][3]=aU*rQ-aQ*rU;
  M4[1][0]=-M4[0][1];
  M4[1][1]=pow(aQ,2)+pow(rQ,2)-(alphasqrt+rsqrt)/2.;
  M4[1][2]=aQ*aU+rQ*rU;
  M4[1][3]=aV*aQ+rV*rQ;
  M4[2][0]=-M4[0][2];
  M4[2][1]=M4[1][2];
  M4[2][2]=pow(aU,2)+pow(rU,2)-(alphasqrt+rsqrt)/2.;
  M4[2][3]=aU*aV+rU*rV;
  M4[3][0]=-M4[0][3];
  M4[3][1]=M4[1][3];
  M4[3][2]=M4[2][3];
  M4[3][3]=pow(aV,2)+pow(rV,2)-(alphasqrt+rsqrt)/2.;

  // Filling O matrix, output
  for (int ii=0;ii<4;ii++){
    for (int jj=0;jj<4;jj++){
      Onu[ii][jj]=exp(-aI*dsem)*(\
        (cosh(lambda1*dsem)+cos(lambda2*dsem))*M1[ii][jj]/2. \
        -sin(lambda2*dsem)*M2[ii][jj]/Theta \
        -sinh(lambda1*dsem)*M3[ii][jj]/Theta \
        +(cosh(lambda1*dsem)-cos(lambda2*dsem))*M4[ii][jj]/Theta);
    }
  }

}