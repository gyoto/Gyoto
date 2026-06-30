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

#ifdef GYOTO_USE_CFITSIO
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); GYOTO_ERROR(ermsg); }
#endif

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
GYOTO_PROPERTY_DOUBLE(Plasmoid, RadiusMax, radiusMax,
		          "Maximun radius of the Plasmoid")
GYOTO_PROPERTY_END(Plasmoid, UniformSphere::properties)

Plasmoid::Plasmoid() : 
  FitsRW(), 
  UniformSphere("Plasmoid"),
  posIni_(NULL),
  fourveldt_(NULL),
  flag_("None"),
  posSet_(false),
  t_inj_(1.),
  radiusMax_(1.),
  varyRadius_("None"),
  beta_(1.),
  filename_("None"),
  emis_polar_array_(NULL),
  abs_polar_array_(NULL),
  rot_polar_array_(NULL),
  freq_array_(NULL),
  time_array_(NULL),
  angle_array_(NULL),
  nb_time_(0.),
  nb_freq_(0.),
  nb_angle_(0.)
{
  kind_="Plasmoid";
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif

  posIni_= new double[4];
  fourveldt_= new double[4];

  emis_polar_array_ = new double*[4];
  abs_polar_array_  = new double*[4];
  rot_polar_array_  = new double*[3];
  
  int ncells_coefs = nb_time_*nb_freq_*nb_angle_;
  for (int ii=0; ii<4; ii++){
    emis_polar_array_[ii] = new double[ncells_coefs];
    abs_polar_array_[ii] = new double[ncells_coefs];
    if (ii<3)
      rot_polar_array_[ii] = new double[ncells_coefs];
  }

  spectrumKappaSynch_ = new Spectrum::KappaDistributionSynchrotron();
}

Plasmoid::Plasmoid(const Plasmoid& orig) :
  FitsRW(orig),
  UniformSphere(orig),
  posIni_(NULL),
  fourveldt_(NULL),
  flag_(orig.flag_),
  posSet_(orig.posSet_),
  t_inj_(orig.t_inj_),
  radiusMax_(orig.radiusMax_),
  varyRadius_(orig.varyRadius_),
  beta_(orig.beta_),
  filename_(orig.filename_),
  emis_polar_array_(NULL),
  abs_polar_array_(NULL),
  rot_polar_array_(NULL),
  freq_array_(NULL),
  time_array_(NULL),
  angle_array_(NULL),
  nb_time_(orig.nb_time_),
  nb_freq_(orig.nb_freq_),
  nb_angle_(orig.nb_angle_)
{

  if(orig.posIni_){
	  posIni_= new double[4];
	  memcpy(posIni_,orig.posIni_, 4*sizeof(double));
  }

  if(orig.fourveldt_){
	  fourveldt_= new double[4];
	  memcpy(fourveldt_,orig.fourveldt_, 4*sizeof(double));
  }
  
  int ncells_coefs = nb_time_*nb_freq_*nb_angle_;
  if (orig.emis_polar_array_){
    emis_polar_array_ = new double*[4];
    for (int ii=0; ii<4; ii++){
      emis_polar_array_[ii] = new double[ncells_coefs];
      memcpy(emis_polar_array_[ii],orig.emis_polar_array_[ii], ncells_coefs*sizeof(double));
    }
  }

  if (orig.abs_polar_array_){
    abs_polar_array_ = new double*[4];
    for (int ii=0; ii<4; ii++){
      abs_polar_array_[ii] = new double[ncells_coefs];
      memcpy(abs_polar_array_[ii],orig.abs_polar_array_[ii], ncells_coefs*sizeof(double));
    }
  }

  if (orig.rot_polar_array_){
    rot_polar_array_ = new double*[3];
    for (int ii=0; ii<3; ii++){
      rot_polar_array_[ii] = new double[ncells_coefs];
      memcpy(rot_polar_array_[ii],orig.rot_polar_array_[ii], ncells_coefs*sizeof(double));
    }
  }

  if(orig.time_array_){
    time_array_= new double[nb_time_];
    memcpy(time_array_,orig.time_array_, nb_time_*sizeof(double));
  }

  if(orig.freq_array_){
    freq_array_= new double[nb_freq_];
    memcpy(freq_array_,orig.freq_array_, nb_freq_*sizeof(double));
  }

  if(orig.angle_array_){
    angle_array_= new double[nb_angle_];
    memcpy(angle_array_,orig.angle_array_, nb_angle_*sizeof(double));
  }

  if (orig.spectrumKappaSynch_()) spectrumKappaSynch_=orig.spectrumKappaSynch_->clone();
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
  if (filename_=="None")
      GYOTO_ERROR("In Plamsoid RadiativeQ : filename_ not defined, please use file(string)");
  
  double tcur=coord_ph[0]; // in M units
  double t0 = posIni_[0];  // t0 M units

  double vel[4]; // 4-velocity of emitter
  for (int ii=0;ii<4;ii++){
    vel[ii]=coord_obj[ii+4];
  }

  // Defining jnus, anus
  double jnu[nbnu];
  double anu[nbnu];
  
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initializing to <0 value to create errors if not updated
    // [ exp(-anu*ds) will explose ]
    jnu[ii]=-1.;
    anu[ii]=-1.;
  }

  int X_params[3] = {nb_time_, nb_freq_, nb_angle_};
  double** X;
  X = new double*[3];
  X[0] = time_array_;
  X[1] = freq_array_;
  X[2] = angle_array_;

  std::string cond_limits[3] = {"Linear", "Linear", "Constant"}; // time, frequency, theta_mag
  double theta_mag = 10.; // avg
  if (nb_angle_ != 1.){ // Meaning the coefficients are not averaged over pitch angle
    cond_limits[2] = "None";
    double B4vect[4]={0.,0.,0.,0.};
    double Btor[4]={0.,0.,0.,0.};
    double Bpol[4]={0.,0.,0.,0.};
    computeB4vect(Btor, "Toroidal", coord_obj, coord_ph);
    computeB4vect(Bpol, "Poloidal", coord_obj, coord_ph);    
    for (int i=0;i<4;i++){
      B4vect[i] = Btor[i]/(beta_+1) + Bpol[i]*beta_/(beta_+1);
    }
    double norm=sqrt(gg_->ScalarProd(&coord_ph[0], B4vect, B4vect));
    if (fabs(norm-1.)>GYOTO_DEFAULT_ABSTOL) GYOTO_ERROR("Bad mf normalization");
    //gg_->multiplyFourVect(B4vect,1./norm);
    //computeB4vect(B4vect, magneticConfig_, coord_obj, coord_ph);
    theta_mag = get_theta_mag(B4vect, coord_ph, vel);
  }

  // KAPPA SYNCHRO
  double jnu_synch_kappa[nbnu], anu_synch_kappa[nbnu];
  double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*15.
    /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq
  cout << "nu0 : " << nu0 << endl;
  double hypergeom = Gyoto::hypergeom(4., 50.);
  spectrumKappaSynch_->kappaindex(4.);
  spectrumKappaSynch_->angle_averaged(1);
  spectrumKappaSynch_->angle_B_pem(theta_mag); // avg so we don't care
  spectrumKappaSynch_->cyclotron_freq(nu0);
  spectrumKappaSynch_->thetae(50.);
  spectrumKappaSynch_->hypergeometric(hypergeom);

  

  // COMPUTE VALUES IN FUNCTION OF PHASE
  if (tcur<=t0){
    for (size_t ii=0; ii<nbnu; ++ii){
      jnu[ii]=0;
      anu[ii]=0;
    }
  }
  else{
    double tt=(tcur-t0); // in M unit
    for (size_t ii=0; ii<nbnu; ++ii){
      double Xq[3] = {tt, nu_ems[ii], theta_mag};
      jnu[ii]=interpolate(3, emis_polar_array_[0], Xq, X, X_params, cond_limits);
      anu[ii]=interpolate(3, abs_polar_array_[0],  Xq, X, X_params, cond_limits);
      if ((tt-t_inj_)<1.){
        double ne = min(5.e6,5.e6*tt/t_inj_);
        spectrumKappaSynch_->numberdensityCGS(ne);
        spectrumKappaSynch_->radiativeQ(jnu_synch_kappa,anu_synch_kappa,nu_ems,nbnu);
        cout << "tt, nu_em, theta_mag : " << tt << ", " << nu_ems[ii] << ", " << theta_mag << endl;
        cout << "Compare jnu : " << jnu[ii] << ", " << jnu_synch_kappa[ii] << endl;
        cout << "Compare anu : " << anu[ii] << ", " << anu_synch_kappa[ii] << endl;
      }
    }
  }


  // RETURNING TOTAL INTENSITY AND TRANSMISSION
  for (size_t ii=0; ii<nbnu; ++ii){
    double jnu_tot = jnu[ii], anu_tot = anu[ii];
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

void Plasmoid::radiativeQ(double *Inu, double *Qnu, double *Unu, double *Vnu,
           Eigen::Matrix4d *Onu, double const *nuem , size_t nbnu,
           double dsem, state_t const &coord_ph, double const *coord_obj) const {
  cout << "Hello !" << endl;
  if (coord_ph.size()!=16){
    // Onu is the transmission matrix which contains in particular the non-polarised transmission
    // So we need to update it.
    Eigen::Matrix4d identity;
    identity << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;

    double* Taunu;
    Taunu = new double[nbnu];
    for (size_t ii=0; ii<nbnu; ++ii) {
      radiativeQ(Inu, Taunu, nuem, nbnu, dsem, coord_ph, coord_obj);
      Onu[ii] = Taunu[ii]*identity;
    }
  } else {
    if (filename_=="None")
      GYOTO_ERROR("In Plamsoid RadiativeQ : filename_ not defined, please use file(string)");

    cout << "Hello !" << endl;
    double tcur=coord_ph[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C; // in sec
    double t0 = posIni_[0]*GYOTO_G_OVER_C_SQUARE*gg_->mass()/GYOTO_C;  // t0 in sec

    double vel[4]; // 4-velocity of emitter
    for (int ii=0;ii<4;ii++){
      vel[ii]=coord_obj[ii+4];
    }

    double B4vect[4]={0.,0.,0.,0.};
    double Btor[4]={0.,0.,0.,0.};
    double Bpol[4]={0.,0.,0.,0.};
    computeB4vect(Btor, "Toroidal", coord_obj, coord_ph);
    computeB4vect(Bpol, "Poloidal", coord_obj, coord_ph);    
    for (int i=0;i<4;i++){
      B4vect[i] = Btor[i]/pow(beta_*beta_+1.,0.5) + Bpol[i]*beta_/pow(beta_*beta_+1.,0.5);
    }
    double norm=sqrt(gg_->ScalarProd(&coord_ph[0], B4vect, B4vect));
    if (fabs(norm-1.)>GYOTO_DEFAULT_ABSTOL) GYOTO_ERROR("Bad mf normalization");
    //gg_->multiplyFourVect(B4vect,1./norm);
    //computeB4vect(B4vect, magneticConfig_, coord_obj, coord_ph);

    double Chi=getChi(B4vect, coord_ph, vel); // this is EVPA

    double theta_mag = get_theta_mag(B4vect, coord_ph, vel);

    if (theta_mag<0. or theta_mag>M_PI) throwError("Blob: bad B angle");

    // Setup interpolation parameters
    int X_params[3] = {nb_time_, nb_freq_, nb_angle_};
    double** X;
    X = new double*[3];
    X[0] = time_array_;
    X[1] = freq_array_;
    X[2] = angle_array_;
    std::string cond_limits[3] = {"Linear", "Linear", "None"}; // time, frequency, theta_mag

    // Defining emission, absoprtion and rotation coefficients for the transmission matrix
    double jInu[nbnu], jQnu[nbnu], jUnu[nbnu], jVnu[nbnu];
    double aInu[nbnu], aQnu[nbnu], aUnu[nbnu], aVnu[nbnu];
    double rotQnu[nbnu], rotUnu[nbnu], rotVnu[nbnu];
    
    cout << "Hello !" << endl;
    for (size_t ii=0; ii<nbnu; ++ii){
      // Initialze them to -1 to create error if not updated
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

    if (tcur<=t0){
      for (size_t ii=0; ii<nbnu; ++ii){
        jInu[ii]=0.;
        jQnu[ii]=0.;
        jUnu[ii]=0.;
        jVnu[ii]=0.;
        aInu[ii]=0.;
        aQnu[ii]=0.;
        aUnu[ii]=0.;
        aVnu[ii]=0.;
        rotQnu[ii]=0.;
        rotUnu[ii]=0.;
        rotVnu[ii]=0.;
      }
    }
    else{
      double tt=(tcur-t0); // in sec
      for (size_t ii=0; ii<nbnu; ++ii){
        double Xq[3] = {tt, nuem[ii], theta_mag};
        jInu[ii]=interpolate(3, emis_polar_array_[0],  Xq, X, X_params, cond_limits);
        jQnu[ii]=interpolate(3, emis_polar_array_[1],  Xq, X, X_params, cond_limits);
        jUnu[ii]=interpolate(3, emis_polar_array_[2],  Xq, X, X_params, cond_limits);
        jVnu[ii]=interpolate(3, emis_polar_array_[3],  Xq, X, X_params, cond_limits);
        aInu[ii]=interpolate(3, abs_polar_array_[0],   Xq, X, X_params, cond_limits);
        aQnu[ii]=interpolate(3, abs_polar_array_[1],   Xq, X, X_params, cond_limits);
        aUnu[ii]=interpolate(3, abs_polar_array_[2],   Xq, X, X_params, cond_limits);
        aVnu[ii]=interpolate(3, abs_polar_array_[3],   Xq, X, X_params, cond_limits);
        rotQnu[ii]=interpolate(3, rot_polar_array_[0], Xq, X, X_params, cond_limits);
        rotUnu[ii]=interpolate(3, rot_polar_array_[1], Xq, X, X_params, cond_limits);
        rotVnu[ii]=interpolate(3, rot_polar_array_[2], Xq, X, X_params, cond_limits);
      }
    }

    for (size_t ii=0; ii<nbnu; ++ii) {
      Eigen::Vector4d Jstokes=rotateJs(jInu[ii], jQnu[ii], jUnu[ii], jVnu[ii], Chi);
      Eigen::Matrix4d Omat = Omatrix(aInu[ii], aQnu[ii], aUnu[ii], aVnu[ii], rotQnu[ii], rotUnu[ii], rotVnu[ii], Chi, dsem);
      Eigen::Matrix4d Pmat = Pmatrix(aInu[ii], aQnu[ii], aUnu[ii], aVnu[ii], rotQnu[ii], rotUnu[ii], rotVnu[ii], sin(2.*Chi), cos(2.*Chi), dsem);

      // Computing the increment of the Stokes parameters. Equivalent to dInu=exp(-anu*dsem)*jnu*dsem in the non-polarised case.
      Eigen::Vector4d Stokes=Pmat*Jstokes;

      if (Stokes(0) <=0.){
        Inu[ii] = Qnu[ii] =  Unu[ii] = Vnu[ii] = 0.;
        Onu[ii] = Eigen::Matrix4d::Identity();
      } else {
        Inu[ii] = Stokes(0);
        Qnu[ii] = Stokes(1);
        Unu[ii] = Stokes(2);
        Vnu[ii] = Stokes(3);
        Onu[ii] = Omat;
      }
      
      if (Inu[ii]<0.)
        GYOTO_ERROR("In Plasmoid::radiativeQ(): Inu<0");
      if (isnan(Inu[ii]) or isnan(Qnu[ii]) or isnan(Unu[ii]) or isnan(Vnu[ii]) or isnan(Onu[ii](0,0)))
        GYOTO_ERROR("In Plasmoid::radiativeQ(): Snu or Taunu is nan");
      if (isinf(Inu[ii]) or isinf(Qnu[ii]) or isinf(Unu[ii]) or isinf(Vnu[ii]) or isinf(Onu[ii](0,0)))
        GYOTO_ERROR("In Plasmoid::radiativeQ(): Snu or Taunu is infinite");
    }
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
    //cout << "t, r, theta, phi = " << tt << ", " << r << ", " << theta << ", " << phi << endl;

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
  
  if (flag_=="Helical") // Nathanail+20
  {
  	vel[0] = 1.;
  	vel[1] = fourveldt_[1];
  	vel[2] = 0.;
  	vel[3] = fourveldt_[3];
  	gg_->normalizeFourVel(pos, vel);
  }
  else if (flag_=="Helical-Ejection"){ // El Mellah+23
    vel[0] = 1.;
    vel[1] = max(0.8,fourveldt_[1]*pow(pos[1]/posIni_[1],0.5));
    vel[2] = 0.;
    vel[3] = fourveldt_[3];
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

	double radiusMin = 0.4;
  double t0 = posIni_[0];  // t0 in M unit

  size_t sz = ph -> parallelTransport()?16:8;
  state_t p1(sz);
  ph->getCoord(index, p1);
  double tcur = p1[0]; //tcur in M unit

  
  if (varyRadius_== "Varying")
  {
    if (tcur<=t0) radius(radiusMin);
    else radius(radiusMin+0.001*(tcur-t0)); // Nathanail+20
  }
	else if (varyRadius_== "Constant") radius(radiusMax_);
  else{
    GYOTO_ERROR("In Plasmoid::Impact operation on radius not recognized. Use Radius('Constant' or 'Varying')");
  }

	return Standard::Impact(ph, index, data);
}


void Plasmoid::file(std::string const &f) {
  # ifdef GYOTO_USE_CFITSIO
    fitsRead(f);
  # else
    GYOTO_ERROR("This Gyoto has no FITS i/o");
  # endif
}

#ifdef GYOTO_USE_CFITSIO
void Plasmoid::fitsRead(string filename) {
  // Remove first char if it is "!"
  if (filename.substr(0,1)=="!")
    filename.erase(0,1);

  GYOTO_MSG << "Plasmoid reading FITS file: " << filename << endl;

  filename_ = filename;
  fitsfile* fptr      = NULL;
  
  fptr = FitsRW::fitsOpen(filename);

  nb_time_   = FitsRW::fitsReadKey(fptr, "NB_TIME");
  nb_freq_   = FitsRW::fitsReadKey(fptr, "NB_FREQ");
  nb_angle_  = FitsRW::fitsReadKey(fptr, "NB_ANGL");
  t_inj_     = FitsRW::fitsReadKey(fptr, "T_INJ");

  //time_array_ = new double[nb_time_];
  //freq_array_ = new double[nb_freq_];
  //angle_array_= new double[nb_angle_];

  emis_polar_array_[0] = FitsRW::fitsReadHDUData(fptr, "J_I", nullptr, 0);
  emis_polar_array_[1] = FitsRW::fitsReadHDUData(fptr, "J_Q", nullptr, 0);
  emis_polar_array_[2] = FitsRW::fitsReadHDUData(fptr, "J_U", nullptr, 0);
  emis_polar_array_[3] = FitsRW::fitsReadHDUData(fptr, "J_V", nullptr, 0);
  abs_polar_array_[0]  = FitsRW::fitsReadHDUData(fptr, "ALPHA_I", nullptr, 0);
  abs_polar_array_[1]  = FitsRW::fitsReadHDUData(fptr, "ALPHA_Q", nullptr, 0);
  abs_polar_array_[2]  = FitsRW::fitsReadHDUData(fptr, "ALPHA_U", nullptr, 0);
  abs_polar_array_[3]  = FitsRW::fitsReadHDUData(fptr, "ALPHA_V", nullptr, 0);
  rot_polar_array_[0]  = FitsRW::fitsReadHDUData(fptr, "R_Q", nullptr, 0);
  rot_polar_array_[1]  = FitsRW::fitsReadHDUData(fptr, "R_U", nullptr, 0);
  rot_polar_array_[2]  = FitsRW::fitsReadHDUData(fptr, "R_V", nullptr, 0);
  time_array_          = FitsRW::fitsReadHDUData(fptr, "TIME", nullptr, 0);
  freq_array_          = FitsRW::fitsReadHDUData(fptr, "FREQUENCY", nullptr, 0);
  angle_array_         = FitsRW::fitsReadHDUData(fptr, "ANGLE", nullptr, 0);

  FitsRW::fitsClose(fptr);
  return;

}
#endif

void Plasmoid::beta(double beta){
  beta_=beta;
}

double Plasmoid::beta() const{
  return beta_;
}
