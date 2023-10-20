/*
    Copyright 2011 Thibaut Paumard, Frederic Vincent

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
#include "GyotoWorldline.h"
#include "GyotoStar.h"
#include "GyotoProperty.h"
#include "GyotoPhoton.h"
#include "GyotoPowerLawSpectrum.h"
#include "GyotoBlackBodySpectrum.h"
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
GYOTO_PROPERTY_START(Gyoto::Astrobj::Star,
 "UniformSphere following a time-like Gyoto::Worldline.")
// Star only need to implement the Worldline interface on top of the 
// UniformSphere interface, which is trivially tone with this macro:
GYOTO_WORLDLINE_PROPERTY_END(Star, UniformSphere::properties)

#define POLAR_TEST 0

// XML I/O
// We also need to parse and write Position+Velocity in addition to
// InitCoord, which is done by overriding setParameter(), setParameters()
// and fillProperty()
int Star::setParameter(std::string name,
			    std::string content,
			    std::string unit) {
  double coord[8];
  if (name=="InitialCoordinate") {
    name="InitCoord";
    return UniformSphere::setParameter(name, content, unit);
  } else if (name=="Position") {
    if (FactoryMessenger::parseArray(content, coord, 4) != 4)
      GYOTO_ERROR("Worldline \"Position\" requires exactly 4 tokens");
    if (init_vel_) {
      setInitCoord(coord, init_vel_);
      delete[] init_vel_; init_vel_=NULL;
    } else setPosition(coord);
    wait_pos_ = 0;
  } else if (name=="Velocity") {
    if (FactoryMessenger::parseArray(content, coord, 3) != 3)
      GYOTO_ERROR("Worldline \"Velocity\" requires exactly 3 tokens");
    if (wait_pos_) {
      if (init_vel_) delete [] init_vel_;
      init_vel_ = new double[3];
      memcpy(init_vel_, coord, 3*sizeof(double));
    } else setVelocity(coord);
  }
  else return UniformSphere::setParameter(name, content, unit);
  return 0;
}

#ifdef GYOTO_USE_XERCES
void Star::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const {
  if (p.name == "InitCoord") {
    if (imin_ <= imax_) {
      state_t coord;
      getInitialCoord(coord);
      // For massive particule, express initial condition with 3-velocity
      double vel[3] = {coord[5]/coord[4], coord[6]/coord[4], coord[7]/coord[4]};
      fmp -> setParameter ("Position", &coord[0], 4);
      fmp -> setParameter ("Velocity", vel, 3);
    }
    return;
  }
  UniformSphere::fillProperty(fmp, p);
}

void Star::setParameters(FactoryMessenger* fmp) {
  wait_pos_ = 1;
  metric(fmp->metric());
  UniformSphere::setParameters(fmp);
  wait_pos_ = 0;
  if (init_vel_) {
    delete[] init_vel_; init_vel_=NULL;
    GYOTO_ERROR("Worldline::setParameters(): "
	       "Velocity was found but not Position");
  }
}
#endif
///

Star::Star() :
  UniformSphere("Star"),
  Worldline(),
  spectrumThermalSynch_(NULL)
{
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
}

Star::Star(SmartPointer<Metric::Generic> met, double rad,
	   double const pos[4],
	   double const v[3]) :
  UniformSphere("Star"),
  Worldline(),
  spectrumThermalSynch_(NULL)
{
  if (debug()) {
    cerr << "DEBUG: Star Construction " << endl
	 << "       POS=[" << pos[0];
    for (int i=1; i<4; ++i) cerr << ", " << pos[i];
    cerr << "]\n       VEL=[" << v[0] ;
    for (int i=1; i<3; ++i) cerr << ", " << v[i];
    cerr << "]\n       RADIUS=" << rad << endl;

  }

  metric(met);
  setInitCoord(pos, v);
  radius(rad);
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
}

Star::Star(const Star& orig) :
  UniformSphere(orig), Worldline(orig), spectrumThermalSynch_(NULL)
{
  GYOTO_DEBUG << endl;
  // we have two distinct clones of the metric, not good...
  Worldline::metric(UniformSphere::metric());
  if (orig.spectrumThermalSynch_()) spectrumThermalSynch_=orig.spectrumThermalSynch_->clone();
}

Star* Star::clone() const { return new Star(*this); }

Star::~Star() {
  if (debug()) cerr << "DEBUG: Star::~Star()\n";
}

string Star::className() const { return  string("Star"); }
string Star::className_l() const { return  string("star"); }

SmartPointer<Metric::Generic> Star::metric() const { return gg_; }
void Star::metric(SmartPointer<Metric::Generic> gg) {
  UniformSphere::metric(gg);
  Worldline::metric(gg);
}

void Star::setInitialCondition(double const coord[8]) {
  if (!metric_) GYOTO_ERROR("Please set metric before calling Star::setInitialCondition(double*)");
  Worldline::setInitialCondition(metric_, coord, 0);
}

double Star::getMass() const {return 1. ;}

void Star::getVelocity(double const pos[4], double vel[4]) {
  getCoord(pos, 1, NULL, NULL, NULL, vel, vel+1, vel+2, vel+3);
}

void Star::getCartesian(double const * const t,
			size_t const n,
			double* const x, double*const y, double*const z,
			double*const xp, double*const yp, double*const zp) {
  Worldline::getCartesian(t, n, x, y, z, xp, yp, zp);
}


double Star::rMax() {
  if (rmax_==DBL_MAX && i0_>=imin_ && i0_<=imax_) {
    size_t i;
    rmax_=x1_[i0_];
    int ck=gg_->coordKind();
    for (i=imin_;i<=imax_;++i) {
      if (x1_[i]>rmax_) rmax_=x1_[i];
      if (ck==GYOTO_COORDKIND_CARTESIAN) {
	if (x2_[i]>rmax_) rmax_=x2_[i];
	if (x3_[i]>rmax_) rmax_=x3_[i];
      }
    }
    rmax_ *= 3.;
  }
  return rmax_;
}

void Star::radiativeQ(double *Inu, double *Qnu, double *Unu, double *Vnu,
       Eigen::Matrix4d *Onu, double const *nuem , size_t nbnu, double dsem,
       state_t const &cph, double const *co) const {

  if (POLAR_TEST==1)
  {
    // The following part was only for test purpose with ad-hoc choice of magnetic field (ipole formalism)
    // and emission process

    double rr, rcyl, theta, phi, xx, yy, zz=0.;
    switch (gg_->coordKind()) {
    case GYOTO_COORDKIND_SPHERICAL:
      rr = cph[1];
      theta = cph[2];
      phi = cph[3];
      
      rcyl = rr*sin(theta);
      xx = rcyl*cos(phi);
      yy = rcyl*sin(phi);
      zz   = rr*cos(theta);
      break;
    case GYOTO_COORDKIND_CARTESIAN:
      rcyl = pow(cph[1]*cph[1]+cph[2]*cph[2], 0.5);
      rr = sqrt(cph[1]*cph[1]+cph[2]*cph[2]
          +cph[3]*cph[3]);
      theta   = acos(cph[3]/rr);
      xx = cph[1];
      yy = cph[2];
      zz   = cph[3];
      break;
    default:
      GYOTO_ERROR("In Star::radiativeQ(): Unknown coordinate system kind");
    }

    double vel[4]; // 4-velocity of emitter
    gg_->circularVelocity(co, vel);
    
    // Check the size of the photon coordinate : 8 means no polarisation, 16 means polarisation.
    if (cph.size()==16) { // POLARISED CASE
      double B4vect[4]={0.,0.,0.,0.};
      double B_1=0.,B_2=0.,B_3=0;
      double spin = 0.;
      double gtt = gg_->gmunu(&cph[0],0,0),
             grr = gg_->gmunu(&cph[0],1,1),
             gthth = gg_->gmunu(&cph[0],2,2),
             gpp = gg_->gmunu(&cph[0],3,3);
      double dx1=0.025,
             dx2=0.025;

      double g_det = sqrt(M_PI*M_PI*pow(rr,6)*pow(sin(theta),2));

      double F11 = exp(log(rr)-dx1)*sin(theta-dx2*M_PI),
             F12 = exp(log(rr)-dx1)*sin(theta+dx2*M_PI),
             F21 = exp(log(rr)+dx1)*sin(theta-dx2*M_PI),
             F22 = exp(log(rr)+dx1)*sin(theta+dx2*M_PI);
      B_1 = -(F11-F12+F21-F22)/(2.*dx2*g_det);
      B_2 =  (F11+F12-F21-F22)/(2.*dx1*g_det);
      B_3 = 0.;
      

      // compute contravariant velocity in KS' from BL
      double dtKS_drBL   = 2. * rr / (rr*rr - 2.*rr + spin*spin);
      double dphiKS_drBL = spin / (rr*rr - 2.*rr + spin*spin);
      double Ucon_KSm[4]={0.,0.,0.,0.};
      Ucon_KSm[0]=vel[0]+vel[1]*dtKS_drBL;
      Ucon_KSm[1]=vel[1]/rr;
      Ucon_KSm[2]=vel[2]/M_PI;
      Ucon_KSm[3]=vel[3]+vel[1]*dphiKS_drBL;

      // Compute KS' metric
      double gcov_ksm[4][4];
      double sin2=sin(theta)*sin(theta), rho2=rr*rr+spin*spin*cos(theta)*cos(theta);
      double gcov_ks[4][4];
      for(int mu=0;mu<4;mu++)
        for(int nu=0;nu<4;nu++)
          gcov_ks[mu][nu]=0.;

      gcov_ks[0][0] = -1. + 2. * rr / rho2 ;
      gcov_ks[0][1] = 2. * rr / rho2 ;
      gcov_ks[0][3] = -2. * spin * rr * sin(theta)*sin(theta) / rho2;
      gcov_ks[1][0] = gcov_ks[0][1];
      gcov_ks[1][1] = 1. + 2. * rr / rho2 ;
      gcov_ks[1][3] = -spin * sin(theta)*sin(theta) * (1. + 2. * rr / rho2);
      gcov_ks[2][2] = rho2 ;
      gcov_ks[3][0] = gcov_ks[0][3];
      gcov_ks[3][1] = gcov_ks[1][3];
      gcov_ks[3][3] = sin(theta)*sin(theta) * (rho2 + spin * spin * sin(theta)*sin(theta) * (1. + 2. * rr / rho2));

      // convert from ks metric to a modified one using Jacobian
      double dxdX[4][4];
      double hslope=0.;
      for(int mu=0;mu<4;mu++)
        for(int nu=0;nu<4;nu++)
          dxdX[mu][nu]=0.;

      dxdX[0][0] = 1.;
      dxdX[1][1] = rr;
      dxdX[2][2] = M_PI  + hslope*2*M_PI*cos(2*theta); 
      dxdX[3][3] = 1.;

      for(int mu=0;mu<4;mu++){
        for(int nu=0;nu<4;nu++){
          gcov_ksm[mu][nu] = 0;
          for (int lam = 0; lam < 4; ++lam) {
            for (int kap = 0; kap < 4; ++kap) {
              gcov_ksm[mu][nu] += gcov_ks[lam][kap] * dxdX[lam][mu] * dxdX[kap][nu];
            }
          }
        }
      }

      // Compute covariante velocity in KS'
      double Ucov_KSm[4]={0.,0.,0.,0.};
      for(int mu=0;mu<4;mu++){
        for(int nu=0;nu<4;nu++){
          Ucov_KSm[mu] += gcov_ksm[mu][nu]*Ucon_KSm[nu];
        }
      }

      // Copute Magnetic field in KS'
      //cout << "r sth, velBL, ucov KS'= " << co[1] << " " << sin(co[2]) << " " << vel[0] << " " << vel[3] << " " << Ucov_KSm[1] << " " << Ucov_KSm[2] << " " << Ucov_KSm[3] << endl;
      //throwError("test disk");
      double B0=B_1*Ucov_KSm[1]+B_2*Ucov_KSm[2]+B_3*Ucov_KSm[3],
        B1=(B_1+B0*Ucon_KSm[1])/Ucon_KSm[0],
        B2=(B_2+B0*Ucon_KSm[2])/Ucon_KSm[0],
        B3=(B_3+B0*Ucon_KSm[3])/Ucon_KSm[0];

      // Conversion Magnetic field from KS' -> BL
      double Delta = pow(rr,2)-2.*rr+pow(spin,2.);
      B4vect[0]=B0-B1*2.*pow(rr,2)/Delta;
      B4vect[1]=B1*rr;
      B4vect[2]=B2*M_PI;
      B4vect[3]=B3-B1*spin*rr/Delta;

      double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*10.
      /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq

      Eigen::Matrix4d Omat;
      Omat << 1, 0, 0, 0,
          0, 1, 0, 0,
          0, 0, 1, 0,
          0, 0, 0, 1;

      double Theta = 50.;
      double temperature=GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS*Theta/GYOTO_BOLTZMANN_CGS;

      // Computing the angle between the parallel transported observer polarisation basis and local emission basis.
      double Chi=getChi(B4vect, cph, vel);

      // Computing the angle theta_mag between the magnetic field vector and photon tgt vector in the rest frame of the emitter
      gg_->projectFourVect(&cph[0],B4vect,vel); //Projection of the 4-vector B to 4-velocity to be in the rest frame of the emitter
      double photon_emframe[4]; // photon tgt vector projected in comoving frame
      for (int ii=0;ii<4;ii++){
        photon_emframe[ii]=cph[ii+4];
      }
      gg_->projectFourVect(&cph[0],photon_emframe,vel);
      double bnorm = gg_->norm(&cph[0],B4vect);
      double lnorm = gg_->norm(&cph[0],photon_emframe);
      double lscalb = gg_->ScalarProd(&cph[0],photon_emframe,B4vect);
      double theta_mag = acos(lscalb/(lnorm*bnorm));

      // Defining emission, absoprtion and rotation coefficients for the transmission matrix
      double jInu[nbnu], jQnu[nbnu], jUnu[nbnu], jVnu[nbnu];
      double aInu[nbnu], aQnu[nbnu], aUnu[nbnu], aVnu[nbnu];
      double rotQnu[nbnu], rotUnu[nbnu], rotVnu[nbnu];
      
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

      double besselK2 = bessk(2, 1./Theta);
      // THERMAL SYNCHROTRON
      spectrumThermalSynch_->temperature(temperature);
      spectrumThermalSynch_->numberdensityCGS(1.e6);
      spectrumThermalSynch_->angle_averaged(0); //  no angle avg of course
      spectrumThermalSynch_->angle_B_pem(theta_mag);   
      spectrumThermalSynch_->cyclotron_freq(nu0);
      spectrumThermalSynch_->besselK2(besselK2);
      spectrumThermalSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu,
                                        aInu, aQnu, aUnu, aVnu,
                                        rotQnu, rotUnu, rotVnu, nuem, nbnu);

      for (size_t ii=0; ii<nbnu; ++ii) {
      //cout << "In Star: jInu, jQnu, jUnu, jVnu: " << jInu[ii] << ", " << jQnu[ii] << ", " << jUnu[ii] << ", " << jVnu[ii] << endl;
      //cout << "In Star: aInu, aQnu, aUnu, aVnu: " << aInu[ii] << ", " << aQnu[ii] << ", " << aUnu[ii] << ", " << aVnu[ii] << endl;
      //cout << "In Star: rQnu, rUnu, rVnu: " << rotQnu[ii] << ", " << rotUnu[ii] << ", " << rotVnu[ii] << endl;
      Eigen::Vector4d Jstokes=rotateJs(jInu[ii], jQnu[ii], jUnu[ii], jVnu[ii], Chi)*dsem*gg_->unitLength();
      //cout << Jstokes << endl;
      Omat = Omatrix(aInu[ii], aQnu[ii], aUnu[ii], aVnu[ii], rotQnu[ii], rotUnu[ii], rotVnu[ii], Chi, dsem);
      //cout << Omat << endl;
      // Computing the increment of the Stokes parameters. Equivalent to dInu=exp(-anu*dsem)*jnu*dsem in the non-polarised case.
      Eigen::Vector4d Stokes=Omat*Jstokes;
      //cout << Stokes << endl;
      Inu[ii] = Stokes(0);
      Qnu[ii] = Stokes(1);
      Unu[ii] = Stokes(2);
      Vnu[ii] = Stokes(3);
      Onu[ii] = Omat;

      //cout << "In Star: r,th,ph, Inu, Qnu, Unu, Vnu, dsem, LP: " << rr << " " << theta << " " << phi << " " << Inu[ii] << ", " << Qnu[ii] << ", " << Unu[ii] << ", " << Vnu[ii] << ", " << dsem << ", " << pow(Qnu[ii]*Qnu[ii]+Unu[ii]*Unu[ii],0.5)/Inu[ii] << endl;

      if (Inu[ii]<0.)
        GYOTO_ERROR("In Star::radiativeQ(): Inu<0");
      if (Inu[ii]!=Inu[ii] or Onu[ii](0,0)!=Onu[ii](0,0))
        GYOTO_ERROR("In Star::radiativeQ(): Inu or Taunu is nan");
      if (Inu[ii]==Inu[ii]+1. or Onu[ii](0,0)==Onu[ii](0,0)+1.)
        GYOTO_ERROR("In Star::radiativeQ(): Inu or Taunu is infinite");
      }

    } else { // NON POLARISED CASE

      // Onu is the transmission matrix which contains in particular the non-polarised transmission
      // So we need to update it.
      Eigen::Matrix4d identity;
      identity << 1, 0, 0, 0,
                  0, 1, 0, 0,
                  0, 0, 1, 0,
                  0, 0, 0, 1;

      for (size_t ii=0; ii<nbnu; ++ii) {
        // If RadiativeQ(non-polarised) is implemented, it should be called instead of emission and transmission.
        Inu[ii] = emission(nuem[ii], dsem, cph, co);
        double Tau = transmission(nuem[ii], dsem, cph, co);
        //cout << "Tau= " << Tau << endl;
        Onu[ii] = Tau*identity;
      }
    }
  }
  else{
    // NON TESTING POLARISED RADIATIVEQ TO BE DONE
    Eigen::Matrix4d identity;
    identity << 1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0,
                0, 0, 0, 1;

    for (size_t ii=0; ii<nbnu; ++ii) {
      Inu[ii] = emission(nuem[ii], dsem, cph, co);
      double Tau = transmission(nuem[ii], dsem, cph, co);
      Qnu[ii] = 0.;
      Unu[ii] = 0.;
      Vnu[ii] = 0.;
      Onu[ii] = Tau*identity;
    }
  }
}