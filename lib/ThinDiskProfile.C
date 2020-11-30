/*
    Copyright 2020 Frederic Vincent, Thibaut Paumard

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

#include "GyotoPhoton.h"
#include "GyotoThinDiskProfile.h"
#include "GyotoProperty.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"


#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <string>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

GYOTO_PROPERTY_START(ThinDiskProfile)
GYOTO_PROPERTY_END(ThinDiskProfile, ThinDisk::properties)


ThinDiskProfile::ThinDiskProfile() :
  ThinDisk("ThinDiskProfile")
{
  if (debug()) cerr << "DEBUG: ThinDiskProfile Construction" << endl;
}

ThinDiskProfile::ThinDiskProfile(const ThinDiskProfile& o) :
  ThinDisk(o)
{
  if (o.gg_()) gg_=o.gg_->clone();
  Generic::gg_=gg_;
}
ThinDiskProfile* ThinDiskProfile::clone() const
{ return new ThinDiskProfile(*this); }

bool ThinDiskProfile::isThreadSafe() const {
  return ThinDisk::isThreadSafe();
}

ThinDiskProfile::~ThinDiskProfile() {
  if (debug()) cerr << "DEBUG: ThinDiskProfile Destruction" << endl;
}

double ThinDiskProfile::emission(double nu, double,
			    state_t const &,
			    double const coord_obj[8]) const{
  double rr = coord_obj[1];
  // Gralla+20 model for M87
  double spin=0.94, a2=spin*spin;
  double rhor=1.+sqrt(1.-a2), rminus=1.-sqrt(1.-a2),
    z1 = 1. + pow(1.-a2,1./3.)*(pow(1.+spin,1./3.) + pow(1.-spin,1./3.)),
    z2 = sqrt(3.*a2 + z1*z1),
    risco = 3. + z2 - sqrt((3.-z1)*(3.+z1+2.*z2));

  // Choose profile here:
  double gamma=-3./2., mu=rminus, sigG=1./2.;
  //double gamma=-3., mu=risco-0.33, sigG=0.25;
  
  double tmp = gamma+asinh((rr-mu)/sigG);
  double emiss = exp(-0.5*tmp*tmp)/sqrt((rr-mu)*(rr-mu)+sigG*sigG);
  return 1e-6*emiss; // 1e-6 just for getting a reasonable flux
}

void ThinDiskProfile::getVelocity(double const pos[4], double vel[4])
{
  
  double risco = 0.;
  if (gg_->kind()!="Minkowski" && gg_->kind()!="Hayward")
    risco=gg_->getRms(); // prograde Kerr ISCO
  // let risco=0 if metric is Minko; then ISCO not defined
  // let also risco=0 for Hayward as we would need to
  // compute it numerically and give it in xml Metric field,
  // not implemented so far

  //cout << "in velo, r isco= " << pos[1] << " " << risco << endl;
  
  if (pos[1] > risco){
    // Keplerian velocity above ISCO
    gg_ -> circularVelocity(pos, vel, 1);
  }else{
    double gpp = gg_->gmunu(pos,3,3), gtt = gg_->gmunu(pos,0,0),
      gtp = gg_->gmunu(pos,0,3), grr = gg_->gmunu(pos,1,1);
    double utZAMO = sqrt(-gpp/(gtt*gpp-gtp*gtp)),
      uphiZAMO = -utZAMO*gtp/gpp;

    double V = 0.61; // velo norm as observed by ZAMO
    double Gamma = 1./sqrt(1.-V*V);
    double rr = pos[1];
    double rhor=1.6;
    double Vphi_over_V = (rr-rhor)/(risco-rhor);

    double Vphi = Vphi_over_V*V / sqrt(gpp),
      Vr = sqrt(1-Vphi_over_V*Vphi_over_V)*V / sqrt(grr);
    
    vel[0] = Gamma*utZAMO;
    vel[1] = -Gamma*Vr; // minus sign coz matter is going towards BH
    vel[2] = 0.;
    vel[3] = Gamma*(uphiZAMO + Vphi);

    //cout << "V2= " << gg_->gmunu(pos,1,1)*Vr*Vr + gg_->gmunu(pos,3,3)*Vphi*Vphi << endl;
    //cout << "u2 = " << gg_->ScalarProd(pos,vel,vel) << endl;
  }

  //cout << "4vel= " << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3]<< endl;
}

void ThinDiskProfile::processHitQuantities(Photon* ph,
					   state_t const &coord_ph_hit,
					   double const *coord_obj_hit,
					   double dt,
					   Properties* data) const {
#if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
#endif
  double dlambda = dt/coord_ph_hit[4]; //dlambda = dt/tdot
  double ggredm1 = -gg_->ScalarProd(&coord_ph_hit[0],coord_obj_hit+4,
				    &coord_ph_hit[4]);// / 1.; 
                                       //this is nu_em/nu_obs
  if (noredshift_) ggredm1=1.;
  double ggred = 1./ggredm1;           //this is nu_obs/nu_em
  double dsem = dlambda*ggredm1; // *1.
  double inc =0.;

  if (data){ // this check is necessary as process can be called
    // with data replaced by NULL (see ComplexAstrobj and
    // Photon nb_cross_eqplane)
    if (data->user4) {
      // *** CAREFUL!! ***
      /*
	I have to include here the "fudge factor" of Gralla+.
	Do not forget to remove it to consider some other disk profile.
      */
      double max_cross_eqplane = ph->maxCrossEqplane();
      if (max_cross_eqplane==DBL_MAX) cout << "WARNING: in ThinDiskProfile::process: max_cross_eqplane is DBL_MAX and probably should not be" << endl;
      int nb_cross_eqplane = ph->nb_cross_eqplane();
      double fudge_Gralla=1.;
      if (nb_cross_eqplane>0) fudge_Gralla=1.5;
      inc = fudge_Gralla
	* (emission(ggredm1, dsem, coord_ph_hit, coord_obj_hit))
	* (ph -> getTransmission(size_t(-1)))
	* ggred*ggred*ggred*ggred; // I/nu^4 invariant
      *data->user4 += inc;
#if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(*data->user4);
#endif
      
    }else GYOTO_ERROR("unimplemented data");
  }
}
