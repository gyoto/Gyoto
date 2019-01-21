/*
  Copyright 2015 Frederic Vincent, 2017 Thibaut Paumard

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

// Lorene headers
#include "metric.h"
#include "nbr_spx.h"
#include "utilitaires.h"
#include "graphique.h"

//Gyoto headers
#include "GyotoUtils.h"
#include "GyotoPhoton.h"
#include "GyotoNeutronStar.h"
#include "GyotoFactoryMessenger.h"

//Std headers
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <cstring>
#include <sstream>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;
using namespace Lorene;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(NeutronStar,
		     "Neutron star emitting at its surface.")
GYOTO_PROPERTY_END(NeutronStar, Standard::properties)

NeutronStar::NeutronStar() : Standard("NeutronStar"), gg_(NULL) {
  GYOTO_DEBUG << endl;
  Generic::gg_=gg_;
}

NeutronStar::NeutronStar(std::string kin) :
  Standard(kin), gg_(NULL)
{
  GYOTO_DEBUG << endl;
  Generic::gg_=gg_;
}

NeutronStar::NeutronStar(const NeutronStar& o) :
  Standard(o),
  gg_(NULL)
{
  GYOTO_DEBUG << endl;
  if (o.gg_()) gg_=o.gg_->clone();
  Generic::gg_=gg_;
}
NeutronStar * NeutronStar::clone() const { return new NeutronStar(*this); }

NeutronStar::~NeutronStar() {
  GYOTO_DEBUG << endl;
}

SmartPointer<Gyoto::Metric::Generic> NeutronStar::metric() const {
  GYOTO_DEBUG << endl;
  return gg_;
}

void NeutronStar::metric(SmartPointer<Gyoto::Metric::Generic> met) {
  GYOTO_DEBUG << endl;
  SmartPointer<Metric::NumericalMetricLorene> smptr =
    SmartPointer<Metric::NumericalMetricLorene>(met);
  if (met && !smptr) {
    // The above cast will yield a null pointer if met is not a
    // NumericalMetricLorene. It's an error if smptr is null but met
    // is not.
    GYOTO_ERROR("NeutronStar::metric(): metric should "
    "be a NumericalMetricLorene");
  }
  gg_ = smptr;
  Generic::metric(met);
}

double NeutronStar::operator()(double const coord[4]) {
  GYOTO_DEBUG << endl;
  if (gg_->coordKind() != GYOTO_COORDKIND_SPHERICAL){
    GYOTO_ERROR("In NeutronStar::operator(): so far only spherical coord");
  }

  double rcur = coord[1], thcur=coord[2], phcur=coord[3];
  Valeur* ns_surf = gg_->getNssurf_tab()[0];
  ns_surf->std_base_scal();
  double rstar = ns_surf->val_point(0,0.,thcur,phcur);

  //cout << "rcur rstar in NS= " << rcur << " " << rstar << endl;
  return rcur-rstar;
}

void NeutronStar::getVelocity(double const pos[4], double uu[4]){
  GYOTO_DEBUG << endl;
  double rr=pos[1], th=pos[2], phi=pos[3];
  double rsinth = rr*sin(th);
  if (rr==0.) GYOTO_ERROR("In NeutronStar.C::computeVelSurf r is 0!");
  if (rsinth==0.) GYOTO_ERROR("In NeutronStar.C::computeVelSurf on z axis!");
  double rm1 = 1./rr, rm2 = rm1*rm1, sm1 = 1./sin(th),
    sm2 = sm1*sm1, rsm1 = rm1*sm1;

  const Vector& v_i = *(gg_->getVsurf_tab()[0]); // [0] means at t=0 (stationary spacetime here!)
  double v_r = v_i(1).val_point(rr,th,phi),
    v_t = rr*v_i(2).val_point(rr,th,phi),
    v_p = rr*sin(th)*v_i(3).val_point(rr,th,phi);

  const Sym_tensor& g_up_ij = *(gg_->getGamcon_tab()[0]);
  double grr=g_up_ij(1,1).val_point(rr,th,phi), 
    gtt=rm2*g_up_ij(2,2).val_point(rr,th,phi),
    gpp=rm2*sm2*g_up_ij(3,3).val_point(rr,th,phi);
  double vr = v_r*grr, vt = v_t*gtt, vp = v_p*gpp; //contravariant 3-velocity
  
  //cout << "3v= " << vr << " " << vt << " " << vp << endl;
  
  Scalar* lorentz_scal = gg_->getLorentz_tab()[0];
  double lorentz = lorentz_scal->val_point(rr,th,phi);
  const Vector& shift = *(gg_->getShift_tab()[0]);
  double betar = shift(1).val_point(rr,th,phi),
    betat = rm1*shift(2).val_point(rr,th,phi),
    betap = rsm1*shift(3).val_point(rr,th,phi);
  //cout << "beta= " << betar << " " << betat << " " << betap << endl;
  Scalar* lapse_scal = gg_->getLapse_tab()[0];
  double lapse = lapse_scal->val_point(rr,th,phi);

  uu[0] = lorentz/lapse;
  uu[1] = lorentz*(vr - betar/lapse);
  uu[2] = lorentz*(vt - betat/lapse);
  uu[3] = lorentz*(vp - betap/lapse); //emitter 4-velocity
  //cout << "NS vel v3, u3= " << vp << " " << uu[3] << endl;
  /*cout << "4-vel= " ;
  for (int ii=0;ii<4;ii++) cout << uu[ii] << " ";
  cout << endl;*/

}

