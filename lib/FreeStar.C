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
#include "GyotoFreeStar.h"
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
GYOTO_PROPERTY_START(Gyoto::Astrobj::FreeStar,
 "UniformSphere with a user defined orbit (not time-like).")
GYOTO_PROPERTY_VECTOR_DOUBLE(FreeStar, InitPosition, initPosition,
              "(t,r,theta,phi) initial position of freeStar")
GYOTO_PROPERTY_VECTOR_DOUBLE(FreeStar, InitVelocity, initVelocity,
              "(dr/dt,dtheta/dt,dphi/dt) initial 3-velocity of freeStar")
GYOTO_PROPERTY_END(FreeStar, UniformSphere::properties)

FreeStar::FreeStar() : 
  UniformSphere("FreeStar"),
  posSet_(false),
  posIni_(NULL),
  fourveldt_(NULL)
{
  kind_="FreeStar";
# ifdef GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done." << endl;
# endif

  posIni_= new double[4];
  fourveldt_= new double[4];
}

FreeStar::FreeStar(const FreeStar& orig) :
  UniformSphere(orig),
  posSet_(orig.posSet_),
  posIni_(NULL),
  fourveldt_(NULL)
{

  if(orig.posIni_){
      posIni_= new double[4];
      memcpy(posIni_,orig.posIni_, 4*sizeof(double));
  }

  if(orig.fourveldt_){
      fourveldt_= new double[4];
      memcpy(fourveldt_,orig.fourveldt_, 4*sizeof(double));
  }
}


FreeStar* FreeStar::clone() const { return new FreeStar(*this); }

FreeStar::~FreeStar() {
  if (debug()) cerr << "DEBUG: FreeStar::~FreeStar()\n";
}

string FreeStar::className() const { return  string("FreeStar"); }
string FreeStar::className_l() const { return  string("freeStar"); }


void FreeStar::initPosition(std::vector<double> const &v) {
  posIni_[0] = v[0];
  posIni_[1] = v[1];
  posIni_[2] = v[2];
  posIni_[3] = v[3];
  posSet_=true;
}

std::vector<double> FreeStar::initPosition() const {
  std::vector<double> v (4, 0.);
  v[0] = posIni_[0];
  v[1] = posIni_[1];
  v[2] = posIni_[2];
  v[3] = posIni_[3];
  return v;
}

void FreeStar::initVelocity(std::vector<double> const &v) {
  if (!posSet_)
    GYOTO_ERROR("In FreeStar::initVelocity initial Position not defined");
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
    GYOTO_ERROR("In FreeStar::initVelocity Initial Velocity over C");

  gg_->normalizeFourVel(posIni_, fourveldt_);

}

std::vector<double> FreeStar::initVelocity() const {
  std::vector<double> v (3, 0.);
  v[0] = fourveldt_[1];
  v[1] = fourveldt_[2];
  v[2] = fourveldt_[3];
  return v;
}

void FreeStar::initCoord(std::vector<double> const &v) {
  posIni_[0] = v[0];
  posIni_[1] = v[1];
  posIni_[2] = v[2];
  posIni_[3] = v[3];
  fourveldt_[0] = v[4];
  fourveldt_[1] = v[5];
  fourveldt_[2] = v[6];
  fourveldt_[3] = v[7];
}

std::vector<double> FreeStar::initCoord() const {
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

void FreeStar::getCartesian(double const * const dates, size_t const n_dates,
          double * const x, double * const y, double * const z, 
          double * const xprime, double * const yprime, double * const zprime){
  // this yields the position of the center of the UnifSphere
  // at time t
  // fourveldt_ is the initial 3-velocity dxi/dt
  // vel is the 4-velocity dxnu/dtau

  if (n_dates!=1)
    GYOTO_ERROR("In FreeStar::getCartesian n_dates!=1");

  double tt=dates[0];
  
  double r, theta, phi; // spherical coordinates
  double vel[4];
  getVelocity(posIni_, vel);

  r = posIni_[1]+vel[1]/vel[0]*(tt-posIni_[0]);
  theta = posIni_[2]+vel[2]/vel[0]*(tt-posIni_[0]);
  phi = posIni_[3] + vel[3]/vel[0]*(tt-posIni_[0]);

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

void FreeStar::getVelocity(double const pos[4], double vel[4]){
  if (!gg_)
    GYOTO_ERROR("In FreeStar::getVelocity Metric not set");
    
  vel[0] = fourveldt_[0];
  vel[1] = fourveldt_[1];
  vel[2] = fourveldt_[2];
  vel[3] = fourveldt_[3];  
}