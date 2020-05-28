/*
    Copyright 2016, 2018-2020 Frederic Vincent, Thibaut Paumard

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
#include "GyotoEquatorialHotSpot.h"
#include "GyotoPhoton.h"
#include "GyotoPageThorneDisk.h"
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
#include <cstring>
#include <time.h> 

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

/// Properties
#include "GyotoProperty.h"
GYOTO_PROPERTY_START(EquatorialHotSpot, "Equatorial hot spot with beaming")
GYOTO_PROPERTY_DOUBLE(EquatorialHotSpot, SpotRadSize, spotRadSize)
GYOTO_PROPERTY_STRING(EquatorialHotSpot, BeamingKind, beaming,
		      "One of: IsotropicBeaming, NormalBeaming, RadialBeaming, "
		      "IsotropicConstant (emission is isotropic and constant"
		      "equals to 1)")
GYOTO_PROPERTY_DOUBLE(EquatorialHotSpot, BeamAngle, beamAngle)
GYOTO_WORLDLINE_PROPERTY_END(EquatorialHotSpot, ThinDisk::properties)

// accessors
void EquatorialHotSpot::spotRadSize(double t) {sizespot_=t;}
double EquatorialHotSpot::spotRadSize() const {return sizespot_;}

void EquatorialHotSpot::beaming(std::string const &b) {
  if (b=="IsotropicBeaming") beaming_=IsotropicBeaming;
  else if (b=="NormalBeaming") beaming_=NormalBeaming;
  else if (b=="RadialBeaming") beaming_=RadialBeaming;
  else if (b=="IsotropicConstant") beaming_=IsotropicConstant;
  else GYOTO_ERROR("Unknown beaming kind");
}
std::string EquatorialHotSpot::beaming() const {
  string b;
  switch (beaming_) {
  case IsotropicBeaming: b="IsotropicBeaming"; break;
  case NormalBeaming:    b="NormalBeaming";    break;
  case RadialBeaming:    b="RadialBeaming";    break;
  case IsotropicConstant: b="IsotropicConstant"; break;
  default: GYOTO_ERROR("Unknown beaming kind");
  }
  return b;
}

void EquatorialHotSpot::beamAngle(double t) {beamangle_=t;}
double EquatorialHotSpot::beamAngle() const {return beamangle_;}

// Needed for legacy XML files
int EquatorialHotSpot::setParameter(string name, string content, string unit) {
  double coord[8];
  char* tc = const_cast<char*>(content.c_str());
  if (name=="InitialCoordinate") {
    name="InitCoord";
    return ThinDisk::setParameter(name, content, unit);
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
  } else if (name=="NormalBeaming") {
    GYOTO_WARNING << "<" << name << "/> is deprecated, please use "
      "<BeamingKind> " << name << " </BeamingKind> instead";
    beaming(name);
  } else if (name=="NormalBeaming" || name=="RadialBeaming") {
    GYOTO_WARNING << "<" << name << "/> is deprecated, please use \n";
    GYOTO_WARNING << "<BeamingKind> " << name << " </BeamingKind>" << endl;
    GYOTO_WARNING << "<BeamAngle> " << content << "</BeamAngle>" << endl;
    GYOTO_WARNING <<" instead";
    beaming(name);
    setParameter("BeamAngle", content, unit);
  } else return ThinDisk::setParameter(name, content, unit);
  return 0;
}

// Needed for wait_pos_
#ifdef GYOTO_USE_XERCES
void EquatorialHotSpot::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const {
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
  ThinDisk::fillProperty(fmp, p);
}

void EquatorialHotSpot::setParameters(FactoryMessenger* fmp) {
  wait_pos_ = 1;
  ThinDisk::setParameters(fmp);
  wait_pos_ = 0;
  if (init_vel_) {
    delete[] init_vel_; init_vel_=NULL;
    GYOTO_ERROR("Worldline::setParameters(): "
	       "Velocity was found but not Position");
  }
}
#endif
///

Gyoto::Astrobj::EquatorialHotSpot::EquatorialHotSpot()
  : ThinDisk("EquatorialHotSpot"), Worldline(), 
    sizespot_(0.), beaming_(IsotropicBeaming), beamangle_(0.)
{
  GYOTO_DEBUG  << "Building EquatorialHotSpot";
}

Gyoto::Astrobj::EquatorialHotSpot::EquatorialHotSpot(const EquatorialHotSpot &o)
  : ThinDisk(o), Worldline(o),
    sizespot_(o.sizespot_), beaming_(o.beaming_), beamangle_(o.beamangle_)
{
  GYOTO_DEBUG  << "Copying EquatorialHotSpot";
}
EquatorialHotSpot * EquatorialHotSpot::clone() const { 
  return new EquatorialHotSpot(*this); }

Gyoto::Astrobj::EquatorialHotSpot::~EquatorialHotSpot()
{
  GYOTO_DEBUG  << "Destroying EquatorialHotSpot";
}

double EquatorialHotSpot::getMass() const {return 1. ;}

void EquatorialHotSpot::metric(SmartPointer<Metric::Generic> gg) {
  ThinDisk::metric(gg);
  Worldline::metric(gg);
}

void EquatorialHotSpot::setInitialCondition(double coord[8]) {
  if (!metric_) GYOTO_ERROR("Please set metric before calling "
			   "EquatorialHotSpot::setInitialCondition(double*)");
  Worldline::setInitialCondition(metric_, coord, 1);
}

void EquatorialHotSpot::getVelocity(double const pos[4], double vel[4]) {
  double coord_spot[4]={pos[0]};
  const_cast<EquatorialHotSpot*>(this)
    ->getCoord(coord_spot, 1, coord_spot+1, coord_spot+2, coord_spot+3);
  gg_ -> circularVelocity(coord_spot, vel, dir_);
}

double EquatorialHotSpot::emission(double nu_em, double dsem,
				   state_t const &coord_ph,
				   double const coord_obj[8]) const{
  double coord_spot[4]={coord_obj[0]};
  const_cast<EquatorialHotSpot*>(this)
    ->getCartesian(coord_spot, 1, coord_spot+1, coord_spot+2, coord_spot+3);
  //above: nasty trick to deal with constness of emission
  double xspot=coord_spot[1], yspot=coord_spot[2];
  double rr=coord_obj[1], phi=coord_obj[3];
  double difx=(rr*cos(phi)-xspot),
    dify=(rr*sin(phi)-yspot);
  double d2 = difx*difx+dify*dify;
  double ds2=sizespot_*sizespot_;
  if (d2 < 16*ds2){ // we are within 4*rspot,
                    // same as in Schnittman & Bertschinger 2004

    // computing the angle (normal,photon tangent)
    double cosalpha=0.;
    if (beaming_ == NormalBeaming or beaming_ == RadialBeaming){
      double gthth=gg_->gmunu(&coord_ph[0],2,2);
      double pth=coord_ph[6];
      double uemitter[4];
      const_cast<EquatorialHotSpot*>(this)
	->getVelocity(&coord_ph[0],uemitter);
      double pscalu=fabs(gg_->ScalarProd(&coord_ph[0],&coord_ph[4],
					 uemitter));
      if (pscalu==0.) GYOTO_ERROR("Undefined cosalpha!");
      cosalpha = 1./pscalu*sqrt(gthth)*fabs(pth); // = |cos(alpha)|

      if (fabs(cosalpha)>1.)
	GYOTO_ERROR("cosalpha>1!");
    }
      
    // emission Gaussian width
    double sigma2=ds2 ; // following choice of Schnittman & Bertschinger 2004:
                        // sigma = Rspot

    switch (beaming_) {
    case IsotropicBeaming:
      return exp(-d2/(2*sigma2));
    case IsotropicConstant:
      return 1.;
    case NormalBeaming:
      return cosalpha*cosalpha*exp(-d2/(2*sigma2));
    case RadialBeaming:
      return (1.-cosalpha)*(1.-cosalpha)*exp(-d2/(2*sigma2));
    default:
      GYOTO_ERROR("In EquatorialHotSpot::emission:"
		 " incorrect beaming argument");
    }
  }
  // else
  return 0.;
}
