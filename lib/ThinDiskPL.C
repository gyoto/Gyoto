/*
    Copyright 2011 Frederic Vincent, Thibaut Paumard

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
#include "GyotoThinDiskPL.h"
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

ThinDiskPL::ThinDiskPL() :
  ThinDisk("ThinDiskPL"),
  PLSlope_(0.), PLRho_(1.), PLRadRef_(1.),
  spectrumBB_(NULL)
{
  if (debug()) cerr << "DEBUG: ThinDiskPL Construction" << endl;
  spectrumBB_ = new Spectrum::BlackBody(); 
}

ThinDiskPL::ThinDiskPL(const ThinDiskPL& o) :
  ThinDisk(o),
  PLSlope_(o.PLSlope_), PLRho_(o.PLRho_), PLRadRef_(o.PLRadRef_),
  spectrumBB_(NULL)
{
  if (o.gg_()) gg_=o.gg_->clone();
  if (o.spectrumBB_()) spectrumBB_=o.spectrumBB_->clone();
  Generic::gg_=gg_;
}
ThinDiskPL* ThinDiskPL::clone() const
{ return new ThinDiskPL(*this); }

ThinDiskPL::~ThinDiskPL() {
  if (debug()) cerr << "DEBUG: ThinDiskPL Destruction" << endl;
}

double ThinDiskPL::emission(double nu, double,
				    double *,
				    double coord_obj[8]) const{

  //cout << "param= " << PLSlope_ << " " << PLRho_ << " " << PLRadRef_ 
  //    << " " << getInnerRadius() << " " << getOuterRadius() << endl;

  double Iem = emissionBB(nu,coord_obj);
  /* Only BB emission so far, other ways of
   computing emission from rho can be added */
  
  /*if (flag_radtransf_)
    throwError("In ThinDiskPL::emission() "
    "optically thin integration not provided");*/

  return Iem;

}

double ThinDiskPL::emissionBB(double nu,
			      double co[8]) const{

  double rcur=projectedRadius(co);
  double rho_si = PLRho_*pow(rcur/PLRadRef_,PLSlope_);
  //Assuming: pressure = kappa*(mass density)^gamma, gamma=5/3 (eq of state)
  // and: pressure = (mass density)*R/Mm*T (Mm = molar mass)
  
  double gamma=5./3.;
  double Mm=6e-4;//Navogadro*Munit/gamma
  double kappa=3e10;//pressure coef: p = kappa*rho^gamma, rho=mass density
  
  double cs2=kappa*gamma*pow(rho_si,gamma-1.);
  double TT=Mm/GYOTO_GAS_CST*cs2;//Temperature in SI
  //cout << "TT after rl= " << TT << endl;
  //cout << "r,rho,T= " << rcross << " " << rho_si << " " << TT << endl;
  spectrumBB_->setTemperature(TT);
  return (*spectrumBB_)(nu);
}

int ThinDiskPL::setParameter(std::string name,
			     std::string content,
			     std::string unit) {
  if      (name=="PLSlope") PLSlope_=atof(content.c_str());
  else if (name=="PLRho") PLRho_=atof(content.c_str());
  else if (name=="PLRadRef") PLRadRef_=atof(content.c_str());
  else if (name=="Rmin") setInnerRadius(atof(content.c_str()));
  else if (name=="Rmax") setOuterRadius(atof(content.c_str()));
  else return ThinDisk::setParameter(name, content, unit);
  return 0;
}


#ifdef GYOTO_USE_XERCES
void ThinDiskPL::fillElement(FactoryMessenger *fmp) const {
  if (PLSlope_) fmp->setParameter("PLSlope", PLSlope_);
  if (PLRho_) fmp->setParameter("PLRho", PLRho_);
  if (PLRadRef_) fmp->setParameter("PLRadRef", PLRadRef_);
  ThinDisk::fillElement(fmp);
}
#endif
