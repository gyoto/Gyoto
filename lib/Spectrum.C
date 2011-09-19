/*
    Copyright 2011 Thibaut Paumard

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

#include "GyotoSpectrum.h"
#include "GyotoRegister.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"

#include <cmath>
#include <iostream>
using namespace Gyoto;
using namespace std;

Spectrum::Generic::Generic(const string kind) : kind_(kind) {}
Spectrum::Generic * Spectrum::Generic::clone() const 
{
  string msg = "Spectrum::clone() called: "
    "cloning unimplemented for Spectrum kind ";
  msg += kind_;
  throwError(msg);
  return const_cast<Spectrum::Generic*>(this);
              // avoid warning, we won't get to that point
}
const string Spectrum::Generic::getKind() const { return kind_; }

double Spectrum::Generic::integrate(double nu1, double nu2) {
  double nu;

  if (nu1>nu2) {nu=nu1; nu1=nu2; nu2=nu;}

  double Inu1 = operator()(nu1), Inu2=operator()(nu2);
  double dnux2 = ((nu2-nu1)*2.);
  double Icur = (Inu2+Inu1)*dnux2*0.25;
  double Iprev;

  if (debug())
      cerr << "DEBUG: Spectrum::Generic::integrate(): "
	   << "Icur=" << Icur << endl;

  do {
    Iprev = Icur; 
    dnux2 *= 0.5;
    for (nu = nu1 + 0.5*dnux2; nu < nu2; nu += dnux2) {
      Icur += operator()(nu) * dnux2;
    }
    Icur *= 0.5;
    if (debug())
      cerr << "DEBUG: Spectrum::Generic::integrate(): "
	   << "Icur=" << Icur << endl;
  } while( fabs(Icur-Iprev) > (1e-2 * Icur) );

  if (debug())
    cerr << "DEBUG: Spectrum::Generic::integrate(): "
	 << "dnu=" << dnux2*0.5
	 << "=(nu2-nu1)/" << (nu2-nu1)/(dnux2*0.5) << endl;

  return Icur;
}


double Spectrum::Generic::integrate(double nu1, double nu2,
				    const Spectrum::Generic *  opacity,
				    double dsem) {
  double nu;

  if (nu1>nu2) {nu=nu1; nu1=nu2; nu2=nu;}

  double Inu1 = operator()(nu1, (*opacity)(nu1), dsem),
    Inu2=operator()(nu2, (*opacity)(nu2), dsem);
  double dnux2 = ((nu2-nu1)*2.);
  double Icur = (Inu2+Inu1)*dnux2*0.25;
  double Iprev;

  if (debug())
      cerr << "DEBUG: Spectrum::Generic::integrate(): "
	   << "Icur=" << Icur << endl;

  do {
    Iprev = Icur; 
    dnux2 *= 0.5;
    for (nu = nu1 + 0.5*dnux2; nu < nu2; nu += dnux2) {
      Icur += operator()(nu, (*opacity)(nu), dsem) * dnux2;
    }
    Icur *= 0.5;
    if (debug())
      cerr << "DEBUG: Spectrum::Generic::integrate(): "
	   << "Icur=" << Icur << endl;
  } while( fabs(Icur-Iprev) > (1e-2 * Icur) );

  if (debug())
    cerr << "DEBUG: Spectrum::Generic::integrate(): "
	 << "dnu=" << dnux2*0.5
	 << "=(nu2-nu1)/" << (nu2-nu1)/(dnux2*0.5) << endl;

  return Icur;
}

double Spectrum::Generic::operator()(double nu, double opacity, double ds)
  const {
  double thickness;
  // assume emissivity = opacity
  if ((thickness=(opacity*ds))) return operator()(nu) * (1. - exp (-thickness)) ;
  return 0.;
}

#ifdef GYOTO_USE_XERCES
// do nothing... for now
void Spectrum::Generic::fillElement(FactoryMessenger *fmp ) const {
  fmp->setSelfAttribute("kind", kind_);
}
void Spectrum::Generic::setGenericParameter(std::string, std::string) {}

Register::Entry* Gyoto::Spectrum::Register_ = NULL;

void Gyoto::Spectrum::Register(std::string name,
			       Gyoto::Spectrum::Subcontractor_t* scp)
{
  Register::Entry* ne =
    new Register::Entry(name,
			(Gyoto::SmartPointee::Subcontractor_t*)scp,
			Gyoto::Spectrum::Register_);
  Gyoto::Spectrum::Register_ = ne;
}

Spectrum::Subcontractor_t* Spectrum::getSubcontractor(std::string name) {
  if (!Gyoto::Spectrum::Register_) throwError("No Spectrum kind registered!");
  return (Spectrum::Subcontractor_t*)Gyoto::Spectrum::Register_
    -> getSubcontractor(name);
}

void Gyoto::Spectrum::initRegister() {
  if (Gyoto::Spectrum::Register_) delete Gyoto::Spectrum::Register_;
  Gyoto::Spectrum::Register_ = NULL;
}

#endif

