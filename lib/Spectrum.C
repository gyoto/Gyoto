/*
    Copyright 2011-2016 Thibaut Paumard

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

/// Properties

// There is no generic properties for spectra. Nevertheless, we define
// this to derived classes can point to Spectrum::Generic::properties
// rather than Object::properties

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(Spectrum::Generic)
GYOTO_PROPERTY_END(Spectrum::Generic, Object::properties)

///

Spectrum::Generic::Generic()
: SmartPointee(), Object() {}
Spectrum::Generic::Generic(const string kin)
: SmartPointee(), Object(kin) {}
Spectrum::Generic::Generic(const Generic& o)
: SmartPointee(o), Object(o) {}
Spectrum::Generic * Spectrum::Generic::clone() const 
{
  string msg = "Spectrum::clone() called: "
    "cloning unimplemented for Spectrum kind ";
  msg += kind_;
  GYOTO_ERROR(msg);
  return const_cast<Spectrum::Generic*>(this);
              // avoid warning, we won't get to that point
}
Spectrum::Generic::~Generic() { GYOTO_DEBUG << endl; }

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

GYOTO_GETSUBCONTRACTOR(Spectrum)

void Gyoto::Spectrum::initRegister() {
  if (Gyoto::Spectrum::Register_) delete Gyoto::Spectrum::Register_;
  Gyoto::Spectrum::Register_ = NULL;
}


