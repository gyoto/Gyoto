/*
  Copyright 2017 Frederic Vincent, Thibaut Paumard

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
#include "GyotoNeutronStarAnalyticEmission.h"
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
GYOTO_PROPERTY_START(NeutronStarAnalyticEmission,
		     "Neutron star emitting at its surface.")
GYOTO_PROPERTY_SPECTRUM(NeutronStarAnalyticEmission, Spectrum, spectrum,
			"Emission law.")
GYOTO_PROPERTY_END(NeutronStarAnalyticEmission, NeutronStar::properties)

NeutronStarAnalyticEmission::NeutronStarAnalyticEmission() :
NeutronStar("NeutronStarAnalyticEmission"), spectrum_(NULL) {
  GYOTO_DEBUG << endl;
}

NeutronStarAnalyticEmission::NeutronStarAnalyticEmission(const NeutronStarAnalyticEmission& o) :
  NeutronStar(o), spectrum_(NULL) {
  GYOTO_DEBUG << endl;
  if (o.spectrum_()) spectrum_ = o.spectrum_->clone();
}
NeutronStarAnalyticEmission * NeutronStarAnalyticEmission::clone() const {
  return new NeutronStarAnalyticEmission(*this); }

NeutronStarAnalyticEmission::~NeutronStarAnalyticEmission() {
  GYOTO_DEBUG << endl;
}

SmartPointer<Spectrum::Generic> NeutronStarAnalyticEmission::spectrum() const { return spectrum_; }
  void NeutronStarAnalyticEmission::spectrum(SmartPointer<Spectrum::Generic> sp) {spectrum_=sp;}

double NeutronStarAnalyticEmission::emission(double nu_em, double, state_t const &,
			      double const *) const{
  GYOTO_DEBUG << endl;
  if (flag_radtransf_)
    GYOTO_ERROR("Radiative transfer not implemented for NeutronStarAnalyticEmission.");

  return (*spectrum_)(nu_em);
}

