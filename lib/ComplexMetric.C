/*
    Copyright 2020 Thibaut Paumard

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

/*
  TODO: what about mass?
 */


#include "GyotoComplexMetric.h"
#include "GyotoFactoryMessenger.h"

using namespace Gyoto;
using namespace Gyoto::Metric;

Complex::Complex() :
  WIP("Gyoto::Metric::Complex"),
  Generic(GYOTO_COORDKIND_UNSPECIFIED, "Complex"),
  cardinal_(0),
  elements_(NULL)
{}

Complex::Complex(const Complex& o) :
  Metric::Generic(o),
  cardinal_(o.cardinal_),
  elements_(NULL)
{
  coordKind(o.coordKind());
  if (cardinal_) {
    elements_ = new SmartPointer<Generic> [cardinal_];
    for (size_t i=0; i< cardinal_; ++i) {
      elements_[i] = o[i]->clone();
    }
  }
}
Complex *Complex::clone() const {return new Complex(*this); }

Complex::~Complex()
{
  if (cardinal_) for (size_t i=0; i< cardinal_; ++i) elements_[i] = NULL;
}

bool Complex::isThreadSafe() const {
  bool safe = Generic::isThreadSafe();
  for (size_t i=0; i < cardinal_; ++i) safe &= elements_[i] -> isThreadSafe();
  return safe;
}

void Complex::append(SmartPointer<Generic> e)
{
  GYOTO_DEBUG << std::endl;
  if (cardinal_+1 == 0) GYOTO_ERROR("Complex::append(): OVERFLOW");
  if (cardinal_ && (e->coordKind() != coordKind())) GYOTO_ERROR("inconsistent coord kind");
  SmartPointer<Generic> * orig = elements_;
  elements_ = new SmartPointer<Generic> [cardinal_+1];
  for (size_t i=0; i< cardinal_; ++i) {
    elements_[i] = orig[i];
    orig[i] = NULL;
  }
  delete [] orig; orig = NULL;
  elements_[cardinal_] = e;
  ++cardinal_;
  coordKind(e->coordKind());
  GYOTO_DEBUG << "done" << std::endl;
}

SmartPointer<Generic>& Complex::operator[](size_t i)
{
  if (i >= cardinal_)
    GYOTO_ERROR("Complex::operator[](size_t i): no such element");
  return elements_[i];
}

SmartPointer<Generic> const& Complex::operator[](size_t i) const
{
  if (i >= cardinal_)
    GYOTO_ERROR("Complex::operator[](size_t i): no such element");
  return elements_[i];
}

void Complex::remove(size_t i) {
  if (i >= cardinal_)
    GYOTO_ERROR("Complex::remove(size_t i): no such element");
  SmartPointer<Generic> * orig = elements_;
  if (--cardinal_) elements_ = new SmartPointer<Generic> [cardinal_];
  else elements_ = NULL;
  size_t k, j=0;
  for (k=0; k<= cardinal_; ++k) {
    if (k != i) elements_[j++] = orig[k];
    orig[k] = NULL;
  }
  delete [] orig;
  if (!cardinal_) coordKind(GYOTO_COORDKIND_UNSPECIFIED);
}

size_t Complex::getCardinal() const {return cardinal_; }

void Complex::gmunu(double g[4][4], const double x[4]) const {
  size_t mu, nu, k;
  double cur[4][4];
  for (mu=0; mu<4 ; ++mu)
    for (nu=0; nu<4 ; ++nu)
      g[mu][nu] = 0.;
  for (k=0; k<cardinal_; ++k) {
    elements_[k]->gmunu(cur, x);
    for (mu=0; mu<4 ; ++mu)
      for (nu=0; nu<4 ; ++nu)
	g[mu][nu] += cur[mu][nu];
  }
}

void Complex::jacobian(double jac[4][4][4], const double x[4]) const {
  size_t a, mu, nu, k;
  double cjac[4][4][4];
  for (mu=0; mu<4 ; ++mu)
    for (nu=0; nu<4 ; ++nu)
      for (a=0; a<4; ++a) jac[a][mu][nu] = 0.;
  for (k=0; k<cardinal_; ++k) {
    elements_[k]->jacobian(cjac, x);
    for (mu=0; mu<4 ; ++mu)
      for (nu=0; nu<4 ; ++nu)
	for (a=0; a<4; ++a) jac[a][mu][nu] += cjac[a][mu][nu];
  }
}

#ifdef GYOTO_USE_XERCES
void Complex::fillElement(FactoryMessenger *fmp) const {
  FactoryMessenger * childfmp=NULL;

  for (size_t i=0; i<cardinal_; ++i) {
    childfmp = fmp -> makeChild ( "SubMetric" );
    elements_[i] -> fillElement(childfmp);
    delete childfmp;
  }

  Metric::Generic::fillElement(fmp);
}

void Complex::setParameters(FactoryMessenger *fmp) {
  GYOTO_DEBUG << std::endl;

  std::string name="", content="", unit="";
  std::vector<std::string> plugin;
  FactoryMessenger * child = NULL;

  while (fmp->getNextParameter(&name, &content, &unit)) {
    GYOTO_DEBUG_EXPR(name);
    if (name=="SubMetric") {
      content = fmp -> getAttribute("kind");
      plugin  = split(fmp -> getAttribute("plugin"), ",");
      child   = fmp -> getChild();
      append ((*Metric::getSubcontractor(content, plugin))(child, plugin));
      delete child;
    } else setParameter(name, content, unit);
  }

  GYOTO_DEBUG << "done" << std::endl;
}
#endif

int Complex::isStopCondition(double const coord[8]) const {
  for (size_t k=0; k<cardinal_; ++k)
    if (elements_[k]-> isStopCondition(coord))
      return 1;
  return 0;
}
