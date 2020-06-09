/*
    Copyright 2011-2014, 2016 Thibaut Paumard

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

#include "GyotoComplexAstrobj.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoMetric.h"
#include "GyotoPhoton.h"

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

Complex::Complex() :
  Generic("Complex"),
  cardinal_(0),
  elements_(NULL),
  step_max_(GYOTO_DEFAULT_DELTA)
{
  rmax_=0.; // By default rMax is computed from the sub-astrobjs
}

Complex::Complex(const Complex& o) :
  Astrobj::Generic(o),
  cardinal_(o.cardinal_),
  elements_(NULL),
  step_max_(o.step_max_)
{
  if (cardinal_) {
    elements_ = new SmartPointer<Generic> [cardinal_];
    for (size_t i=0; i< cardinal_; ++i) {
      elements_[i] = o[i]->clone();
    }
  }
  metric(gg_); // to set the same metric in all elements
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

void Complex::metric(SmartPointer<Metric::Generic> gg)
{
  Generic::metric(gg);
  for (size_t i=0; i<cardinal_; ++i) {
    if (debug()) {
      cerr << "DEBUG: Complex::metric(gg): ";
      cerr << "elements_["<<i<<"] is a ";
      cerr << elements_[i]->kind();
      cerr << ". Setting metric." << endl;
    }
    elements_[i]->metric(gg_);
  }
}

void Complex::append(SmartPointer<Generic> e)
{
  if (debug())
    cerr << "DEBUG: in Complex::append(SmartPointer<Generic> e)" << endl;
  if (cardinal_+1 == 0) GYOTO_ERROR("Complex::append(): OVERFLOW");
  SmartPointer<Generic> * orig = elements_;
  elements_ = new SmartPointer<Generic> [cardinal_+1];
  for (size_t i=0; i< cardinal_; ++i) {
    elements_[i] = orig[i];
    orig[i] = NULL;
  }
  delete [] orig; orig = NULL;
  elements_[cardinal_] = e;
  ++cardinal_;
  if (gg_) e->metric(gg_);
  else gg_ = e->metric();
  if (debug())
    cerr << "DEBUG: out Complex::append(SmartPointer<Generic> e)" << endl;
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
}

size_t Complex::getCardinal() const {return cardinal_; }

double Complex::deltaMax(double coord[8]) {
  double h1max=DBL_MAX, tmp;
  for (size_t i=0; i<cardinal_; ++i)
    if (h1max> (tmp=elements_[i]->deltaMax(coord))) h1max=tmp;
  return h1max;
}

double Complex::rMax() {
  double rmax = Generic::rMax(), rmaxnew=0.;
  for (size_t i=0; i<cardinal_; ++i)
    if (rmax < (rmaxnew=elements_[i] -> rMax())) rmax=rmaxnew;
  return rmax;
}

int Complex::Impact(Photon* ph, size_t index, Properties *data)
{
  int res=0, *impact = new int[cardinal_];
  size_t n_impact = 0;
  for (size_t i=0; i<cardinal_; ++i)
    n_impact += impact[i] = elements_[i] -> Impact(ph, index, NULL);

  if (debug())
    cerr << "DEBUG: Complex::Impact(...): " <<n_impact <<" sub-impacts" << endl;

  if (n_impact==1) {
    res = 1;
    for (size_t i=0; i<cardinal_; ++i)
      if (impact[i])
	elements_[i] -> Impact(ph, index, data);
  } else if (n_impact >= 2) {
    res = 1;
    if (debug())
      cerr << "DEBUG: Complex::Impact(...): refining Photon" << endl;
    Photon::Refined refine (ph, index, 1, step_max_);
    size_t n_refine = refine . get_nelements();
    if (debug())
      cerr << "DEBUG: Complex::Impact(...): n_refine=="<<n_refine << endl;
    for (size_t n=n_refine-2; n!=size_t(-1); --n) {
      for (size_t i=0; i<cardinal_; ++i)
	if (impact[i]) {
	  if (debug())
	    cerr << "DEBUG: Complex::Impact(...): calling Impact for elements_["
		 << i << "] (" << elements_[i]->kind() << ")" << endl;
	  elements_[i]->Impact(&refine, n, data);
	}
    }
  }

  delete [] impact;
  return res;
}

#ifdef GYOTO_USE_XERCES
void Complex::fillElement(FactoryMessenger *fmp) const {
  FactoryMessenger * childfmp=NULL;

  fmp -> metric (metric()) ;

  for (size_t i=0; i<cardinal_; ++i) {
    childfmp = fmp -> makeChild ( "SubAstrobj" );
    elements_[i] -> fillElement(childfmp);
    delete childfmp;
  }

  Astrobj::Generic::fillElement(fmp);
}

void Complex::setParameters(FactoryMessenger *fmp) {
  if (debug())
    cerr << "DEBUG: in Complex::setParameters()" << endl;

  string name="", content="", unit="";
  vector<string> plugin;
  FactoryMessenger * child = NULL;

  metric( fmp->metric() );

  while (fmp->getNextParameter(&name, &content, &unit)) {
    if (debug())
      cerr << "DEBUG: Astrobj::Complex::Subcontractor(): name=" << name << endl;
    if (name=="SubAstrobj") {
      content = fmp -> getAttribute("kind");
      plugin  = split(fmp -> getAttribute("plugin"), ",");
      child   = fmp -> getChild();
      append ((*Astrobj::getSubcontractor(content, plugin))(child, plugin));
      delete child;
    } else setParameter(name, content, unit);
  }

  if (debug())
    cerr << "DEBUG: out Complex::setParameters()" << endl;
}
#endif
