/*
    Copyright 2011 Thibaut Paumard, Frederic Vincent

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
  setMetric(gg_); // to set the same metric in all elements
}
Complex *Complex::clone() const {return new Complex(*this); }


Complex::~Complex()
{
  if (cardinal_) for (size_t i=0; i< cardinal_; ++i) elements_[i] = NULL;
}

void Complex::setMetric(SmartPointer<Metric::Generic> gg)
{
  Generic::setMetric(gg);
  for (size_t i=0; i<cardinal_; ++i) {
    if (debug()) {
      cerr << "DEBUG: Complex::setMetric(gg): ";
      cerr << "elements_["<<i<<"] is a ";
      cerr << elements_[i]->getKind();
      cerr << ". Setting metric." << endl;
    }
    elements_[i]->setMetric(gg_);
  }
}

void Complex::append(SmartPointer<Generic> e)
{
  if (debug())
    cerr << "DEBUG: in Complex::append(SmartPointer<Generic> e)" << endl;
  if (cardinal_+1 == 0) throwError("Complex::append(): OVERFLOW");
  SmartPointer<Generic> * orig = elements_;
  elements_ = new SmartPointer<Generic> [cardinal_+1];
  for (size_t i=0; i< cardinal_; ++i) {
    elements_[i] = orig[i];
    orig[i] = NULL;
  }
  delete [] orig; orig = NULL;
  elements_[cardinal_] = e;
  ++cardinal_;
  if (gg_) e->setMetric(gg_);
  else gg_ = e->getMetric();
  if (debug())
    cerr << "DEBUG: out Complex::append(SmartPointer<Generic> e)" << endl;
}

SmartPointer<Generic> Complex::operator[](size_t i)
{
  if (i > cardinal_)
    throwError("Complex::operator[](size_t i): no such element");
  return elements_[i];
}

SmartPointer<Generic> const Complex::operator[](size_t i) const
{
  if (i > cardinal_)
    throwError("Complex::operator[](size_t i): no such element");
  return elements_[i];
}

void Complex::remove(size_t i) {
  if (i >= cardinal_)
    throwError("Complex::remove(size_t i): no such element");
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
    Photon refine (ph, index, 1, step_max_);
    double trans_orig = refine . getTransmission(size_t(-1));
    size_t n_refine = refine . get_nelements();
    if (debug())
      cerr << "DEBUG: Complex::Impact(...): n_refine=="<<n_refine << endl;
    for (size_t n=n_refine-2; n!=size_t(-1); --n) {
      for (size_t i=0; i<cardinal_; ++i)
	if (impact[i]) {
	  if (debug())
	    cerr << "DEBUG: Complex::Impact(...): calling Impact for elements_["
		 << i << "] (" << elements_[i]->getKind() << ")" << endl;
	  elements_[i]->Impact(&refine, n, data);
	}
    }
    ph->transmit(size_t(-1),
		 trans_orig?refine . getTransmission(size_t(-1))/trans_orig:0.);
  }

  delete [] impact;
  return res;
}

#ifdef GYOTO_USE_XERCES
void Complex::fillElement(FactoryMessenger *fmp) const {
  FactoryMessenger * childfmp=NULL;

  fmp -> setMetric (getMetric()) ;

  for (size_t i=0; i<cardinal_; ++i) {
    childfmp = fmp -> makeChild ( "SubAstrobj" );
    elements_[i] -> fillElement(childfmp);
    delete childfmp;
  }

  Astrobj::Generic::fillElement(fmp);
}

SmartPointer<Astrobj::Generic> Gyoto::Astrobj::Complex::Subcontractor(FactoryMessenger* fmp) {

  if (debug())
    cerr << "DEBUG: in Complex::Subcontractor()" << endl;

  string name="", content="";
  FactoryMessenger * child = NULL;
  SmartPointer<Complex> cplx = new Complex();

  cplx -> setMetric( fmp->getMetric() );

  while (fmp->getNextParameter(&name, &content)) {
    if (debug())
      cerr << "DEBUG: Astrobj::Complex::Subcontractor(): name=" << name << endl;
    if (name=="SubAstrobj") {
      content = fmp -> getAttribute("kind");
      child = fmp -> getChild();
      cplx -> append ((*Astrobj::getSubcontractor(content))(child));
      delete child;
    }
    else
      cplx -> setGenericParameter(name, content);
  }

  if (debug())
    cerr << "DEBUG: out Complex::Subcontractor()" << endl;
  return cplx;
}

void Gyoto::Astrobj::Complex::Init() {
  Gyoto::Astrobj::Register("Complex", &Gyoto::Astrobj::Complex::Subcontractor);
}
#endif
