/*
    Copyright 2013-2014, 2016 Thibaut Paumard

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

#include "GyotoComplexSpectrometer.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoMetric.h"
#include "GyotoPhoton.h"

#include <cstring>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Spectrometer;

kind_t const Complex::Kind = "Complex";

Complex::Complex() :
  Generic(Complex::Kind),
  cardinal_(0),
  elements_(NULL)
{

}

Complex::Complex(const Complex& o) :
  Generic(o),
  cardinal_(o.cardinal_),
  elements_(NULL)
{
  if (cardinal_) {
    elements_ = new SmartPointer<Generic> [cardinal_];
    for (size_t i=0; i< cardinal_; ++i) {
      elements_[i] = o[i]->clone();
    }
  }
}
Complex *Complex::clone() const {return new Complex(*this); }

bool Complex::isThreadSafe() const {
  bool safe = Generic::isThreadSafe();
  for (size_t i=0; i < cardinal_; ++i) safe &= elements_[i] -> isThreadSafe();
  return safe;
}

Complex::~Complex()
{
  if (cardinal_) for (size_t i=0; i< cardinal_; ++i) elements_[i] = NULL;
}

void Complex::append(SmartPointer<Generic> e)
{
  GYOTO_DEBUG << "DEBUG: in Complex::append(SmartPointer<Generic> e)" << endl;
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
  e->hook(this);
  GYOTO_DEBUG << "DEBUG: out Complex::append(SmartPointer<Generic> e)" << endl;
  tell(this);
}

SmartPointer<Generic> & Complex::operator[](size_t i)
{
  if (i > cardinal_)
    GYOTO_ERROR("Complex::operator[](size_t i): no such element");
  return elements_[i];
}

SmartPointer<Generic> const & Complex::operator[](size_t i) const
{
  if (i > cardinal_)
    GYOTO_ERROR("Complex::operator[](size_t i): no such element");
  return elements_[i];
}

void Complex::remove(size_t i) {
  if (i >= cardinal_)
    GYOTO_ERROR("Complex::remove(size_t i): no such element");
  elements_[i]->unhook(this);
  SmartPointer<Generic> * orig = elements_;
  if (--cardinal_) elements_ = new SmartPointer<Generic> [cardinal_];
  else elements_ = NULL;
  size_t k, j=0;
  for (k=0; k<= cardinal_; ++k) {
    if (k != i) elements_[j++] = orig[k];
    orig[k] = NULL;
  }
  delete [] orig;
  tell(this); // tell self!
}

size_t Complex::getCardinal() const {return cardinal_; }

#ifdef GYOTO_USE_XERCES
void Complex::fillElement(FactoryMessenger *fmp) const {
  FactoryMessenger * childfmp=NULL;

  for (size_t i=0; i<cardinal_; ++i) {
    childfmp = fmp -> makeChild ( "SubSpectrometer" );
    elements_[i] -> fillElement(childfmp);
    delete childfmp;
  }

  Spectrometer::Generic::fillElement(fmp);
}

void Complex::setParameters(FactoryMessenger *fmp) {
  if (debug())
    cerr << "DEBUG: in Complex::setParameters()" << endl;

  string name="", content="", unit="";
  std::vector<std::string> plugin;
  FactoryMessenger * child = NULL;

  while (fmp->getNextParameter(&name, &content, &unit)) {
    if (debug())
      cerr << "DEBUG: Spectrometer::Complex::Subcontractor(): name=" << name << endl;
    if (name=="SubSpectrometer") {
      content = fmp -> getAttribute("kind");
      child   = fmp -> getChild();
      plugin  = split(fmp -> getAttribute("plugin"), ",");
      append ((*Spectrometer::getSubcontractor(content, plugin))(child, plugin));
      delete child;
    } else setParameter(name, content, unit);
  }

  if (debug())
    cerr << "DEBUG: out Complex::setParameters()" << endl;
}
#endif

void Complex::tell(Gyoto::Hook::Teller *) {
  // This is called each time an element is added or mutated
  // This is suboptimal, but most straightforward
  nboundaries_=nsamples_=0;
  for (size_t i=0; i<cardinal_; ++i) {
    nsamples_ += elements_[i]->nSamples();
    nboundaries_ += elements_[i]->getNBoundaries();
  }
  if (boundaries_) delete [] boundaries_;
  if (widths_)     delete [] widths_;
  if (midpoints_)  delete [] midpoints_;
  if (chanind_)    delete [] chanind_;
  boundaries_ = new double [nboundaries_];
  widths_     = new double [nsamples_];
  midpoints_  = new double [nsamples_];
  chanind_    = new size_t [2*nsamples_];
  
  size_t boffset=0, offset=0;
  size_t const * chanind=0;
  for (size_t i=0; i<cardinal_; ++i) {
    size_t enb = elements_[i]->getNBoundaries();
    size_t ens = elements_[i]->nSamples();
    memcpy(boundaries_+boffset,
	   elements_[i]->getChannelBoundaries(),
	   enb*sizeof(double));
    memcpy(widths_+offset,
	   elements_[i]->getWidths(),
	   ens*sizeof(double));
    memcpy(midpoints_+offset,
	   elements_[i]->getMidpoints(),
	   ens*sizeof(double));
    chanind=elements_[i]->getChannelIndices();
    for (size_t j=0; j<2*ens; ++j)
      chanind_[2*offset+j]=chanind[j]+boffset;
    boffset += enb;
    offset  += ens ;
  }
  tellListeners();
}
