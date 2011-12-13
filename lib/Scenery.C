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

#include "GyotoUtils.h"
#include "GyotoScenery.h"
#include "GyotoPhoton.h"
#include "GyotoFactoryMessenger.h"

#include <cmath>
#include <cfloat>
#include <cstring>
#include <cstdlib>

#define DEFAULT_TLIM 0.

using namespace Gyoto;
using namespace std;

/*Scenery::Scenery() :
  gg_(NULL), screen_(NULL), obj_(NULL),
  deltatau_(0.01) {}
*/
Scenery::Scenery() :
  gg_(NULL), screen_(NULL), obj_(NULL), delta_(0.01),
  quantities_(0), ph_(), tlim_(DEFAULT_TLIM) {}

Scenery::Scenery(const Scenery& o) :
  SmartPointee(o),
  gg_(NULL), screen_(NULL), obj_(NULL), delta_(o.delta_),
  quantities_(o.quantities_), ph_(o.ph_), tlim_(o.tlim_)
{
  // We have up to 3 _distinct_ clones of the same Metric.
  // Keep only one.
  if (o.gg_()) gg_=o.gg_->clone();
  if (o.screen_()) {
    screen_=o.screen_->clone();
    screen_->setMetric(gg_);
  }
  if (o.obj_()) {
    obj_=o.obj_->clone();
    obj_->setMetric(gg_);
  }
}
Scenery * Scenery::clone() const { return new Scenery(*this); }

/*Scenery::Scenery(SmartPointer<Metric::Generic> met, SmartPointer<Screen> screen, SmartPointer<Astrobj::Generic> obj) :
  gg_(met), screen_(screen), obj_(obj),
  deltatau_(0.01)
{}
*/
Scenery::Scenery(SmartPointer<Metric::Generic> met, SmartPointer<Screen> screen, SmartPointer<Astrobj::Generic> obj) :
  gg_(met), screen_(screen), obj_(obj), delta_(0.01),
  quantities_(0)
{
  if (screen_) screen_->setMetric(gg_);
  if (obj_) obj_->setMetric(gg_);
}

Scenery::~Scenery() {
  if (debug())
    cerr << "DEBUG: in Scenery::~Scenery()\n"
	 << "DEBUG: Scenery::~Scenery(): freeing metric\n";
  gg_ = NULL;

  if (debug())
    cerr << "DEBUG: Scenery::~Scenery(): freeing screen\n";
  screen_ = NULL;

  if (debug())
    cerr << "DEBUG: Scenery::~Scenery(): freeing astrobj\n";
  obj_ = NULL;
 }

SmartPointer<Metric::Generic> Scenery::getMetric() { return gg_; }

void Scenery::setMetric(SmartPointer<Metric::Generic> met) {
  gg_ = met;
  if (!screen_) screen_ = new Screen ();
  screen_ -> setMetric(gg_);
  if (obj_) obj_ -> setMetric(gg_);
}

SmartPointer<Screen> Scenery::getScreen() { return screen_; }

void Scenery::setScreen(SmartPointer<Screen> screen) {
  screen_ = screen;
  if (gg_) screen_ -> setMetric (gg_) ;
}

SmartPointer<Astrobj::Generic> Scenery::getAstrobj() { return obj_; }
void Scenery::setAstrobj(SmartPointer<Astrobj::Generic> obj) {
  obj_ = obj;
  if (gg_) obj_ -> setMetric (gg_) ;
}

double Scenery::getDelta() const { return delta_; }
void Scenery::setDelta(double d) { delta_ = d; }

void Scenery::rayTrace(size_t imin, size_t imax,
		       size_t jmin, size_t jmax,
		       Astrobj::Properties *data,
		       double * impactcoords) {

  //  if (debug()) cout << "screen dist beg ray trace= " << screen_ -> getDistance() << endl;

  const size_t npix = screen_->getResolution();
  //ofstream pixels(filename);//illuminated pixels on screen
  SmartPointer<Spectrometer> spr = screen_->getSpectrometer();
  ph_.setSpectrometer(spr);
  size_t nbnuobs = spr() ? spr -> getNSamples() : 0;
  ph_.setTlim(tlim_);

  double coord[8];

  screen_->computeBaseVectors();
         // Necessary for KS integration, computes relation between
         // observer's x,y,z coord and KS X,Y,Z coord. Will be used to
         // compute photon's initial tangent vector.
     // Note : this is a BUG if this is required, should be done automagically.

  imax=(imax<=(npix)?imax:(npix));
  jmax=(jmax<=(npix)?jmax:(npix));

  for (size_t j=jmin;j<=jmax;j++) {

    if (verbose() >= GYOTO_QUIET_VERBOSITY && !impactcoords)
      cout << "\rj = " << j << " / " << jmax << " " << flush;

    for (size_t i=imin;i<=imax;i++) {
      ph_.setDelta(delta_);

      if (debug()) 
		cerr << "DEBUG: Scenery::rayTrace(): i = " << i
		     << ", j = " << j << endl;
      
      data -> init(nbnuobs); // Initialize requested quantities to 0. or DBL_MAX

      if (impactcoords) {
	if (impactcoords[0] != DBL_MAX) {
	  ph_.setInitialCondition(gg_, obj_, impactcoords+8);
	  ph_.resetTransmission();
	  obj_ -> processHitQuantities(&ph_,
				       impactcoords+8,impactcoords,0.,data);
	}
      } else {
	screen_ -> getRayCoord(i,j, coord);
	ph_.setInitialCondition(gg_, obj_, coord);
	ph_.hit(data);
      }

      if (data) ++(*data); // operator++ is overloaded
      if (impactcoords) impactcoords += 16 ;

    }
  }
}

void Scenery::operator() (size_t i, size_t j,
			  Astrobj::Properties *data, double * impactcoords) {
  double coord[8];
  SmartPointer<Spectrometer> spr = screen_->getSpectrometer();
  ph_.setSpectrometer(spr);
  size_t nbnuobs = spr() ? spr -> getNSamples() : 0;
  ph_.setDelta(delta_);
  ph_.setTlim(tlim_);

  if (data) data -> init(nbnuobs);

  if (impactcoords) {
    if(impactcoords[0] != DBL_MAX) {
      ph_.setInitialCondition(gg_, obj_, impactcoords+8);
      ph_.resetTransmission();
      obj_ -> processHitQuantities(&ph_,impactcoords+8,impactcoords,0.,data);
    }
  } else {
    screen_ -> getRayCoord(i,j, coord);
    ph_.setInitialCondition(gg_, obj_, coord);
    ph_.hit(data);
  }
}

void Scenery::setRequestedQuantities(Gyoto::Quantity_t quant)
{quantities_=quant;}
void Scenery::setRequestedQuantities(std::string squant) {
  quantities_=0;
  char * tk = strtok(const_cast<char*>(squant.c_str()), " \t\n");
  while (tk != NULL) {
    if (!strcmp(tk, "Intensity"))
      quantities_ |= GYOTO_QUANTITY_INTENSITY;
    else if (!strcmp(tk, "EmissionTime"))
      quantities_ |= GYOTO_QUANTITY_EMISSIONTIME;
    else if (!strcmp(tk, "MinDistance"))
      quantities_ |= GYOTO_QUANTITY_MIN_DISTANCE;
    else if (!strcmp(tk, "FirstDistMin"))
      quantities_ |= GYOTO_QUANTITY_FIRST_DMIN;
    else if (!strcmp(tk, "Redshift"))
      quantities_ |= GYOTO_QUANTITY_REDSHIFT;
    else if (!strcmp(tk, "ImpactCoords"))
      quantities_ |= GYOTO_QUANTITY_IMPACTCOORDS;
    else if (!strcmp(tk, "Spectrum"))
      quantities_ |= GYOTO_QUANTITY_SPECTRUM;
    else if (!strcmp(tk, "BinSpectrum"))
      quantities_ |= GYOTO_QUANTITY_BINSPECTRUM;
    else if (!strcmp(tk, "User1"))
      quantities_ |= GYOTO_QUANTITY_USER1;
    else if (!strcmp(tk, "User2"))
      quantities_ |= GYOTO_QUANTITY_USER2;
    else if (!strcmp(tk, "User3"))
      quantities_ |= GYOTO_QUANTITY_USER3;
    else if (!strcmp(tk, "User4"))
      quantities_ |= GYOTO_QUANTITY_USER4;
    else if (!strcmp(tk, "User5"))
      quantities_ |= GYOTO_QUANTITY_USER5;
    else throwError("ScenerySubcontractor(): unkwon quantity"); 
    tk = strtok(NULL, " \t\n");
  }
  if (debug())
    cerr << "DEBUG: Scenery::setRequestedQuantities("<<squant<<"): "
	 << "quantities_=" << quantities_ << endl;
}
Gyoto::Quantity_t Scenery::getRequestedQuantities() const {
  return quantities_?quantities_:(obj_()?obj_->getDefaultQuantities():0);
}

std::string Scenery::getRequestedQuantitiesString() const {
  string squant = "";
  Quantity_t quantities
    = quantities_?quantities_:(obj_()?obj_->getDefaultQuantities():0);
  if (quantities & GYOTO_QUANTITY_INTENSITY   ) squant+="Intensity ";
  if (quantities & GYOTO_QUANTITY_EMISSIONTIME) squant+="EmissionTime ";
  if (quantities & GYOTO_QUANTITY_MIN_DISTANCE) squant+="MinDistance ";
  if (quantities & GYOTO_QUANTITY_FIRST_DMIN  ) squant+="FirstDistMin ";
  if (quantities & GYOTO_QUANTITY_REDSHIFT    ) squant+="Redshift ";
  if (quantities & GYOTO_QUANTITY_IMPACTCOORDS) squant+="ImpactCoords ";
  if (quantities & GYOTO_QUANTITY_SPECTRUM    ) squant+="Spectrum ";
  if (quantities & GYOTO_QUANTITY_BINSPECTRUM ) squant+="BinSpectrum ";
  if (quantities & GYOTO_QUANTITY_USER1       ) squant+="User1 ";
  if (quantities & GYOTO_QUANTITY_USER2       ) squant+="User2 ";
  if (quantities & GYOTO_QUANTITY_USER3       ) squant+="User3 ";
  if (quantities & GYOTO_QUANTITY_USER4       ) squant+="User4 ";
  if (quantities & GYOTO_QUANTITY_USER5       ) squant+="User5 ";
  return squant;
}

size_t Scenery::getScalarQuantitiesCount() const {
  size_t nquant=0;
  Quantity_t quantities
    = quantities_?quantities_:(obj_()?obj_->getDefaultQuantities():0);
  if (quantities & GYOTO_QUANTITY_INTENSITY   ) ++nquant;
  if (quantities & GYOTO_QUANTITY_EMISSIONTIME) ++nquant;
  if (quantities & GYOTO_QUANTITY_MIN_DISTANCE) ++nquant;
  if (quantities & GYOTO_QUANTITY_FIRST_DMIN  ) ++nquant;
  if (quantities & GYOTO_QUANTITY_REDSHIFT    ) ++nquant;
  //  SPECTRUM is not a SCALAR, don't add the following:
  //  if (quantities & GYOTO_QUANTITY_SPECTRUM    ) ++nquant;
  //  if (quantities & GYOTO_QUANTITY_BINSPECTRUM ) ++nquant;
  //  Idem IMPACTCOORDS:
  //  if (quantities & GYOTO_QUANTITY_IMPACTCOORDS) ++nquant;
  if (quantities & GYOTO_QUANTITY_USER1       ) ++nquant;
  if (quantities & GYOTO_QUANTITY_USER2       ) ++nquant;
  if (quantities & GYOTO_QUANTITY_USER3       ) ++nquant;
  if (quantities & GYOTO_QUANTITY_USER4       ) ++nquant;
  if (quantities & GYOTO_QUANTITY_USER5       ) ++nquant;
  return nquant;
}

double Scenery::getTlim() const { return tlim_; }
void Scenery::setTlim(double tlim) { tlim_ = tlim; }

#ifdef GYOTO_USE_XERCES
void Scenery::fillElement(FactoryMessenger *fmp) {
  if (gg_)     fmp -> setMetric (gg_) ;
  if (screen_) fmp -> setScreen (screen_) ;
  if (obj_)    fmp -> setAstrobj (obj_) ;
  if (delta_ != GYOTO_DEFAULT_DELTA)
    fmp -> setParameter ("Delta", delta_);
  if (getRequestedQuantities()) {
    fmp -> setParameter("Quantities", getRequestedQuantitiesString());
  }
  if (tlim_ != DEFAULT_TLIM) fmp -> setParameter("MinimumTime", tlim_);
}

SmartPointer<Scenery> Gyoto::ScenerySubcontractor(FactoryMessenger* fmp) {

  string name="", content="";
  double delta = GYOTO_DEFAULT_DELTA ;
  SmartPointer<Metric::Generic> gg = NULL;
  SmartPointer<Screen> scr = NULL;
  SmartPointer<Astrobj::Generic> ao = NULL;
  string squant = "";
  double tlim = DEFAULT_TLIM;

  gg = fmp->getMetric();
  scr= fmp->getScreen();
  ao = fmp->getAstrobj();


  while (fmp->getNextParameter(&name, &content)) {
    char* tc = const_cast<char*>(content.c_str());
    if (name=="Delta") delta = atof(tc);
    if (name=="Quantities") squant = content;
    if (name=="TLim") tlim = atof(tc);
  }

  SmartPointer<Scenery> sc = new Scenery(gg, scr, ao);
  sc -> setDelta(delta);
  sc -> setTlim(tlim);
  if (squant!="") sc -> setRequestedQuantities(squant);
  return sc;
}
#endif
