#ifdef GYOTO_USE_XERCES

/**
 * \file GyotoFactory.h
 * \brief XML I/O
 *
 * Whenever implementing a new Astrobj or Metric, it's XML equivalent
 * should be implemented here.
 *
 */
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

#ifndef __GyotoFactory_H_
#define __GyotoFactory_H_

#ifndef XERCES_INCLUDE_WCHAR_H
#define XERCES_INCLUDE_WCHAR_H 0
#endif

#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <GyotoScenery.h>
#include <GyotoPhoton.h>
#include <GyotoSpectrum.h>
#include <sstream>
#include <string>

namespace Gyoto {
  class Factory;
  class FactoryMessenger;
  class Spectrometer;
}

class Gyoto::Factory
{
  friend class Gyoto::FactoryMessenger;

 protected:
  // XERCES MACHINERY
  xercesc::ErrorHandler *reporter_;
  xercesc::DOMDocument *doc_;
  xercesc::DOMElement *root_;
  xercesc::XercesDOMParser *parser_;
  xercesc::DOMXPathNSResolver* resolver_;
  xercesc::DOMImplementation* impl_;

  // Elements which must happen only once in a file
  // but may happen about anywhere
  xercesc::DOMElement *gg_el_;
  xercesc::DOMElement *obj_el_;
  xercesc::DOMElement *ph_el_;

  // GYOTO elements
  SmartPointer<Scenery> scenery_;
  SmartPointer<Metric::Generic> gg_;
  SmartPointer<Screen> screen_;
  SmartPointer<Astrobj> obj_;
  SmartPointer<Photon> photon_;
  SmartPointer<Spectrometer> spectro_;

  // Factory stuff
  std::string filename_;
  std::string kind_;

 public:
  // Constructor for READING
  Factory(char * filename);

  // Constructors for SAVING
  Factory(SmartPointer<Scenery> sc);
  Factory(SmartPointer<Metric::Generic> gg);
  Factory(SmartPointer<Astrobj> ao);
  Factory(SmartPointer<Spectrum::Generic> sp);
  Factory(SmartPointer<Screen> screen);
  Factory(SmartPointer<Photon> photon);
  Factory(SmartPointer<Spectrometer> Spectrometer);

  // Destructor
  ~Factory();

 private:
  // Internal stuff
  void setReporter(xercesc::ErrorHandler*);
  xercesc::DOMElement * getRoot();
  xercesc::DOMDocument* getDoc();

 public:
  // Name of TOP-LEVEL object (Scenery, Metric, Astrobj...)
  const std::string getKind();

  // Building and getting SmartPointer<OBJECTS>
  Gyoto::SmartPointer<Gyoto::Scenery> getScenery();
  Gyoto::SmartPointer<Gyoto::Metric::Generic>  getMetric();
  Gyoto::SmartPointer<Gyoto::Screen>  getScreen();
  Gyoto::SmartPointer<Gyoto::Astrobj> getAstrobj();
  Gyoto::SmartPointer<Gyoto::Photon>  getPhoton();
  Gyoto::SmartPointer<Gyoto::Spectrum::Generic>  getSpectrum();
  Gyoto::SmartPointer<Gyoto::AstrobjProperties> getAstrobjProperties();

  // XML OUTPUT
  void write(const char* const fname=0);
  std::string format();

  // Setting elements
  void setMetric(SmartPointer<Metric::Generic> gg, xercesc::DOMElement *el);
  void setAstrobj(SmartPointer<Astrobj> ao, xercesc::DOMElement *el);
  void setScreen(SmartPointer<Screen> scr, xercesc::DOMElement *el);
  void setContent(std::string content, xercesc::DOMElement *el);
  void setParameter(std::string name, xercesc::DOMElement *pel);
  void setParameter(std::string name, double value, xercesc::DOMElement *pel);
  void setParameter(std::string name, int value, xercesc::DOMElement *pel);
  void setParameter(std::string name, unsigned int value, xercesc::DOMElement *pel);
  void setParameter(std::string name, long value, xercesc::DOMElement *pel);
  void setParameter(std::string name, unsigned long value, xercesc::DOMElement *pel);
  void setParameter(std::string name, std::string sval, xercesc::DOMElement*);
  void setParameter(std::string , double val[], size_t, xercesc::DOMElement*,
		    FactoryMessenger **child = NULL);

  /**
   * Input is path relative to XML file. Output is absolute path to same file.
   */
  std::string fullPath(std::string relpath);
};

#endif
#endif
