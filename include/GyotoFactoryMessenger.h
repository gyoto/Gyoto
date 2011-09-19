#ifdef GYOTO_USE_XERCES
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

#ifndef __GyotoFactoryMessenger_H_
#define __GyotoFactoryMessenger_H_


#ifndef XERCES_INCLUDE_WCHAR_H
#define XERCES_INCLUDE_WCHAR_H 0
#endif
#include <xercesc/dom/DOMElement.hpp>
#include <string>
#include <GyotoDefs.h>
#include <GyotoSmartPointer.h>

namespace Gyoto {
  class Factory;
  class FactoryMessenger;
  class Metric;
  class Astrobj;
  namespace Spectrum { class Generic ; }
  class Screen;
  class Photon;
}

class Gyoto::FactoryMessenger {
 private:
  Gyoto::Factory* employer_;
  xercesc::DOMElement *element_;
  xercesc::DOMNodeList* children_;
  XMLSize_t nodeCount_;
  XMLSize_t curNodeIndex_;
 public:
  FactoryMessenger(Gyoto::Factory*, xercesc::DOMElement*);
  FactoryMessenger(const FactoryMessenger& parent, std::string) ;

  void reset(); // get back to first parameter
  int getNextParameter(std::string* name, std::string* content);
  std::string getSelfAttribute(std::string attrname) const ;
  std::string getAttribute(std::string attrname) const ;
  std::string getFullContent() const ;
  FactoryMessenger * getChild() const ;

  /**
   * The child messenger is allocated with new, you need to delete it
   * after use.
   */
  FactoryMessenger* makeChild(std::string name);
  void setSelfAttribute(std::string attrname, std::string value) ;
  void setSelfAttribute(std::string attrname, unsigned long value) ;
  void setSelfAttribute(std::string attrname, unsigned int value) ;
  void setSelfAttribute(std::string attrname, double value) ;
  void setFullContent(std::string value) ; ///< Low level, prefer setParameter
  void setParameter(std::string name);
  void setParameter(std::string name, double value);
  void setParameter(std::string name, long int value);
  void setParameter(std::string name, unsigned int value);
  void setParameter(std::string name, unsigned long value);
  void setParameter(std::string name, int value);
  void setParameter(std::string name, std::string value);
  /**
   * If child is not NULL, a new FactoryMessenger is created to access
   * the new parameter element e.e. to set attributes in it. You then
   * need to delete the child.
   */
  void setParameter(std::string name, double val[], size_t n,
		    FactoryMessenger** child= NULL);
  void setAstrobj(SmartPointer<Astrobj>);
  void setScreen(SmartPointer<Screen>);
  void setMetric(SmartPointer<Metric>);
  SmartPointer<Metric>  getMetric  () ;
  SmartPointer<Screen>  getScreen  () ;
  SmartPointer<Photon>  getPhoton  () ;
  SmartPointer<Astrobj> getAstrobj () ;
  std::string fullPath(std::string) ;
};

#endif
#endif
