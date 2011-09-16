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

#ifndef __GyotoRegister_H_
#define __GyotoRegister_H_


#ifndef XERCES_INCLUDE_WCHAR_H
#define XERCES_INCLUDE_WCHAR_H 0
#endif
#include <xercesc/dom/DOMElement.hpp>
#include <string>
#include <GyotoDefs.h>
#include <GyotoSmartPointer.h>

namespace Gyoto {
  namespace Register {
    class Entry;
    void init( char const * pluglist = NULL );
    void list();
  }
  class Factory;
  class factoryMessenger;
  class Metric;
  class Astrobj;
  namespace Spectrum { class Generic ; }
  class Screen;
  class Photon;
  void loadPlugin(   char const * const plugname, int nofail = 0);
}

class Gyoto::Register::Entry {
  friend void Register::list ();
protected:
  std::string name_;
  void* subcontractor_;
  int type_;
  Register::Entry* next_;
public:
  Entry(std::string name,
		void* subcontractor,
		Entry* next);
  ~Entry();
  void* getSubcontractor(std::string);
};

class Gyoto::factoryMessenger {
 private:
  Gyoto::Factory* employer_;
  xercesc::DOMElement *element_;
  xercesc::DOMNodeList* children_;
  XMLSize_t nodeCount_;
  XMLSize_t curNodeIndex_;
 public:
  factoryMessenger(Gyoto::Factory*, xercesc::DOMElement*);
  factoryMessenger(const factoryMessenger& parent, std::string) ;

  void reset(); // get back to first parameter
  int getNextParameter(std::string* name, std::string* content);
  std::string getSelfAttribute(std::string attrname) const ;
  std::string getAttribute(std::string attrname) const ;
  std::string getFullContent() const ;
  factoryMessenger * getChild() const ;

  /**
   * The child messenger is allocated with new, you need to delete it
   * after use.
   */
  factoryMessenger* makeChild(std::string name);
  void setSelfAttribute(std::string attrname, std::string value) ;
  void setSelfAttribute(std::string attrname, unsigned long value) ;
  void setSelfAttribute(std::string attrname, size_t value) ;
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
   * If child is not NULL, a new factoryMessenger is created to access
   * the new parameter element e.e. to set attributes in it. You then
   * need to delete the child.
   */
  void setParameter(std::string name, double val[], size_t n,
		    factoryMessenger** child= NULL);
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
