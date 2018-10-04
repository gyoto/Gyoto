/**
 * \file GyotoFactoryMessenger.h
 * \brief Factory / SmartPointee::Subcontractor_t interface
 */
#ifdef GYOTO_USE_XERCES
/*
    Copyright 2011-2014, 2016, 2018 Thibaut Paumard

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

/// Internal to Xerces
/**
 * For obscure reasons, this needs to be set to 0 instead of just
 * undefined.
 */
#ifndef XERCES_INCLUDE_WCHAR_H
#define XERCES_INCLUDE_WCHAR_H 0
#endif
#include <xercesc/dom/DOMElement.hpp>
#include <vector>
#include <string>
#include <GyotoDefs.h>
#include <GyotoSmartPointer.h>

namespace Gyoto {
  class Factory;
  class FactoryMessenger;
  namespace Metric { class Generic; }
  namespace Astrobj { class Generic; }
  namespace Spectrum { class Generic ; }
  class Screen;
  class Photon;
}

/**
 * \class Gyoto::FactoryMessenger
 * \brief Factory / SmartPointee::Subcontractor_t interface
 *
 * A FactoryMessenger instance is like an employee passing messages
 * between its employer (the Factory) and a subcontractor (a function
 * of the SmartPointee::Subcontractor_t type).
 *
 * The FactoryMessenger also communicate with the fillElement method
 * of some classes (Astrobj::Generic::fillElement(),
 * Metric::Generic::fillElement(), Spectrum::Generic::fillElement()).
 *
 * A subcontractor function typically loops calling getNextParameter()
 * to read all the parameters provided for it in an XML file. If BASE
 * is one of Astrobj, Metric or Spectrum, and MyClass is an
 * implementation of BASE::Generic, the subcontractor static member
 * function often looks like this:
 *
 * \code
 *  SmartPointer<Gyoto::BASE::Generic> MyClass::Subcontractor(Gyoto::FactoryMessenger *messenger) {
 *  SmartPointer<Gyoto::BASE::MyClass> deliverable = new MyClass();
 *  while (messenger->getNextParameter(name, content) {
 *   if (name=="SomeProperty") deliverable -> setSomeProperty(content);
 *   else if (name=="AnotherProperty") deliverable -> setAnotherProperty(content);
 *  }
 *  return deliverable;
 * }
 * \endcode
 *
 * Other get* methods are provided to cope with more complex syntax
 * (e.g. when XML attributes are used, as in &lt;ParameterName
 * attribute="attrvalue"&gt;ParameterValue&lt;/ParameterName&gt;
 *
 * Conversely, the Factory asks certain Gyoto classes through their
 * fillElement() method how they should be printed or saved to an XML
 * file.  Those fillElement() methods use the FactoryMessenger::set*()
 * methods (in particular setParameter()) as well as, occasionally,
 * makeChild() to describe themselves to the Factory.
 */

class Gyoto::FactoryMessenger {
 private:
  Gyoto::Factory* employer_;
  ///< The Factory that sent this messenger
  xercesc::DOMElement *element_;
  ///< The XML element concerned by this transaction
  xercesc::DOMNodeList* children_;
  ///< The children of the XML element concerned by this transaction
  XMLSize_t nodeCount_;
  ///< The number of children of the XML element concerned by this transaction
  XMLSize_t curNodeIndex_;
  ///< Current child
 public:
  FactoryMessenger(Gyoto::Factory*, xercesc::DOMElement*);
  ///< Constructor called before subcontracting
  FactoryMessenger(const FactoryMessenger& parent, std::string) ;
  ///< Constructor called before fillElement

  void reset();
  ///< Get back to first parameter

  ///// GET METHODS, CALLED FROM A SUBCONTRACTOR

  /**
   * An Gyoto XML file may contain at most a single Metric section
   * and it may be present about anywhere in the XML tree. Individual
   * subcontractors should not try to interpret this section directly,
   * but should call metric() to find and interpret the Metric
   * section.
   */
  SmartPointer<Metric::Generic>  metric  () ;
  ///< Build and get the Metric described in this XML file 

  /**
   * An Gyoto XML file may contain at most a single Screen section
   * and it may be present about anywhere in the XML tree. Individual
   * subcontractors should not try to interpret this section directly,
   * but should call screen() to find and interpret the Screen
   * section.
   */
  SmartPointer<Screen>  screen  () ;
  ///< Build and get the Screen described in this XML file 

  /**
   * An Gyoto XML file may contain at most a single Photon section
   * and it may be present about anywhere in the XML tree. Individual
   * subcontractors should not try to interpret this section directly,
   * but should call photon() to find and interpret the Photon
   * section.
   */
  SmartPointer<Photon>  photon  () ;
  ///< Build and get the Photon described in this XML file 

  /**
   * An Gyoto XML file may contain at most a single Astrobj section
   * and it may be present about anywhere in the XML tree. Individual
   * subcontractors should not try to interpret this section directly,
   * but should call astrobj() to find and interpret the Astrobj
   * section.
   */
  SmartPointer<Astrobj::Generic> astrobj () ;
  ///< Build and get the Astrobj described in this XML file 

  /**
   * On each call, return a pair name-content of one of the
   * children_. Usually, "name" is the name of a parameter and
   * "content" is the string representation of the corresponding
   * value. For instance:
   * \code
   * <Name>Content</Name>
   * \endcode
   *
   * \param name upon output, name of the child
   * \param content of the child
   * \param unit= propertty of the child
   * \return 1 if there remains parameters to retrieve, 0 otherwise.
   */
  int getNextParameter(std::string* name,
		       std::string* content,
		       std::string* unit=NULL);
  ///< Get name and value of next parameter

  /**
   * For instance a Spectrometer description looks like this
   * \code
   * <Spectrometer kind="wave" nsamples="10"> 2.0e-6 2.4e-6</Astrobj>
   * \endcode
   * and the Spectrometer builder uses getSelfAttribute() to retrieve
   * the attributes "kind" and "nsamples".
   * 
   * \param attrname name of the attribute
   * \return attrvalue
   */
  std::string getSelfAttribute(std::string attrname) const ;
  ///< Get attribute of FactoryMessenger::element_

  /**
   * For instance
   * \code
   * <ParameterName attrname="attrvalue">ParameterContent</ParameterName>
   * \endcode
   * 
   * \param attrname name of the attribute
   * \return attrvalue
   */
  std::string getAttribute(std::string attrname) const ;
  ///< Get attribute of a last retrieved parameter

  /**
   * In exceptional circumstances, it may be necessary to get the
   * entire text content of the topmost element
   * FactoryMessenger::element_ instead or getting only the individual
   * FactoryMessenger::children_ .
   * 
   * For instance a Spectrometer description looks like this:
   * \code
   * <Spectrometer kind="wave" nsamples="10"> 2.0e-6 2.4e-6</Astrobj>
   * \endcode
   * and the Spectrometer builder uses getFullContent() to retrieve
   * the spectral boundaries (2.0e-6 and 2.4e-6 here).
   */
  std::string getFullContent() const ;
  ///< Get full content of element_

  /**
   * If one of the FactoryMessenger::children_ is complex (for
   * instance the complete description of a Gyoto::Spectrum), it is
   * possible to initialize a new FactoryMessenger and call the
   * correct subcontractor:
   * \code
   * SmartPointer<Spectrum::Generic> spectrum = NULL;
   * while (messenger->getNextParameter(name, content) {
   *  if (name=="Spectrum") {
   *   content = messenger->getAttribute("kind");
   *   FactoryMessenger* child = messenger->getChild();
   *   deliverable->spectrum( (*Spectrum::getSubcontractor(content))(child) );
   *   delete child;
   *  }
   * }
   * \endcode
   * The child is allocated with new and must be deleted after use.
   */
  FactoryMessenger * getChild() const ;
  ///< Get another FactoryMessenger instance initialized to current child 


  /**
   * This function takes a relative path (e.g. ../foo/bar.data) and
   * transforms it into a full path (starting with "/"). It is not
   * guaranteed to be portable (we assume that the path separator is
   * "/" and that absolute paths start with "/").
   *
   * \param relpath path relative to the directory where the XML file
   * on which the Factory works is located.
   *
   * \return fullpath at full path specification to the same point pon
   * the file-system.
   */
  std::string fullPath(std::string relpath) ;
  ///< Transform path into full path specification


  ///////// SET METHODS, CALLED FROM FILLELEMENT

  /**
   * At most one Metric section may be present in a give Gyoto XML file.
   *
   * When an object's fillElement() method is called, if this object
   * is connected to a Metric, it should call metric() with this
   * Metric. Very often, the Metric will already have been set
   * previously. The Factory will check that all the objects in the
   * hierarchy are attached to the same Metric instance, and save this
   * instance only once. Trying to set the Metric to something else
   * than the already set Metric instance is an error condition.
   *
   * To make things clearer: Assume "scenery" is a fully filled
   * Scenery. scenery->fillElement(messenger) will call:
   * \code
   * messenger->metric(Scenery::gg_)
   * messenger->screen(Scenery::screen_)
   * messenger->astrobj(Scenery::obj_);
   * \endcode
   *
   * The Factory will then call screen_->fillElement(child_messenger)
   * and obj_->fillElement(child_messenger), each of which will also
   * call metric(). If the same Metric is connected to the Astrobj,
   * to the Screen and to the Scenery, all is well. Else, you have a
   * bug to fix.
   */
  void metric(SmartPointer<Metric::Generic>);
  ///< Set the Metric

  /**
   * Same as metric(), but for the Astrobj.
   */
  void astrobj(SmartPointer<Astrobj::Generic>);
  ///< Set the Astrobj

  /**
   * Same as metric(), but for the Screen.
   */
  void screen(SmartPointer<Screen>);
  ///< Set the Screen


  /**
   * Create child XML element of the form
   * \code
   * <name/>
   * \endcode
   * for instance when "name" is boolean (present or absent), or only
   * takes attributes (see FactoryMessenger::setAttribute()). As an
   * example, Astrobj::Generic::fillElement() uses
   * setParameter() to set either Opticallythin or OpticallyThick.
   */
  void setParameter(std::string name);
  ///< Output parameter

  /**
   * Convert value to striing "svalue" and create an XML child element
   * of the form
   * \code
   * <name>svalue</name>
   * \endcode
   */
  void setParameter(std::string name, double value);
  ///< Output parameter

  /**
   * Convert value to striing "svalue" and create an XML child element
   * of the form
   * \code
   * <name>svalue</name>
   * \endcode
   */
  void setParameter(std::string name, long int value);
  ///< Output parameter

  /**
   * Convert value to striing "svalue" and create an XML child element
   * of the form
   * \code
   * <name>svalue</name>
   * \endcode
   */
  void setParameter(std::string name, unsigned int value);
  ///< Output parameter

  /**
   * Convert value to striing "svalue" and create an XML child element
   * of the form
   * \code
   * <name>svalue</name>
   * \endcode
   */
  void setParameter(std::string name, unsigned long value);
  ///< Output parameter

  /**
   * Convert value to string "svalue" and create an XML child element
   * of the form
   * \code
   * <name>svalue</name>
   * \endcode
   */
  void setParameter(std::string name, int value);
  ///< Output parameter

  /**
   * Create an XML child element of the form
   * \code
   * <name>value</name>
   * \endcode
   */
  void setParameter(std::string name, std::string value);
  ///< Output parameter

  /**
   * For instance:
   * \code
   * double val[4] = {1., 2., 3., 4.};
   * messenger->setParameter("MyArray", val, 4);
   * \endcode
   * will result in something like this:
   * \code
   * <MyArray>1.000000 2.000000 3.000000 4.000000</MyArray>
   * \endcode
   *
   * The exact format is unspecified, determined at compile time, and
   * by default, unlike in the example above, outputs a large number
   * of digits for each double (about 20).
   *
   * \param name the name of the parameter
   * \param val[] an array of doubles
   * \param n number of cells in val[]
   * \param child (optional) if not NULL, a new FactoryMessenger is
   * created to access the new parameter element e.g. to set
   * attributes in it (using setSelfAttribute()). You then need to
   * delete the child.
   */
  void setParameter(std::string name, double val[], size_t n,
		    FactoryMessenger** child= NULL);
  ///< Output an array of parameters

  void setParameter(std::string name, std::vector<double> const &val,
		    FactoryMessenger** child= NULL);
  void setParameter(std::string name, std::vector<unsigned long> const &val,
		    FactoryMessenger** child= NULL);
  ///< Output a vector of parameters

  /**
   * For instance Spectrometer::fillElement() sets its "kind"
   * attribute somewhat like this:
   * \code
   * messenger->setSelfAttribute("kind", "wave");
   * \endcode
   * to produce something like this:
   * \code
   * <Spectrometer kind="wave"/>
   * \endcode
   */
  void setSelfAttribute(std::string attrname, std::string value) ;
  ///< Set attribute in FactoryMessenger::element_

  /**
   * See setSelfAttribute(std::string attrname, std::string value)
   */
  void setSelfAttribute(std::string attrname, unsigned long value) ;
  ///< Set attribute in FactoryMessenger::element_

  /**
   * See setSelfAttribute(std::string attrname, std::string value)
   */
  void setSelfAttribute(std::string attrname, unsigned int value) ;
  ///< Set attribute in FactoryMessenger::element_

  /**
   * See setSelfAttribute(std::string attrname, std::string value)
   */
  void setSelfAttribute(std::string attrname, double value) ;
  ///< Set attribute in FactoryMessenger::element_

  /**
   * Exceptionnaly, a class instance may be best described by setting
   * the entire content of the corresponding element than by setting a
   * bunch of "parameters". This is the case of the spectrometer,
   * which sets a couple of attributes and reserves the full content
   * for the spectral boundaries (see Spectrometer::fillElement()).
   */
  void setFullContent(std::string value) ;
  ///< Low level, prefer setParameter()

  /**
   * To be used from fillElement() methods. For instance, the
   * Star::fillElement() method calls makeChild() to save the Star's
   * Spectrum and Opacity members somewhat like this:
   * \code
   * FactoryMessenger * child;
   * child = messenger-> makeChild("Spectrum");
   * spectrum_ -> fillElement(child);
   * delete child;
   * child = messenger-> makeChild("Opacity");
   * opacity_ -> fillElement(child);
   * delete child;
   * child=NULL;
   * \endcode
   *
   * The child messenger is allocated with new, you need to delete it
   * after use.
   */
  FactoryMessenger* makeChild(std::string name);
  ///< Create child FactoryMessenger

  /**
   * \brief Parse string into array
   *
   * Parse at most max_tokens tokens from string src into
   * pre-allocated array dst. Returns the number of tokens actually
   * found (interpreted using atof). dst must be at least of size
   * max_tokens.
   */
  static size_t parseArray(std::string src, double dst[], size_t max_tokens);

  /**
   * \brief Parse string into array
   *
   * Parse tokens from string src, returns them into a
   * std::vector<double>.
   */
  static std::vector<double> parseArray(std::string src);

  /**
   * \brief Parse string into array
   *
   * Parse tokens from string src, returns them into a
   * std::vector<unsigned long>.
   */
  static std::vector<unsigned long> parseArrayULong(std::string src);
};

#endif
#endif
