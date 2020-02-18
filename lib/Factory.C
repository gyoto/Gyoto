/*
    Copyright 2011-2016, 2018-2020 Thibaut Paumard

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

#include "GyotoConfig.h"
#ifdef GYOTO_USE_XERCES

#include "GyotoFactory.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoUtils.h"

#include <string>
#include <libgen.h>
#include <unistd.h>
#include <cstdlib>
#include <locale>

// Let's imbue 'C' locale to every stream to make sure decimal_point
// is always actually a '.'
static std::locale Cloc("C");

#include "GyotoMetric.h"
#include "GyotoAstrobj.h"
#include "GyotoSpectrum.h"
#include "GyotoSpectrometer.h"

#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/dom/DOMErrorHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/framework/MemBufFormatTarget.hpp>
#include <xercesc/framework/MemBufInputSource.hpp>

using namespace Gyoto;
using namespace xercesc;
using namespace std;

// max length should be at least 26: " -0.1234567890123456E+123 \0"
#define dvalLength 30
#define dfmt " %.16g "

// support DBL_MIN, DBL_MAX, and put right format
#define d2txt(txt, val)					\
  if      (val== DBL_MAX) strcpy(txt,  "DBL_MAX");	\
  else if (val==-DBL_MAX) strcpy(txt, "-DBL_MAX");	\
  else if (val== DBL_MIN) strcpy(txt,  "DBL_MIN");	\
  else if (val==-DBL_MIN) strcpy(txt, "-DBL_MIN");	\
  else {						\
    ostringstream ss;					\
    ss.imbue(Cloc);					\
    ss << setprecision(GYOTO_PREC)			\
       << setw(GYOTO_WIDTH) << val;			\
    strcpy(txt, ss.str().c_str());			\
  }

#define ifmt " %li "
#define i2txt(txt, val) sprintf( txt, ifmt, val);

class XStr
{
public :
    // -----------------------------------------------------------------------
    //  Constructors and Destructor
    // -----------------------------------------------------------------------
    XStr(const char* const toTranscode) : fASCIIForm(NULL)
    {
        // Call the private transcoding method
        fUnicodeForm = XMLString::transcode(toTranscode);
    }

    XStr(const XMLCh* const toTranscode) : fUnicodeForm(NULL)
    {
        // Call the private transcoding method
        fASCIIForm = XMLString::transcode(toTranscode);
    }

    ~XStr()
    {
      if (fUnicodeForm) XMLString::release(&fUnicodeForm);
      if (fASCIIForm)  XMLString::release(&fASCIIForm);
    }


    // -----------------------------------------------------------------------
    //  Getter methods
    // -----------------------------------------------------------------------
    const XMLCh* getXForm() const
    {
        return fUnicodeForm;
    }

    const char* getCForm() const
    {
        return fASCIIForm;
    }

private :
    // -----------------------------------------------------------------------
    //  Private data members
    //
    //  fUnicodeForm
    //      This is the Unicode XMLCh format of the string.
    // -----------------------------------------------------------------------
    XMLCh*   fUnicodeForm;
    char*    fASCIIForm;
};

#define X(str) XStr(str).getXForm()

string Cs (const XMLCh* str) { XStr bidule (str); string truc(bidule.getCForm()); return truc;}
//#define C(str) string(XStr(str).getCForm()).c_str()
#define C(str) Cs(str).c_str()

///// First setup the reporter class. It doesn't need to be in the API.
namespace Gyoto {
  class DOMErrorReporter;
}

class Gyoto::DOMErrorReporter : public ErrorHandler
{
public:
    // -----------------------------------------------------------------------
    //  Constructors and Destructor
    // -----------------------------------------------------------------------
    DOMErrorReporter() :
       fSawErrors(false)
    {
    }

    ~DOMErrorReporter()
    {
    }


    // -----------------------------------------------------------------------
    //  Implementation of the error handler interface
    // -----------------------------------------------------------------------
    void warning(const SAXParseException& toCatch);
    void error(const SAXParseException& toCatch);
    void fatalError(const SAXParseException& toCatch);
    void resetErrors();

    // -----------------------------------------------------------------------
    //  Getter methods
    // -----------------------------------------------------------------------
    bool getSawErrors() const;

    // -----------------------------------------------------------------------
    //  Private data members
    //
    //  fSawErrors
    //      This is set if we get any errors, and is queryable via a getter
    //      method. Its used by the main code to suppress output if there are
    //      errors.
    // -----------------------------------------------------------------------
    bool    fSawErrors;
};

void DOMErrorReporter::warning(const SAXParseException&)
{
  //y_print();
}

void DOMErrorReporter::error(const SAXParseException& toCatch)
{
    fSawErrors = true;
    GYOTO_ERROR(C(toCatch.getMessage()));
}

void DOMErrorReporter::fatalError(const SAXParseException& toCatch)
{
    fSawErrors = true;
    GYOTO_ERROR(C(toCatch.getMessage()));
}

void DOMErrorReporter::resetErrors()
{
    fSawErrors = false;
}

//// Now for the public API

Factory::Factory(char * data)
  : reporter_(NULL), gg_el_(NULL), obj_el_(NULL), ph_el_(NULL),
    scenery_(NULL), gg_(NULL), screen_(NULL), obj_(NULL), photon_(NULL),
    spectro_(NULL),
    filename_("")
{
  // Initialize Xerces XML parser

  XMLPlatformUtils::Initialize();
  parser_ = new XercesDOMParser();
  parser_->setValidationScheme(XercesDOMParser::Val_Never);
  parser_->setDoNamespaces(true);    // optional
  
  DOMErrorReporter *errReporter = new DOMErrorReporter();
  reporter_=errReporter;
  parser_->setErrorHandler(errReporter);
  // Parse file
  
  // If data start with "<?xml", this is an xml buffer.
  size_t len=strlen(data);
  if (len < 5 || data[0]!='<' || data[1] != '?'
      || data[2]!='x' || data[3]!='m' || data[4]!='l') {
    // we are dealing with a file name
    if (data[0]=='\\') ++data; // allow escaping 
    filename_=data;
    parser_->parse(data);
  } else {
    filename_="Gyoto XML data (in memory)";
    xercesc::MemBufInputSource
      xml_buf(reinterpret_cast<const XMLByte *>(data),
		len,
		filename_.c_str());
    parser_->parse(xml_buf);
  }
  
  doc_ = parser_ -> getDocument();
  root_ = doc_ -> getDocumentElement();
  if (!root_ ) throw(Error( "empty XML document" ));
  resolver_=doc_->createNSResolver(root_);

  kind_=(C(root_->getTagName()));

}

Factory::~Factory() {
  //  XMLString::release(&kind_);
  if (resolver_) delete resolver_;
  if (reporter_) delete reporter_;
  if (parser_) delete parser_;
  else if (doc_) delete doc_; // parser_ takes care of doc_
  // if (impl_) delete impl_; // Terminate takes care of that one
  XMLPlatformUtils::Terminate();
  gg_=NULL;
  obj_=NULL;
  scenery_=NULL;
  photon_=NULL;
  spectro_=NULL;
}

void Factory::setReporter(ErrorHandler* eh) {
  reporter_=eh;
  parser_->setErrorHandler(eh);
}

DOMElement * Factory::getRoot() { return root_; }
DOMDocument * Factory::getDoc() { return doc_; }

SmartPointer<Gyoto::Metric::Generic> Factory::metric() {
  if (!gg_) {
    DOMElement *MetricDOM;

    if (kind_.compare("Metric")){
      DOMXPathResult* result;
      result=doc_->evaluate(
			  X(("/"+kind_+"/Metric").c_str()),
			  root_,
			  resolver_,
			  DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE,
			  NULL);
      if (!result->getSnapshotLength()) {
	delete result;
	return NULL;
      }
      MetricDOM = static_cast< xercesc::DOMElement* >(result -> getNodeValue());
      delete result;
    } else MetricDOM = root_;

    vector<string> Plugin = Gyoto::split(C(MetricDOM->getAttribute(X("plugin"))), ",");
    string Kind =
      C(MetricDOM->getAttribute(X("kind")));
    FactoryMessenger fm(this, MetricDOM);

    gg_= (*Metric::getSubcontractor(Kind, Plugin))(&fm, Plugin);

  }

  return gg_;
}


SmartPointer<Gyoto::Astrobj::Generic> Factory::astrobj(){
    
  if (!obj_) {
    DOMXPathResult* result;
    DOMElement *tmpEl;

    if (kind_=="Astrobj") {
      obj_el_ = tmpEl = root_;
    } else {
      result=doc_->evaluate(
			  X(("/"+kind_+"/Astrobj").c_str()),
			  root_,
			  resolver_,
			  DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE,
			  NULL);
      if (!result->getSnapshotLength()) {
	delete result;
	return NULL;
      }
      tmpEl = static_cast< xercesc::DOMElement* >(result -> getNodeValue());
      delete result;
    }
    vector<string> Plugin = split(C(tmpEl->getAttribute(X("plugin"))), ",");
    string AstrobjKind =
      Cs(tmpEl->getAttribute(X("kind")));
    GYOTO_DEBUG_EXPR(AstrobjKind);

    FactoryMessenger fm(this, tmpEl);

    obj_ = (*Astrobj::getSubcontractor(AstrobjKind, Plugin))(&fm, Plugin);

  }
  return obj_;
}

SmartPointer<Gyoto::Photon> Factory::photon(){
    
  if (!photon_) {
    DOMXPathResult* result;
    DOMElement *tmpEl;

    if (kind_=="Photon") {
      ph_el_ = tmpEl = root_;
    } else {
      result=doc_->evaluate(
			  X(("/"+kind_+"/Photon").c_str()),
			  root_,
			  resolver_,
			  DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE,
			  NULL);
      if (!result->getSnapshotLength()) {
	delete result;
	return NULL;
      }
      tmpEl = static_cast< xercesc::DOMElement* >(result -> getNodeValue());
      delete result;
    }

    FactoryMessenger fm(this, tmpEl);
    photon_ = Photon::Subcontractor(&fm);

  }
  return photon_;
}

SmartPointer<Gyoto::Spectrum::Generic> Factory::spectrum(){
    
    DOMXPathResult* result;
    DOMElement *tmpEl;

    if (kind_=="Spectrum") {
      tmpEl = root_;
    } else {
      result=doc_->evaluate(
			  X(("/"+kind_+"/Spectrum").c_str()),
			  root_,
			  resolver_,
			  DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE,
			  NULL);
      if (!result->getSnapshotLength()) {
	delete result;
	return NULL;
      }
      tmpEl = static_cast< xercesc::DOMElement* >(result -> getNodeValue());
      delete result;
    }
    vector<string> Plugin = split(C(tmpEl->getAttribute(X("plugin"))), ",");
    string Kind =
      Cs(tmpEl->getAttribute(X("kind")));
    GYOTO_DEBUG_EXPR(Kind);

    FactoryMessenger fm(this, tmpEl);
    return (*Spectrum::getSubcontractor(Kind, Plugin))(&fm, Plugin);

}

SmartPointer<Gyoto::Spectrometer::Generic> Factory::spectrometer(){
    
    DOMXPathResult* result;
    DOMElement *tmpEl;

    if (kind_=="Spectrometer") {
      tmpEl = root_;
    } else {
      result=doc_->evaluate(
			  X(("/"+kind_+"/Spectrometer").c_str()),
			  root_,
			  resolver_,
			  DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE,
			  NULL);
      if (!result->getSnapshotLength()) {
	delete result;
	return NULL;
      }
      tmpEl = static_cast< xercesc::DOMElement* >(result -> getNodeValue());
      delete result;
    }
    vector<string> Plugin = split(C(tmpEl->getAttribute(X("plugin"))), ",");
    string Kind =
      Cs(tmpEl->getAttribute(X("kind")));
    GYOTO_DEBUG_EXPR(Kind);

    FactoryMessenger fm(this, tmpEl);
    return (*Spectrometer::getSubcontractor(Kind, Plugin))(&fm, Plugin);

}

/// Scenery

SmartPointer<Scenery> Factory::scenery () {
  if (!scenery_) {
    DOMXPathResult* result;
    DOMElement *tmpEl;

    result=doc_->evaluate(
			    X("/Scenery"),
			    root_,
			    resolver_,
			    DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE,
			    NULL);
    tmpEl = static_cast< xercesc::DOMElement* >(result -> getNodeValue());
  
    FactoryMessenger fm(this, tmpEl);
    scenery_ = Scenery::Subcontractor(&fm);
  
    delete result;
  }
  return scenery_;
}

SmartPointer<Gyoto::Screen> Factory::screen(){
  if (!screen_) {
    DOMXPathResult* result;
    DOMElement *ScreenDOM;
    result=doc_->evaluate(
			  X(("/"+kind_+"/Screen").c_str()),
			  root_,
			  resolver_,
			  DOMXPathResult::ORDERED_NODE_SNAPSHOT_TYPE,
			  NULL);
    if (!result->getSnapshotLength()) {
      delete result;
      return NULL;
    }
    
    ScreenDOM = static_cast< xercesc::DOMElement* >(result -> getNodeValue());
    
    FactoryMessenger fm ( this, ScreenDOM );
    screen_ = Screen::Subcontractor(&fm);
    delete result;
  }
  return screen_;
}

const string Factory::kind() { return kind_ ; }

Factory::Factory(SmartPointer<Scenery> sc)
  : reporter_(NULL), parser_(NULL), resolver_(NULL),
    gg_el_(NULL), obj_el_(NULL), ph_el_(NULL),
    scenery_(sc), gg_(sc->metric()),
    screen_(NULL), obj_(sc->astrobj()),
    photon_(NULL), spectro_(NULL), filename_("")
{
  GYOTO_DEBUG << "Initializing XML stuff" << endl;
  XMLPlatformUtils::Initialize();
  impl_ = DOMImplementationRegistry::getDOMImplementation(X("Core"));
  if (!impl_) GYOTO_ERROR("Problem initializing DOMImplementation");
  doc_ = impl_->createDocument(
                           0,                    // root element namespace URI.
                           X("Scenery"),         // root element name
                           0);                   // document type object (DTD).
  root_ = doc_->getDocumentElement();

  GYOTO_DEBUG << "Creating FactoryMessenger" << endl;
  FactoryMessenger fm(this, root_);

  GYOTO_DEBUG << "scenery_ -> fillElement(&fm);" << endl;
  scenery_ -> fillElement(&fm);

}

Factory::Factory(SmartPointer<Screen> scr)
  : reporter_(NULL), parser_(NULL), resolver_(NULL),
    gg_el_(NULL), obj_el_(NULL), ph_el_(NULL),
    scenery_(NULL), gg_(scr->metric()), screen_(scr), obj_(),
    photon_(NULL), spectro_(NULL), filename_("")
{
  XMLPlatformUtils::Initialize();
  impl_ = DOMImplementationRegistry::getDOMImplementation(X("Core"));
  if (!impl_) GYOTO_ERROR("Problem initializing DOMImplementation");
  doc_ = impl_->createDocument(
                           0,                    // root element namespace URI.
                           X("Screen"),         // root element name
                           0);                   // document type object (DTD).
  root_ = doc_->getDocumentElement();

  FactoryMessenger fm(this, root_);
  screen_ -> fillElement(&fm);  

}

Factory::Factory(SmartPointer<Metric::Generic> gg)
  : reporter_(NULL), parser_(NULL), resolver_(NULL),
    gg_el_(NULL), obj_el_(NULL), ph_el_(NULL),
    scenery_(NULL), gg_(gg), screen_(NULL), obj_(NULL),
    photon_(NULL), spectro_(NULL), filename_("")
{
  XMLPlatformUtils::Initialize();
  impl_ = DOMImplementationRegistry::getDOMImplementation(X("Core"));
  if (!impl_) GYOTO_ERROR("Problem initializing DOMImplementation");
  doc_ = impl_->createDocument(
                           0,                    // root element namespace URI.
                           X("Metric"),         // root element name
                           0);                   // document type object (DTD).
  gg_el_ = root_ = doc_->getDocumentElement();

  FactoryMessenger fm(this, gg_el_);
  gg -> fillElement(&fm);

}

Factory::Factory(SmartPointer<Astrobj::Generic> ao)
  : reporter_(NULL), parser_(NULL), resolver_(NULL), gg_el_(NULL),
    scenery_(NULL), gg_(NULL), obj_(ao), photon_(NULL),
    spectro_(NULL), filename_("")
{
  XMLPlatformUtils::Initialize();
  impl_ = DOMImplementationRegistry::getDOMImplementation(X("Core"));
  if (!impl_) GYOTO_ERROR("Problem initializing DOMImplementation");
  doc_ = impl_->createDocument(
                           0,                    // root element namespace URI.
                           X("Astrobj"),         // root element name
                           0);                   // document type object (DTD).
  obj_el_ = root_ = doc_->getDocumentElement();

  FactoryMessenger fm(this, obj_el_);
  ao -> fillElement(&fm);

}

Factory::Factory(SmartPointer<Spectrum::Generic> sp)
  : reporter_(NULL), parser_(NULL), resolver_(NULL), gg_el_(NULL),
    scenery_(NULL), gg_(NULL), obj_(NULL), photon_(NULL),
    spectro_(NULL), filename_("")
{
  XMLPlatformUtils::Initialize();
  impl_ = DOMImplementationRegistry::getDOMImplementation(X("Core"));
  if (!impl_) GYOTO_ERROR("Problem initializing DOMImplementation");
  doc_ = impl_->createDocument(
                           0,                    // root element namespace URI.
                           X("Spectrum"),         // root element name
                           0);                   // document type object (DTD).
  obj_el_ = root_ = doc_->getDocumentElement();

  FactoryMessenger fm(this, obj_el_);
  sp -> fillElement(&fm);

}

Factory::Factory(SmartPointer<Photon> ph)
  : reporter_(NULL), parser_(NULL), resolver_(NULL), 
    gg_el_(NULL), obj_el_(NULL), ph_el_(NULL),
    scenery_(NULL), gg_(ph->metric()), obj_(ph->astrobj()),
    photon_(ph), spectro_(NULL), filename_("")
{
  XMLPlatformUtils::Initialize();
  impl_ = DOMImplementationRegistry::getDOMImplementation(X("Core"));
  if (!impl_) GYOTO_ERROR("Problem initializing DOMImplementation");
  doc_ = impl_->createDocument(
                           0,                    // root element namespace URI.
                           X("Photon"),          // root element name
                           0);                   // document type object (DTD).
  ph_el_ = root_ = doc_->getDocumentElement();

  FactoryMessenger fm(this, root_);
  photon_ -> fillElement(&fm);

}

Factory::Factory(SmartPointer<Spectrometer::Generic> sp)
  : reporter_(NULL), parser_(NULL), resolver_(NULL), 
    gg_el_(NULL), obj_el_(NULL), ph_el_(NULL),
    scenery_(NULL), gg_(NULL), obj_(NULL),
    photon_(NULL), spectro_(sp), filename_("")
{
  XMLPlatformUtils::Initialize();
  impl_ = DOMImplementationRegistry::getDOMImplementation(X("Core"));
  if (!impl_) GYOTO_ERROR("Problem initializing DOMImplementation");
  doc_ = impl_->createDocument(
                           0,                    // root element namespace URI.
                           X("Spectrometer"),          // root element name
                           0);                   // document type object (DTD).
  root_ = doc_->getDocumentElement();

  FactoryMessenger fm(this, root_);
  spectro_ -> fillElement(&fm);

}

void Factory::metric(SmartPointer<Metric::Generic> gg, DOMElement *el) {
  if (!gg && !gg_) return;
  if (gg_ && gg && gg!= gg_) GYOTO_ERROR("Inconsistent use of Metrics");
  if (gg && !gg_el_) {
  
    gg_ = gg;

    gg_el_ = doc_->createElement(X("Metric"));
    el->appendChild(gg_el_);

    FactoryMessenger fm(this, gg_el_);
    gg -> fillElement(&fm);
  }

}

void Factory::astrobj(SmartPointer<Astrobj::Generic> ao, DOMElement *el) {
  GYOTO_DEBUG << endl;
  if (!ao && !obj_) return;
  if (obj_ && ao && ao!= obj_) GYOTO_ERROR("Inconsistent use of Astrobjs");
  if (ao && !obj_el_) {
    GYOTO_DEBUG <<"obj_ = ao;" << endl;
    obj_ = ao;

    GYOTO_DEBUG <<"XML stuff" << endl;
    obj_el_ = doc_->createElement(X("Astrobj"));
    el->appendChild(obj_el_);

    GYOTO_DEBUG <<"XML stuffnew FactoryMessenger" << endl;
    FactoryMessenger fm(this, obj_el_);

    GYOTO_DEBUG <<"ao -> fillElement(&fm);" << endl;
    ao -> fillElement(&fm);
  }

}

void Factory::screen(SmartPointer<Screen> scr, DOMElement *el) {
  if (!scr && !screen_) return;
  if (screen_ && scr && scr!= screen_)
    GYOTO_ERROR("Inconsistent use of Screens");
  
  if (scr && !screen_) {
    screen_ = scr;

    DOMElement * scr_el = doc_->createElement(X("Screen"));
    el->appendChild(scr_el);

    FactoryMessenger fm(this, scr_el);
    scr -> fillElement(&fm);
  }
}

void Factory::write(const char* const goutputfile) {
  filename_ = goutputfile;
  // write file
  DOMLSSerializer   *theSerializer
    = (static_cast<DOMImplementationLS*>(impl_))->createLSSerializer();
  DOMConfiguration  *serializerConfig
    = theSerializer->getDomConfig();
  DOMLSOutput       *theOutputDesc
    = (static_cast<DOMImplementationLS*>(impl_))->createLSOutput();
  XMLFormatTarget   *myFormTarget;

  if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint,true))
    serializerConfig->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
  if (goutputfile)
    myFormTarget=new LocalFileFormatTarget(goutputfile);
  else
    myFormTarget=new StdOutFormatTarget();
  
  theOutputDesc->setByteStream(myFormTarget);
  theSerializer->write(doc_, theOutputDesc);

  delete myFormTarget;
  theOutputDesc->release();
  theSerializer->release();
}

string Factory::format() {
  // write file
  DOMLSSerializer   *theSerializer
    = (static_cast<DOMImplementationLS*>(impl_))->createLSSerializer();
  DOMConfiguration  *serializerConfig
    = theSerializer->getDomConfig();
  DOMLSOutput       *theOutputDesc
    = (static_cast<DOMImplementationLS*>(impl_))->createLSOutput();
  MemBufFormatTarget   *myFormTarget = new MemBufFormatTarget();

  if (serializerConfig->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint,true))
    serializerConfig->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
  
  theOutputDesc->setByteStream(myFormTarget);
  theSerializer->write(doc_, theOutputDesc);

  string res=(const char*) myFormTarget->getRawBuffer();

  delete myFormTarget;
  theOutputDesc->release();
  theSerializer->release();

  return res;
}

void Factory::setParameter(std::string name, DOMElement *pel) {
  DOMElement*  el = doc_->createElement(X(name.c_str()));
  pel -> appendChild(el);
} 

void Factory::setParameter(std::string name, double value, DOMElement *pel) {
  DOMElement*  el = doc_->createElement(X(name.c_str()));
  pel -> appendChild(el);
  char val_string[dvalLength];
  d2txt(val_string,value);
  el->appendChild(doc_->createTextNode(X(val_string)));
} 

void Factory::setParameter(std::string name, int value, DOMElement *pel) {
  DOMElement*  el = doc_->createElement(X(name.c_str()));
  pel -> appendChild(el);
  char val_string[dvalLength];
  sprintf( val_string, " %i ", value);
  el->appendChild(doc_->createTextNode(X(val_string)));
} 

void Factory::setParameter(std::string name, unsigned int value, DOMElement *pel) {
  DOMElement*  el = doc_->createElement(X(name.c_str()));
  pel -> appendChild(el);
  char val_string[dvalLength];
  sprintf( val_string, " %u ", value);
  el->appendChild(doc_->createTextNode(X(val_string)));
} 

void Factory::setParameter(std::string name, long value, DOMElement *pel) {
  DOMElement*  el = doc_->createElement(X(name.c_str()));
  pel -> appendChild(el);
  char val_string[dvalLength];
  sprintf( val_string, " %li ", value);
  el->appendChild(doc_->createTextNode(X(val_string)));
} 

void Factory::setParameter(std::string name, unsigned long int value, DOMElement *pel) {
  DOMElement*  el = doc_->createElement(X(name.c_str()));
  pel -> appendChild(el);
  char val_string[dvalLength];
  sprintf( val_string, " %lu ", value);
  el->appendChild(doc_->createTextNode(X(val_string)));
} 

void Factory::setParameter(std::string name, std::string val, DOMElement *pel) {
  DOMElement*  el = doc_->createElement(X(name.c_str()));
  pel -> appendChild(el);
  el->appendChild(doc_->createTextNode(X(val.c_str())));
} 

void Factory::setParameter(std::string name, double val[],
			   size_t n, DOMElement *pel, FactoryMessenger **child){

  ostringstream ss;
  ss.imbue(Cloc); // set local to 'C'
  ss << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << val[0];
  for (size_t i=1; i<n; ++i) {
    ss << " " << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << val[i];
  }
  DOMElement*  el = doc_->createElement(X(name.c_str()));
  pel -> appendChild(el);
  el->appendChild( doc_->createTextNode(X(ss.str().c_str())) );
  if (child) *child = new FactoryMessenger(this, el);

}

void Factory::setParameter(std::string name,
			   std::vector<double> const &val,
			   DOMElement *pel, FactoryMessenger **child){

  DOMElement*  el = doc_->createElement(X(name.c_str()));
  pel -> appendChild(el);

  size_t n=val.size();

  if (n) {
    ostringstream ss;
    ss.imbue(Cloc); // set local to 'C'
    ss << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << val[0];
    for (size_t i=1; i<n; ++i) {
      ss << " " << setprecision(GYOTO_PREC) << setw(GYOTO_WIDTH) << val[i];
    }
    el->appendChild( doc_->createTextNode(X(ss.str().c_str())) );
  }
  if (child) *child = new FactoryMessenger(this, el);

}

void Factory::setParameter(std::string name,
			   std::vector<unsigned long> const &val,
			   DOMElement *pel, FactoryMessenger **child){

  DOMElement*  el = doc_->createElement(X(name.c_str()));
  pel -> appendChild(el);

  size_t n=val.size();

  if (n) {
    ostringstream ss;
    ss.imbue(Cloc); // set local to 'C'
    ss << val[0];
    for (size_t i=1; i<n; ++i) {
      ss << " " << val[i];
    }
    el->appendChild( doc_->createTextNode(X(ss.str().c_str())) );
  }
  if (child) *child = new FactoryMessenger(this, el);

}

void Factory::setContent(std::string content, DOMElement *el) {
  el -> appendChild( doc_->createTextNode( X( content.c_str() ) ) );
}

///////// Factory Messenger Class ///////////
FactoryMessenger::FactoryMessenger(Gyoto::Factory* emp,
				   xercesc::DOMElement* el) :
  employer_(emp), element_(el), curNodeIndex_(0)
{
  children_ = element_->getChildNodes();
  nodeCount_ = children_->getLength();
}

FactoryMessenger::FactoryMessenger(const FactoryMessenger& fm,
				   std::string name) :
  employer_(fm.employer_), element_(NULL), curNodeIndex_(0)
{
  element_ = employer_ -> doc_ -> createElement(X(name.c_str()));
  fm.element_ -> appendChild (element_) ;
  children_ = element_->getChildNodes();
  nodeCount_ = children_->getLength();
}

void FactoryMessenger::reset() {
  curNodeIndex_=0;
}

int FactoryMessenger::getNextParameter(std::string* namep,
				       std::string* contp,
				       std::string* unitp)
{


  GYOTO_DEBUG
    << "namep="    << namep
    << ", contp="  << contp << ": "
    << "*namep="   << *namep
    << ", *contp=" << *contp << endl;
    
  if (curNodeIndex_ >= nodeCount_) return 0;

  DOMNode *currentNode = children_->item(curNodeIndex_++);

  if( currentNode->getNodeType() &&  // true is not NULL
      currentNode->getNodeType() == DOMNode::ELEMENT_NODE ) // is element
    {
      // Found node which is an Element. Re-cast node as element
      DOMElement *currentElement
	= static_cast< xercesc::DOMElement* >( currentNode );
      *namep = C(currentElement->getTagName());
      *contp = C(currentElement->getTextContent());
      if (unitp) *unitp = C(currentElement->getAttribute(X("unit")));
      return 1;
    }
  return getNextParameter(namep, contp, unitp);
}

FactoryMessenger* FactoryMessenger::makeChild(std::string name) {
  return new FactoryMessenger(*this, name);
}

void FactoryMessenger::setSelfAttribute(std::string attrname,
					std::string attrvalue) {
  element_->setAttribute(X(attrname.c_str()), X(attrvalue.c_str()));
}

void FactoryMessenger::setSelfAttribute(std::string attrname,
					unsigned long attrvalue) {
  char val_string[dvalLength];
  sprintf( val_string, "%lu", attrvalue);
  element_->setAttribute(X(attrname.c_str()), X(val_string));
}

void FactoryMessenger::setSelfAttribute(std::string attrname,
					unsigned int attrvalue) {
  setSelfAttribute(attrname, (unsigned long)(attrvalue));
}

void FactoryMessenger::setSelfAttribute(std::string attrname,
					double attrvalue) {
  char val_string[dvalLength];
  d2txt(val_string, attrvalue);
  element_->setAttribute(X(attrname.c_str()), X(val_string));
}

string FactoryMessenger::getAttribute(std::string attrname) const {
  DOMNode *currentNode = children_->item(curNodeIndex_-1);
  DOMElement *currentElement
    = static_cast< xercesc::DOMElement* >( currentNode );
  return Cs(currentElement->getAttribute(X(attrname.c_str())));
}

string FactoryMessenger::getSelfAttribute(std::string attrname) const {
  return Cs(element_->getAttribute(X(attrname.c_str())));
}

string FactoryMessenger::getFullContent() const {
  return Cs(element_->getTextContent());
}

void FactoryMessenger::setFullContent(std::string content) {
  employer_ -> setContent (content, element_) ;
}

FactoryMessenger* FactoryMessenger::getChild() const {
  DOMNode *currentNode = children_->item(curNodeIndex_-1);
  DOMElement *currentElement
    = static_cast< xercesc::DOMElement* >( currentNode );
  return new FactoryMessenger(employer_, currentElement);
}

void FactoryMessenger::astrobj(SmartPointer<Astrobj::Generic> gg) {
  employer_ -> astrobj (gg, element_);
}

void FactoryMessenger::screen(SmartPointer<Screen> gg) {
  employer_ -> screen (gg, element_);
}

void FactoryMessenger::metric(SmartPointer<Metric::Generic> gg) {
  employer_ -> metric (gg, element_);
}

SmartPointer<Metric::Generic> FactoryMessenger::metric() {
  return employer_ -> metric ();
}

SmartPointer<Screen> FactoryMessenger::screen() {
  return employer_ -> screen ();
}

SmartPointer<Photon> FactoryMessenger::photon() {
  return employer_ -> photon ();
}

SmartPointer<Astrobj::Generic> FactoryMessenger::astrobj() {
  return employer_ -> astrobj ();
}

void FactoryMessenger::setParameter(std::string name){
  employer_ -> setParameter(name, element_);
}
void FactoryMessenger::setParameter(std::string name, double value){
  employer_ -> setParameter(name, value, element_);
}
void FactoryMessenger::setParameter(std::string name, int value){
  employer_ -> setParameter(name, value, element_);
}
void FactoryMessenger::setParameter(std::string name, long value){
  employer_ -> setParameter(name, value, element_);
}
void FactoryMessenger::setParameter(std::string name, unsigned int value){
  employer_ -> setParameter(name, value, element_);
}
void FactoryMessenger::setParameter(std::string name, unsigned long value){
  employer_ -> setParameter(name, value, element_);
}
void FactoryMessenger::setParameter(std::string name, std::string value){
  employer_ -> setParameter(name, value, element_);
}
void FactoryMessenger::setParameter(std::string name, double val[], size_t n,
				    FactoryMessenger **child){
  employer_ -> setParameter(name, val, n, element_, child);
}
void FactoryMessenger::setParameter(std::string name,
				    std::vector<double> const &val,
				    FactoryMessenger **child){
  employer_ -> setParameter(name, val, element_, child);
}

void FactoryMessenger::setParameter(std::string name,
				    std::vector<unsigned long> const &val,
				    FactoryMessenger **child){
  employer_ -> setParameter(name, val, element_, child);
}

size_t FactoryMessenger::parseArray(std::string content, double val[], size_t max_tokens)
{
  char const * const delim = " \t\n" ;
  char const * const c_content=content.c_str();
  size_t len=strlen(c_content);
  if (len==0) return 0;

  size_t n=0;
  char * const tc = new char[len+1];
  memcpy(tc, c_content, len+1);
  char * sub = strtok(tc, delim);

  while (sub && n<max_tokens) {
    val[n++] = Gyoto::atof(sub);
    sub = strtok(NULL, delim);
  }
  
  if (sub) n=max_tokens+1;

  delete [] tc;

  return n;
}

std::vector<double> FactoryMessenger::parseArray(std::string content)
{
  std::vector<double> result;
  char const * const delim = " \t\n" ;
  char const * const c_content=content.c_str();
  size_t len=strlen(c_content);
  if (len==0) return result;

  char * const tc = new char[len+1];
  memcpy(tc, c_content, len+1);
  char * sub = strtok(tc, delim);

  while (sub) {
    result.push_back(Gyoto::atof(sub));
    sub = strtok(NULL, delim);
  }
  
  delete [] tc;

  return result;
}

std::vector<unsigned long> FactoryMessenger::parseArrayULong(std::string content)
{
  std::vector<unsigned long> result;
  char const * const delim = " \t\n" ;
  char const * const c_content=content.c_str();
  size_t len=strlen(c_content);
  if (len==0) return result;

  char * const tc = new char[len+1];
  memcpy(tc, c_content, len+1);
  char * sub = strtok(tc, delim);

  while (sub) {
    result.push_back(atol(sub));
    sub = strtok(NULL, delim);
  }
  
  delete [] tc;

  return result;
}

std::string FactoryMessenger::fullPath(std::string fname) {
  return employer_ -> fullPath(fname);
}

std::string Factory::fullPath(std::string fname) {
  GYOTO_DEBUG << endl;
  if (!fname.compare(0, 1, "/")) return fname; // fname is already absolute 
  string fpath = "";

  string xmlpath="", curwd="";

  {
    // Make sure we call dirname and getcwd correctly and free memory.
    // We do this in a block to scope the temporary variables.

    char * xmlfile = strdup(filename_.c_str()); // why strdup? because
    xmlpath = dirname(xmlfile);	      // dirname may modify xmlfile.
    free (xmlfile); xmlfile = NULL;

    char * cwd = getcwd(NULL, 0);
    curwd = cwd;
    free (cwd); cwd = NULL;
  }

  string prefix="`pwd`/";
  if (fname.compare(0, prefix.size(), prefix)) {
    // fname does not start with "`pwd`/":
    // it is relative to xmlpath
    if (xmlpath.compare(0, 1, "/")) fpath = curwd + "/" ;
    fpath += xmlpath + "/";
    fpath += fname;
  } else {
    // fname starts with "`pwd`/": relative to working directory
    fpath = curwd + "/";
    fpath += fname.substr(prefix.size());
  }

  GYOTO_DEBUG << "returns " << fpath << endl;
  return fpath;
  
}

#endif
