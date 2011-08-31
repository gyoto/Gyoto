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

#include "GyotoRegister.h"
#include "GyotoUtils.h"
#include "GyotoAstrobj.h"
#include "GyotoMetric.h"
#include "GyotoSpectrum.h"

#include <dlfcn.h>
#include <cstdlib>
#include <string>
#include <iostream>

using namespace Gyoto;
using namespace std;

typedef void GyotoInitFcn();

void Gyoto::loadPlugin(char const*const name, int nofail) {
  string dlfile = "libgyoto-" ;
  dlfile += name ;
  dlfile += "." ;
  dlfile += GYOTO_PLUGIN_SFX ;
  string dlfunc = "__Gyoto";
  dlfunc += name ;
  dlfunc += "Init";
  void* handle = NULL;
  GyotoInitFcn* initfcn = NULL;
  char * err = NULL;

  if (debug()) cerr << "DEBUG: loading plugin: " << name
		    << " from file: " << dlfile << endl;
  handle = dlopen(dlfile.c_str(), RTLD_LAZY | RTLD_GLOBAL);
  if (!handle && nofail) {
    if (verbose() >= GYOTO_DEFAULT_VERBOSITY)
      cerr << "WARNING: unable to load optional plugin " << dlfile << endl;
    return;
  }
  if ( (err=dlerror()) ) throwError(err);
  if (!handle) throwError((string("Failed to load pluging ")+dlfile).c_str());
  if (debug()) cerr << "DEBUG: calling plugin init function " << dlfunc << endl;
  initfcn = (GyotoInitFcn*)dlsym(handle, dlfunc.c_str());
  if ( (err=dlerror()) ) throwError(err);
  (*initfcn)();
}

void Gyoto::Register::init(char*  cpluglist) {

  // Clean registers
  Metric::initRegister();
  Astrobj::initRegister();
  Spectrum::initRegister();

  // Init built-in plug-ins
#ifdef GYOTO_BUILTIN_STDPLUG
  __GyotostdplugInit();
#endif
#ifdef GYOTO_BUILTIN_LORENEPLUG
  __GyotoloreneInit();
#endif

  // Load DL plug-ins

  if (!cpluglist) cpluglist = getenv("GYOTO_PLUGINS");
  if (!cpluglist) cpluglist = GYOTO_DEFAULT_PLUGINS;

  std::string pluglist = cpluglist;

  if (pluglist.length()) {
    size_t first=0, last=0;
    size_t len = pluglist.length() ;
    string curplug="";
    int nofail=0;
    while (pluglist.length()) {
      last=pluglist.find(",");
      nofail=0;
      curplug=pluglist.substr(0, last);
      if (debug())
	cerr << "DEBUG: first: " << first << ", last: " << last
	     << ", pluglist: |" << pluglist << "|"
	     << ", curplug: |" << curplug << "|" << endl;
      if (last <= pluglist.length()) pluglist=pluglist.substr(last+1);
      else pluglist="";
      if (!curplug.compare(0, 7, "nofail:")) {
	curplug = curplug.substr(7);
	nofail=1;
      }

      Gyoto::loadPlugin(curplug.c_str(), nofail);
      nofail=0;
    }
  }

  if (debug()) Register::list();

}

Register::Entry::Entry(std::string name,
			     void* subcontractor,
			     Register::Entry* next)
  : name_(name), subcontractor_((void*)subcontractor), next_(next)
{}

Register::Entry::~Entry() { if (next_) delete next_; }


void* Register::Entry::getSubcontractor(std::string name) {
  if (name_==name) return subcontractor_;
  if (next_) return next_ -> getSubcontractor(name);
  throwError ("Unregistered Metric kind: "+name);
}

void Gyoto::Register::list() {
  Register::Entry* entry = NULL;
  cout << "List of available Metrics:" << endl;
  entry = Metric::Register_;
  while (entry) {
    cout << "    " << entry -> name_ << endl;
    entry = entry -> next_ ;
  }
  
  cout << "List of available Astrobjs:" << endl;
  entry = Astrobj::Register_;
  while (entry) {
    cout << "    " << entry -> name_ << endl;
    entry = entry -> next_ ;
  }
  
  cout << "List of available Spectra:" << endl;
  for (entry = Spectrum::Register_; entry; entry = entry -> next_)
    cout << "    " << entry -> name_ << endl;
    
}
