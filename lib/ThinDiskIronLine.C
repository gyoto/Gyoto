/*
  Copyright 2013 Frederic Vincent
  
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

#include "GyotoThinDiskIronLine.h"

#include <iostream>
#include <cmath>
#include <cfloat>
#include <cstdlib>

using namespace Gyoto;
using namespace Gyoto::Astrobj;
using namespace std;

Gyoto::Astrobj::ThinDiskIronLine::ThinDiskIronLine()
  : ThinDisk("ThinDiskIronLine"), plindex_(0.), linefreq_(0.), 
    cutradius_(-DBL_MAX)
{
  
  GYOTO_DEBUG << "Building ThinDiskIronLine" << endl;
}

Gyoto::Astrobj::ThinDiskIronLine::ThinDiskIronLine(const ThinDiskIronLine &o)
  : ThinDisk(o), plindex_(o.plindex_), linefreq_(o.linefreq_),
    cutradius_(o.cutradius_)
{
  GYOTO_DEBUG << "Copying ThinDiskIronLine" << endl;
}
ThinDiskIronLine * ThinDiskIronLine::clone() const { return new ThinDiskIronLine(*this); }

Gyoto::Astrobj::ThinDiskIronLine::~ThinDiskIronLine()
{
  GYOTO_DEBUG << "Destroying dummy ThinDiskIronLine" << endl;
}

double ThinDiskIronLine::emission(double nu_em, double dsem,
				  double *,
				  double coord_obj[8]) const{
  double rr=coord_obj[1];
  if (rr<cutradius_) return 0.;
  double dfreq=linefreq_/100.;
  /*
    NB: this choice of dfreq is related to the
    spectral resolution of e.g. CHANDRA, which is
    around E / (Delta E) = 100 at 6 keV,
    see e.g. iron line review of Reynolds 2003.
   */
  if (abs(nu_em-linefreq_)>dfreq) return 0.;
  
  double Iem = pow(rr,-plindex_);
  return Iem;
}

void ThinDiskIronLine::getVelocity(double const pos[4], double vel[4]) {
  if (pos[1]<cutradius_){
    //any velocity, emission=0 anyway
    for (int ii=1;ii<4;ii++) vel[ii]=0;
    vel[0] = gg_->SysPrimeToTdot(pos, vel+1);
  }else{
    ThinDisk::getVelocity(pos,vel);
  }
}


// to load from XML, we only need to implement this:
int ThinDiskIronLine::setParameter(std::string name,
			   std::string content,
			   std::string unit) {

  char* tc = const_cast<char*>(content.c_str());
  if (name=="PowerLawIndex") {
    plindex_=atof(tc);
  }
  else if (name=="LineFreq") {
    double freq=atof(tc);
    linefreq_=freq*1e3*1.60217657e-19/GYOTO_PLANCK;
  }
  else if (name=="CutRadius") {
    cutradius_=atof(tc);
  }
  else return ThinDisk::setParameter(name, content, unit);
  return 0;
}

#ifdef GYOTO_USE_XERCES
// to print/save XML (i.e. to prepare an XML from Yorick), we need this:
/*void ThinDiskIronLine::fillElement(FactoryMessenger *fmp) const {
  /* for instance:

     fmp -> setParameter("MyParameter", getMyParameter());

     Then call fillElement *on the direct parent*
  */
/*  ThinDisk::fillElement(fmp);
}*/
#endif


