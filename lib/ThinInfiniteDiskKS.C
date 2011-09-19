/*
    Copyright 2011 Frederic Vincent, Thibaut Paumard

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

#include "GyotoPhoton.h"
#include "GyotoThinInfiniteDiskKS.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cstring>
#include <cmath>
#include <limits>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

ThinInfiniteDiskKS::ThinInfiniteDiskKS(const SmartPointer<Metric::KerrKS>& metric) :
  Generic("ThinInfiniteDiskKS"), gg_(metric), Lr_(0.)
{
  if (debug()) cout << "ThinInfiniteDiskKS Construction" << endl;
  Generic::gg_=gg_;

  double aa=gg_->getSpin(), aa2=aa*aa;
  //ISCO radius, see Bardeen et al. 72, (2.21)
  double z1=1.+pow(1.-aa2,1./3.)*(pow(1.+aa,1./3.)+pow(1.-aa,1./3.)),
    z2= sqrt(3.*aa2+z1*z1);
  rmin_=(3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2)));
  
  
}

ThinInfiniteDiskKS::ThinInfiniteDiskKS(const ThinInfiniteDiskKS& o) :
  Generic(o),
  gg_(NULL), Lr_(o.Lr_), rmin_(o.rmin_)
{
  if (o.gg_()) gg_=o.gg_->clone();
  Generic::gg_=gg_;
}
ThinInfiniteDiskKS* ThinInfiniteDiskKS::clone() const
{ return new ThinInfiniteDiskKS(*this); }

ThinInfiniteDiskKS::~ThinInfiniteDiskKS() {
  if (debug()) cout << "ThinInfiniteDiskKS Destruction" << endl;
}

int ThinInfiniteDiskKS::Impact(Photon *ph, size_t index,
			       Astrobj::Properties *data) {
  double coord_ph_hit[8], coord_obj_hit[8];
  double frac, rcross;
  double coord1[8], coord2[8];
  ph->getCoord(index, coord1);
  ph->getCoord(index+1, coord2);

  double z1 = coord1[3], z2 = coord2[3];
  
  double aa=gg_->getSpin(), aa2=aa*aa;
  //double coord_ph_hit[8];

  if ( (z1 > 0) == (z2 > 0) && z1 != 0. && z2 != 0. ) return 0;
  if ( z2 == 0. ) frac = 1.;
  else frac = z1 / (z1-z2) ; // == (0-z1) / (z2-z1) 

  for (int i=0; i<8; ++i)
    coord_ph_hit[i] = coord1[i] + frac * (coord2[i] - coord1[i]);

  rcross = sqrt ( coord_ph_hit[1]*coord_ph_hit[1]
		  + coord_ph_hit[2]*coord_ph_hit[2] - aa2);
  if ( rcross < rmin_ ) return 0;
  
  for (int i=0;i<4;i++) coord_obj_hit[i]=coord_ph_hit[i];
  double Omega=sqrt(1./(rcross*rcross*rcross));//angular Keplerian velocity
  
  double SigmaBL=rcross*rcross; //theta=pi/2, cos(theta)=0, sin(theta=1)
  double gtt=-(1.-2.*rcross/SigmaBL), //BL values
    gtph=-2.*aa*rcross/SigmaBL,
    gphph=rcross*rcross+aa2+2.*aa2*rcross/SigmaBL;

  if (-1./(gtt+gphph*Omega*Omega+2.*gtph*Omega)<0.) throwError("ThinInfiniteDiskKS.C: impossible condition in tdot determination");
  double tdot=sqrt(-1./(gtt+gphph*Omega*Omega+2.*gtph*Omega));
  
  double xdot=-coord_ph_hit[2]*tdot*Omega, ydot=coord_ph_hit[1]*tdot*Omega;//KS xdot and ydot, see Hameury et al. 94
  coord_obj_hit[4]=tdot;coord_obj_hit[5]=xdot;coord_obj_hit[6]=ydot;coord_obj_hit[7]=0.;
  
  processHitQuantities(ph, coord_ph_hit, coord_obj_hit, 0., data);

  return 1;
}

double ThinInfiniteDiskKS::emission(double, double, double coord_ph[8],
				    double *) const{
  if (flag_radtransf_)
    throwError("Radiative transfer not implemented for ThinInfiniteDiskKS.");

  //See Page & Thorne 74 Eqs. 11b, 14, 15. This is F(r).
  double aa=gg_->getSpin(), aa2=aa*aa;//a* of Page & Thorne
  double  z1 = 1. + pow((1. - aa2),1./3.)*(pow((1. + aa),1./3.) 
					   + pow((1. - aa),1./3.)); 
  double  z2 = pow(3.*aa2 + z1*z1,1./2.);
  double  x0 = sqrt((3. + z2 - pow((3. - z1)*(3. + z1 + 2.*z2),1./2.)));
  double x1 = 2.*cos(1./3.*acos(aa) - M_PI/3.),
    x2 = 2.*cos(1./3.*acos(aa) + M_PI/3.),
    x3 = -2.*cos(1./3.*acos(aa));
  double rcross = sqrt ( coord_ph[1]*coord_ph[1]
			 + coord_ph[2]*coord_ph[2] - aa2); //z=0 at crossing
  double xx=sqrt(rcross);
  double ff=3./(2.)*1./(xx*xx*(xx*xx*xx-3.*xx+2.*aa))
    *( xx-x0-3./2.*aa*log(xx/x0)
       -3.*(x1-aa)*(x1-aa)/(x1*(x1-x2)*(x1-x3))*log((xx-x1)/(x0-x1)) 
       -3.*(x2-aa)*(x2-aa)/(x2*(x2-x1)*(x2-x3))*log((xx-x2)/(x0-x2)) 
       -3.*(x3-aa)*(x3-aa)/(x3*(x3-x1)*(x3-x2))*log((xx-x3)/(x0-x3)) );
                          // f of Page&Thorne

  if (ff!=ff) throwError("In ThinInfiniteDiskKS: emitted flux is nan!");
                        //Maybe an outgoing geodesic has crossed the horizon???

  double Iem=1./(4.*M_PI)*1./(xx*xx)*ff; //with Mdot=1

  Iem*=M_PI;//cf Marck 96; I take Mdot=1 and assume isotropic emission
	    //--> flux is then equal to "pi * specific intensity" (cf
	    //Ribicki Lightman)

  return Iem;
}

#ifdef GYOTO_USE_XERCES
void ThinInfiniteDiskKS::fillElement(FactoryMessenger *fmp) const {
  fmp->setMetric(gg_);
  Generic::fillElement(fmp);
}

SmartPointer<Astrobj::Generic> ThinInfiniteDiskKS::Subcontractor(FactoryMessenger* fmp) {
  string name, content;
  SmartPointer<ThinInfiniteDiskKS> ao =
    new ThinInfiniteDiskKS(fmp->getMetric());

  while (fmp->getNextParameter(&name, &content)) {
    ao -> setGenericParameter(name, content);
  }

  return ao;
}

void Gyoto::Astrobj::ThinInfiniteDiskKS::Init() {
  Astrobj::Register("ThinInfiniteDiskKS", &Subcontractor);
}
#endif
