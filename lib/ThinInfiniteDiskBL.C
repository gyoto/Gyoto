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
#include "GyotoThinInfiniteDiskBL.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <string>
#include <cmath>
#include <limits>
#include <string>

using namespace std;
using namespace Gyoto;

ThinInfiniteDiskBL::ThinInfiniteDiskBL(const SmartPointer<KerrBL>& metric) :
  Astrobj("ThinInfiniteDiskBL"), gg_(metric), Lr_(0.)
{
  if (debug()) cout << "ThinInfiniteDiskBL Construction" << endl;
  Astrobj::gg_=gg_;

  double aa=gg_->getSpin(), aa2=aa*aa;
  //ISCO radius, see Bardeen et al. 72, (2.21)
  double z1=1.+pow(1.-aa2,1./3.)*(pow(1.+aa,1./3.)+pow(1.-aa,1./3.)),
    z2= sqrt(3.*aa2+z1*z1);
  rmin_=(3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2)));
  
  
}

ThinInfiniteDiskBL::ThinInfiniteDiskBL(const ThinInfiniteDiskBL& o) :
  Astrobj(o),
  gg_(NULL), Lr_(o.Lr_), rmin_(o.rmin_)
{
  if (o.gg_()) gg_=o.gg_->clone();
  Astrobj::gg_=gg_;
}
ThinInfiniteDiskBL* ThinInfiniteDiskBL::clone() const
{ return new ThinInfiniteDiskBL(*this); }

ThinInfiniteDiskBL::~ThinInfiniteDiskBL() {
  if (debug()) cout << "ThinInfiniteDiskBL Destruction" << endl;
}

int ThinInfiniteDiskBL::Impact(Photon *ph, size_t index,
			       AstrobjProperties *data) {
  double coord_ph_hit[8], coord_obj_hit[8];
  double frac, rcross;
  double coord1[8], coord2[8];
  ph->getCoord(index, coord1);
  ph->getCoord(index+1, coord2);
  
  double theta1 = coord1[2];
  double theta2 = coord2[2];
  
  if (fabs(theta2-theta1) > M_PI)
    throwError ("ThinInfiniteDiskBL::Impact: fishy heuristic");
  
  theta1 -= M_PI*0.5;
  theta2 -= M_PI*0.5;
  while (theta1 < -M_PI) theta1 += 2.*M_PI;
  while (theta1 >= M_PI) theta1 -= 2.*M_PI;
  while (theta2 < -M_PI) theta2 += 2.*M_PI;
  while (theta2 >= M_PI) theta2 -= 2.*M_PI;
  
  if ( (theta1 > 0.) == (theta2 > 0.) && theta1 != 0. && theta2 != 0. )
    return 0;
  
  if ( theta2 == 0. ) frac = 1.;
  else frac = theta1 / (theta1-theta2) ; // == (0-th1) / (th2-th1)
  
  for (int i=0; i<8; ++i)
    coord_ph_hit[i] = coord1[i] + frac * (coord2[i] - coord1[i]);

  if ((rcross=coord_ph_hit[1]) < rmin_) return 0;

  //double coord_obj_hit[8], sp_em;
  //double sp_em;
  for (int i=0;i<4;i++) coord_obj_hit[i]=coord_ph_hit[i];
  double Omega=sqrt(1./(rcross*rcross*rcross));//angular Keplerian velocity
  
  double gtt=gg_->gmunu(coord_ph_hit,0,0), 
    gtph=gg_->gmunu(coord_ph_hit,0,3), 
    gphph=gg_->gmunu(coord_ph_hit,3,3);
  
  if (-1./(gtt+gphph*Omega*Omega+2.*gtph*Omega)<0.) {
    double rrcur= coord_ph_hit[1];
    double admis=1.1;
    if (rrcur < admis*rmin_){
      if (verbose() >= GYOTO_SEVERE_VERBOSITY) cout << "Thin Disk BL: Stopping at r= " << coord_ph_hit[1] 
						    << ", rmin disk= " << rmin_ 
						    << " : impossible to compute tdot here."
						    << " Intensity will be forced to 0." << endl;
      return 2;
    }else{
      cerr << "At r, r_isco, r/r_isco= " << rrcur << " " << rmin_ << " " << rrcur/rmin_ << endl;
      throwError("ThinInfiniteDiskBL.C: impossible condition in tdot determination");
    }
  }
  
  double tdot=sqrt(-1./(gtt+gphph*Omega*Omega+2.*gtph*Omega));
  coord_obj_hit[4]=tdot;coord_obj_hit[5]=coord_obj_hit[6]=0.;coord_obj_hit[7]=Omega*tdot;
  
  processHitQuantities(ph, coord_ph_hit, coord_obj_hit, 0., data);
  if (data) {
    if (data->rimpact) *data->rimpact=coord_obj_hit[1];
  }

  return 1;
}

double ThinInfiniteDiskBL::emission(double, double,
				    double *,
				    double coord_obj[8]) const{
  if (flag_radtransf_)
    throwError("Radiative transfer not implemented for DiskFromFile.");

  //See Page & Thorne 74 Eqs. 11b, 14, 15. This is F(r).
  double aa=gg_->getSpin(), aa2=aa*aa;//a* of Page & Thorne
  double  z1 = 1. + pow((1. - aa2),1./3.)*(pow((1. + aa),1./3.)
					   + pow((1. - aa),1./3.)); 
  double  z2 = pow(3.*aa2 + z1*z1,1./2.);
  double  x0 = sqrt((3. + z2 - pow((3. - z1)*(3. + z1 + 2.*z2),1./2.)));
  double x1 = 2.*cos(1./3.*acos(aa) - M_PI/3.),
    x2 = 2.*cos(1./3.*acos(aa) + M_PI/3.), 
    x3 = -2.*cos(1./3.*acos(aa));
  double xx=sqrt(coord_obj[1]);
  double ff=3./(2.)*1./(xx*xx*(xx*xx*xx-3.*xx+2.*aa))
    *( xx-x0-3./2.*aa*log(xx/x0)
       -3.*(x1-aa)*(x1-aa)/(x1*(x1-x2)*(x1-x3))*log((xx-x1)/(x0-x1)) 
       -3.*(x2-aa)*(x2-aa)/(x2*(x2-x1)*(x2-x3))*log((xx-x2)/(x0-x2)) 
       -3.*(x3-aa)*(x3-aa)/(x3*(x3-x1)*(x3-x2))*log((xx-x3)/(x0-x3)) );
           // f of Page&Thorne

  double Iem=1./(4.*M_PI)*1./(xx*xx)*ff; //with Mdot=1
  Iem*=M_PI;//assume isotropic emission --> flux is then equal to
  //"pi * specific intensity" (cf Ribicki Lightman)

  return Iem;

}

#ifdef GYOTO_USE_XERCES
void ThinInfiniteDiskBL::fillElement(FactoryMessenger *fmp) const {
  fmp->setMetric(gg_);
  Astrobj::fillElement(fmp);
}

SmartPointer<Astrobj> ThinInfiniteDiskBL::Subcontractor(FactoryMessenger* fmp) {
  string name, content;
  SmartPointer<ThinInfiniteDiskBL> ao =
    new ThinInfiniteDiskBL(fmp->getMetric());

  while (fmp->getNextParameter(&name, &content)) {
    ao -> setGenericParameter(name, content);
  }

  return ao;
}

void ThinInfiniteDiskBL::Init() {
  Astrobj::Register("ThinInfiniteDiskBL", &Subcontractor);
}
#endif
