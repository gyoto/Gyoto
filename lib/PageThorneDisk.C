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
#include "GyotoPageThorneDisk.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"


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
using namespace Gyoto::Astrobj;

PageThorneDisk::PageThorneDisk() :
  ThinDisk("PageThorneDisk"), aa_(0.), aa2_(0.),
  x0_(0.), x1_(0.), x2_(0.), x3_(0.)
{
  if (debug()) cerr << "DEBUG: PageThorneDisk Construction" << endl;
}

PageThorneDisk::PageThorneDisk(const PageThorneDisk& o) :
  ThinDisk(o), aa_(o.aa_), aa2_(o.aa2_),
  x0_(o.x0_), x1_(o.x1_), x2_(o.x2_), x3_(o.x3_)
{
  if (o.gg_()) gg_=o.gg_->clone();
  Generic::gg_=gg_;
}
PageThorneDisk* PageThorneDisk::clone() const
{ return new PageThorneDisk(*this); }

PageThorneDisk::~PageThorneDisk() {
  if (debug()) cerr << "DEBUG: PageThorneDisk Destruction" << endl;
}

void PageThorneDisk::updateSpin() {
  switch (gg_->getCoordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    aa_ = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getSpin();
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    aa_ = static_cast<SmartPointer<Metric::KerrKS> >(gg_) -> getSpin();
    break;
  default:
    throwError("PageThorneDisk::getSpin(): unknown COORDKIND");
  }
  aa2_=aa_*aa_;
  double z1 =1.+pow((1.-aa2_),1./3.)*(pow((1.+ aa_),1./3.)+pow((1.-aa_),1./3.));
  double z2 = pow(3.*aa2_ + z1*z1,0.5);
  double acosaao3= acos(aa_)/3.;

  x0_ = sqrt((3. + z2 - pow((3. - z1)*(3. + z1 + 2.*z2),0.5)));
  x1_ = 2.*cos(acosaao3 - M_PI/3.);
  x2_ = 2.*cos(acosaao3 + M_PI/3.); 
  x3_ = -2.*cos(acosaao3);

  rin_=(3.+z2-sqrt((3.-z1)*(3.+z1+2.*z2)));
}

void PageThorneDisk::setMetric(SmartPointer<Metric::Generic> gg) {
  string kind = gg->getKind();
  if (kind != "KerrBL" && kind != "KerrKS")
    throwError
      ("PageThorneDisk::setMetric(): metric must be KerrBL or KerrKS");
  ThinDisk::setMetric(gg);
  updateSpin();
}

double PageThorneDisk::emission(double, double dsem,
				    double *,
				    double coord_obj[8]) const{
  //See Page & Thorne 74 Eqs. 11b, 14, 15. This is F(r).

  double xx;
  switch (gg_->getCoordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    xx=sqrt(coord_obj[1]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    xx=pow(coord_obj[1]*coord_obj[1]+coord_obj[2]*coord_obj[2]-aa2_, 0.25);
    break;
  }
  double ff=
    3./(2.)*1./(xx*xx*(xx*xx*xx-3.*xx+2.*aa_))
    *( xx-x0_-3./2.*aa_*log(xx/x0_)
    -3.*(x1_-aa_)*(x1_-aa_)/(x1_*(x1_-x2_)*(x1_-x3_))*log((xx-x1_)/(x0_-x1_)) 
    -3.*(x2_-aa_)*(x2_-aa_)/(x2_*(x2_-x1_)*(x2_-x3_))*log((xx-x2_)/(x0_-x2_)) 
    -3.*(x3_-aa_)*(x3_-aa_)/(x3_*(x3_-x1_)*(x3_-x2_))*log((xx-x3_)/(x0_-x3_)));
           // f of Page&Thorne

  double Iem=ff/((4.*M_PI)*(xx*xx)); //with Mdot=1
  Iem*=M_PI;//assume isotropic emission --> flux is then equal to
  //"pi * specific intensity" (cf Ribicki Lightman)

  if (flag_radtransf_) Iem *= dsem;

  return Iem;

}

#ifdef GYOTO_USE_XERCES
void PageThorneDisk::fillElement(FactoryMessenger *fmp) const {
  fmp->setMetric(gg_);
  ThinDisk::fillElement(fmp);
}
#endif
