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

#include "GyotoUtils.h"
#include <iostream>
#include <cstdlib>
#include <GyotoWorldlineIntegState.h>
#include <cmath>
#include <string>
#include <cstring>
#include <ctime>

using namespace std ; 
using namespace Gyoto;

WorldlineIntegState::WorldlineIntegState(){}//seems necessary for Kerr::IntegInit, I don't see why

WorldlineIntegState::WorldlineIntegState(SmartPointer<Metric::Generic> gg,
				   const double coord[8], const double delta) :
  gg_(gg),
  delta_(delta),
  deltainit_(delta)
{
  short i;
  for (i=0;i<8;++i) coord_[i]=coord[i];

  norm_=normref_= gg_->ScalarProd(coord,coord+4,coord+4);
}


int WorldlineIntegState::nextStep(Worldline* line, double coord[8], double del) {

  int j;
  int coutcoord=0;//to display coordinates updated values
  int coutnormdel=0;//to display new norm and integ step

  if (coutcoord){
    cout << "previous coord in Wl ada= " ;
    for (int ii=0;ii<8;ii++) cout << setprecision(10) << coord_[ii] << " ";
    cout << endl;
  }

  if (!del){
    if (gg_ -> myrk4_adaptive(line,coord_,norm_,normref_,coordnew_,delta_,h1_)) return 1;
    if (coutcoord){
      cout << "coordnew in Wl ada= " ;
      for (int ii=0;ii<8;ii++) cout << setprecision(10) << coordnew_[ii] << " ";
      cout << endl;
    }
    if (coutnormdel){
      cout << "norm apres RK ada= " << gg_ -> ScalarProd(coordnew_,coordnew_+4,coordnew_+4) << endl;
      cout << "pas apres RK ada= " << h1_ << endl;
    }

    delta_ = h1_;

    //TO USE TO OBTAIN NICE PLOTS OF ORBITS
    /*double deltamax = 100.; //NB: use this if integrating orbits in weak field (far from BH). Steps would be too big to allow "proper" plots (but still they would be correct!)
    
    if (delta_>deltamax) {
      if (debug()) cout << "NOTE: In WorldlineIntegState.C: delta_ becomes too big, reducing it to " << deltamax << endl;
      delta_=deltamax;
      }*/
  }else{
    delta_=del; // use requested step if non-zero

    if (gg_ -> myrk4(line,coord_,delta_,coordnew_)) return 1; 

    if (coutcoord){
      cout << "coordnew in Wl non ada= " ;
      for (int ii=0;ii<8;ii++) cout << coordnew_[ii] << " ";
      cout << endl;
    }
    if (coutnormdel){
      cout << "norm apres RK non ada= " << gg_ -> ScalarProd(coordnew_,coordnew_+4,coordnew_+4) << endl;
      cout << "pas apres RK non ada= " << delta_ << endl;
    }
  }

  for (j=0;j<8;j++) {
    coord_[j] = coordnew_[j];
    coord[j]  = coordnew_[j];
  }

  norm_=gg_ -> ScalarProd(coord,coord+4,coord+4);

  double normtol=.001;
  if (fabs(norm_-normref_)>normtol) {
    GYOTO_SEVERE << "***WARNING: in WlIntegState.C: norm is drifting"
      " - with norm,x1,x2,x3= " << norm_ << " " << coord[1] 
		 << " " << coord[2] << " " << coord[3] << " " << endl;
  }

  if (delta_==delta_+1) return 1; // delta == Infinity : stop condition

# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
  GYOTO_DEBUG_ARRAY(coord,8);
  GYOTO_DEBUG << "delta="<< del << ", delta_="<<delta_<<endl;
  GYOTO_ENDIF_DEBUG
# endif
  
  return 0;
}

double WorldlineIntegState::get_delta() {

  return delta_;
}

void WorldlineIntegState::set_delta(double delta) {
  delta_ = delta;
}

void WorldlineIntegState::setCoord(double coord[8]){
  for (int indice=0;indice<8;indice++){
    coord_[indice] = coord[indice];
  }
}

WorldlineIntegState::~WorldlineIntegState() {}
