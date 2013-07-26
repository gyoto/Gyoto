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
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); throwError(ermsg); }

#include "GyotoPhoton.h"
#include "GyotoPatternDiskBB.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"


#include <fitsio.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <fstream>
#include <cstring>
#include <cmath>
#include <limits>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

PatternDiskBB::PatternDiskBB() :
  PatternDisk(),
  spectrumBB_(NULL),
  SpectralEmission_(0), PLDisk_(0),
  PLSlope_(0.), PLRho_(0.), rPL_(DBL_MAX)
{
  GYOTO_DEBUG << "PatternDiskBB Construction" << endl;
  spectrumBB_ = new Spectrum::BlackBody(); 
}

PatternDiskBB::PatternDiskBB(const PatternDiskBB& o) :
  PatternDisk(o),
  spectrumBB_(NULL),
  SpectralEmission_(o.SpectralEmission_), PLDisk_(o.PLDisk_),
  PLSlope_(o.PLSlope_),PLRho_(o.PLRho_),rPL_(o.rPL_)
{
  GYOTO_DEBUG << "PatternDiskBB Copy" << endl;
  if (o.spectrumBB_()) spectrumBB_=o.spectrumBB_->clone();
}
PatternDiskBB* PatternDiskBB::clone() const
{ return new PatternDiskBB(*this); }

PatternDiskBB::~PatternDiskBB() {
  GYOTO_DEBUG << "PatternDiskBB Destruction" << endl;
}

double const * PatternDiskBB::getVelocity() const { return PatternDisk::getVelocity(); }

void PatternDiskBB::getVelocity(double const pos[4], double vel[4]) {
  double rcur=projectedRadius(pos);
  double risco;
  switch (gg_->getCoordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    risco = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getRms();
    break;
  default:
    throwError("PatternDiskBB::getVelocity: bad COORDKIND");
    risco=0.;
  }
  
  double const * const radius=getGridRadius();
  size_t i[3]; // {i_nu, i_phi, i_r}
  //Search for indices only in non-power-law region
  if (rcur<rPL_)
    getIndices(i, pos, 0.); //NB: last arg should be nu, don't care here
  
  double rgrid=radius[i[2]];
  
  if ((getOuterRadius()==DBL_MAX && rcur>rPL_) || !getVelocity()){
    //Keplerian circ velocity for power law disk region
    //as well as if PatternDisk::velocity_ not provided
    ThinDisk::getVelocity(pos, vel);
  }else if (rgrid<risco){
    //default velocity, emission will be 0 there anyway
    vel[0]=1.;
    for (int ii=1;ii<4;ii++)
      vel[ii]=0.;
  }else{
    PatternDisk::getVelocity(pos, vel);
  }
}

double PatternDiskBB::emission(double nu, double dsem,
			       double *,
			       double co[8]) const{
  GYOTO_DEBUG << endl;

  double risco;
  switch (gg_->getCoordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    risco = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getRms();
    break;
  default:
    throwError("PatternDiskBB::emission: bad COORDKIND");
    risco=0.;
  }

  double rcur=projectedRadius(co);
  //  if (rcur > rout_ || rcur < risco) return 0.; // no emission in any case above rmax_
  size_t i[3]; // {i_nu, i_phi, i_r}

  //Search for indices only in non-power-law region
  if (rcur<rPL_)
    getIndices(i, co, nu);

  double const * const radius=getGridRadius();
  double rgrid=radius[i[2]];
  if (rgrid > rmax_ || rgrid < risco) return 0.; // no emission in any case above rmax_

  double Iem=0.;
  double const * const emiss = getIntensity();
  size_t naxes[3];
  getIntensityNaxes(naxes);
  size_t nnu=naxes[0], nphi=naxes[1];
  if (!SpectralEmission_){
    if (rPL_<DBL_MAX) 
      throwError("In PatternDisk.C: no power law region without SpectralEmission -> rPL_ should be DBL_MAX");
    Iem = emiss[i[2]*(nphi*nnu)+i[1]*nnu+i[0]];
  }else{ //Spectral emission    
    double TT;
    if (rcur<rPL_){
      // -> If r<rPL_ just read temperature value in emission_
      TT = emiss[i[2]*(nphi*nnu)+i[1]*nnu+i[0]];
      spectrumBB_->setTemperature(TT);
      Iem=(*spectrumBB_)(nu);
    }else if (PLDisk_){
      // -> If r>rPL_ compute temperature from first principles

      double rho_si = PLRho_*pow(rcur/risco,PLSlope_);
      //Assuming: pressure = kappa*(mass density)^gamma, gamma=5/3 (eq of state)
      // and: pressure = (mass density)*R/Mm*T (Mm = molar mass)

      double gamma=5./3.;
      double Mm=6e-4;//Navogadro*Munit/gamma
      double kappa=3e10;//pressure coef: p = kappa*rho^gamma, rho=mass density

      double cs2=kappa*gamma*pow(rho_si,gamma-1.);
      TT=Mm/GYOTO_GAS_CST*cs2;//Temperature in SI
      //cout << "TT after rl= " << TT << endl;
      //cout << "r,rho,T= " << rcross << " " << rho_si << " " << TT << endl;
      spectrumBB_->setTemperature(TT);
      Iem=(*spectrumBB_)(nu);
    }

    //cout << "Iem= " << Iem << endl;
  }

  if (!flag_radtransf_) return Iem;

  double thickness;
  double const * const opacity = getOpacity();
  if (rcur>rPL_)
    throwError("In PatternDiskBB::emission: optically thin integration not supported yet");
  if (opacity && (thickness=opacity[i[2]*(nphi*nnu)+i[1]*nnu+i[0]]*dsem))
    return Iem * (1. - exp (-thickness)) ;
  return 0.;
}

void PatternDiskBB::setMetric(SmartPointer<Metric::Generic> gg) {
  //Metric must be KerrBL (see emission function)
  string kind = gg->getKind();
  if (kind != "KerrBL")
    throwError
      ("PatternDiskBB::setMetric(): metric must be KerrBL");
  ThinDisk::setMetric(gg);
}

int PatternDiskBB::setParameter(std::string name,
				std::string content,
				std::string unit) {
  if      (name=="PLSlope"){
    PLDisk_=1;
    PLSlope_=atof(content.c_str());
    rPL_=getOuterRadius();//radius where power law disk begins
    setOuterRadius(DBL_MAX);//infinite power law disk
  }
  else if (name=="PLRho") PLRho_=atof(content.c_str());
  else if (name=="SpectralEmission") SpectralEmission_=1;
  else return PatternDisk::setParameter(name, content, unit);
  return 0;
}

#ifdef GYOTO_USE_XERCES
void PatternDiskBB::fillElement(FactoryMessenger *fmp) const {
  if (PLSlope_) fmp->setParameter("PLSlope", PLSlope_);
  fmp -> setParameter ( SpectralEmission_? "SpectralEmission" : "BolometricEmission");
  PatternDisk::fillElement(fmp);
}

#endif
