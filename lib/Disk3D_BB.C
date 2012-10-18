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
#include "GyotoDisk3D_BB.h"
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
#include <sstream>
#include <dirent.h>

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

Disk3D_BB::Disk3D_BB() :
  Disk3D(),
  spectrumBB_(NULL),
  tinit_(0.), dt_(1.)
{
  GYOTO_DEBUG << "Disk3D_BB Construction" << endl;
  spectrumBB_ = new Spectrum::BlackBody(); 
}

Disk3D_BB::Disk3D_BB(const Disk3D_BB& o) :
  Disk3D(o),
  spectrumBB_(NULL),
  tinit_(o.tinit_), dt_(o.dt_)
{
  GYOTO_DEBUG << "Disk3D_BB Copy" << endl;
  if (o.spectrumBB_()) spectrumBB_=o.spectrumBB_->clone();
}
Disk3D_BB* Disk3D_BB::clone() const
{ return new Disk3D_BB(*this); }

Disk3D_BB::~Disk3D_BB() {
  GYOTO_DEBUG << "Disk3D_BB Destruction" << endl;
  delete [] temperature_array_;
  delete [] velocity_array_;
}

double const * const Disk3D_BB::getVelocity() const { return Disk3D::getVelocity(); }

void Disk3D_BB::copyQuantities(int iq) {
  if (iq<1 || iq>nb_times_)
    throwError("In Disk3D_BB::copyQuantities: incoherent value of iq");
  double * curem = temperature_array_[iq-1],
    * curvel = velocity_array_[iq-1];

  setEmissquant(curem);
  setVelocity(curvel);
}

void Disk3D_BB::getVelocity(double const pos[4], double vel[4]) {
  double rcur=pos[1];
  double risco;
  switch (gg_->getCoordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    risco = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getRms();
    break;
  default:
    throwError("Disk3D_BB::getVelocity: bad COORDKIND");
  }

  if (rcur<risco){
    //default velocity, emission will be 0 there anyway
    vel[0]=1.;
    for (int ii=1;ii<4;ii++)
      vel[ii]=0.;
  }else{ 
    double time = pos[0], tcomp=tinit_;
    int ifits=1;
    while(time>tcomp && ifits<nb_times_){
      tcomp+=dt_;
      ifits++;
    }
    if (ifits==1 || ifits==nb_times_){
      copyQuantities(ifits);
      Disk3D::getVelocity(pos,vel);
    }else{
      double vel1[4], vel2[4];
      copyQuantities(ifits-1);
      Disk3D::getVelocity(pos,vel1);
      copyQuantities(ifits);
      Disk3D::getVelocity(pos,vel2);
      for (int ii=0;ii<4;ii++){ // 1st order interpol
	double t1 = tinit_+(ifits-2)*dt_;
	vel[ii]=vel1[ii]+(vel2[ii]-vel1[ii])/dt_*(time-t1);
      }
    }
  }
}

double Disk3D_BB::emission1date(double nu, double dsem,
			       double *,
			       double co[8]) const{
  GYOTO_DEBUG << endl;

  double * temperature = const_cast<double*>(getEmissquant());

  double risco;
  switch (gg_->getCoordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    risco = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getRms();
    break;
  default:
    throwError("Disk3D_BB::emission1date(): bad COORDKIND"
	       ", should be BL corrdinates");
  }

  double rcur=co[1];
  double th=co[2];

  if (rcur > rout() || rcur < risco) return 0.;

  size_t i[4]; // {i_nu, i_phi, i_z, i_r}
  getIndices(i,co,nu);
  size_t naxes[4];
  getEmissquantNaxes(naxes);
  size_t nnu=naxes[0], nphi=naxes[1], nz=naxes[2];
  double TT = temperature[i[3]*nphi*nz*nnu+i[2]*nphi*nnu+i[1]*nnu+i[0]];
  //This is local temperature in K

  spectrumBB_->setTemperature(TT);
  double Iem=(*spectrumBB_)(nu);

  double Ires=0.;
  if (!flag_radtransf_){
    Ires=Iem;
  }else{
    //SI value of cylindrical r coordinate:
    double dist_unit = GYOTO_G_OVER_C_SQUARE*gg_->getMass();
    double r_si=rcur*dist_unit; //spherical
    //double r_si=rcur*sin(th)*dist_unit; //cylindrical
    //cout << "cyl: " << rcur << " " << th << " " << rcur*sin(th) << " " << risco << endl;
    double risco_si=risco*dist_unit;

    //Local density in cgs:
    /*
      Following fact (in SI units) is defined by:
      fact=RR/(Mm*kappa*gamma)
      with:
       RR=8.3144621; // Perfect gas constant in SI
       Mm=6e-4; //This is N_avogadro*M_atomicmassunit/gamma in kg/mol
       kappa=3e10; // p = kappa*density^gamma
       gamma=1.66667;
      [see DynaDisk3D.i]
     */
    double fact=2.77149e-07;
    //See DynaDisk3D.i, or paper, for relation between TT and density
    double density=pow(fact*TT,1.5);// 1.5 is 1/(gamma-1)
    //density is in SI units, kg/m^3
    
    /*** Computing emission coef: ***/

    //Emission coef jnu for thermal bremsstrahlung (see RybickiLightman 5.14a)

    /* 
       fact2=1/4pi * 2^5*pi*e^6/(3*me*c^3) * sqrt(2pi/(3*kB*me)) * 1/mu^2 
       in SI, with: me=electron mass, e=electron charge, kB=boltzman,
                    mu=atomic mass unit
       Anyway this factor has no importance, 
       we are interested in relative values
    */
    double fact2=7.83315e-12;
    double hok=4.79924e-11; //planck cst / boltzman cst
    double jnu = fact2 * 1./sqrt(TT) * density*density
      * exp(-hok*nu/TT);

    //Elementary intensity added by current dsem segment of worldline
    //in SI units:
    Ires=jnu*dsem*dist_unit;
    //cout << "at end emission1date nu jnu deltaI= " << nu << " "  << jnu << " " << Ires << endl;

    //cout << rcur << " " << TT << endl;
    //cout << Iem << " " << Sem/Vem << " " << dsem << " " << Ires << endl;
    //cout << "stuff 3D= " << TT << " " << Iem << " " << Ires << " " << jnu << " " << Sem/Vem << " " << dsem << endl;
  }

  return Ires;

}

double Disk3D_BB::emission(double nu, double dsem,
			       double *,
			       double co[8]) const {
  GYOTO_DEBUG << endl;
  double time = co[0], tcomp=tinit_;
  int ifits=1;
  while(time>tcomp && ifits<nb_times_){
    tcomp+=dt_;
    ifits++;
  }
  double* fake;
  if (ifits==1 || ifits==nb_times_){
    const_cast<Disk3D_BB*>(this)->copyQuantities(ifits); //awful trick to avoid problems with constness of function emission -> to improve
    return emission1date(nu,dsem,fake,co);
  }else{
    double I1, I2;
    const_cast<Disk3D_BB*>(this)->copyQuantities(ifits-1);
    I1=emission1date(nu,dsem,fake,co);
    const_cast<Disk3D_BB*>(this)->copyQuantities(ifits);
    I2=emission1date(nu,dsem,fake,co);
    double t1 = tinit_+(ifits-2)*dt_;
    return I1+(I2-I1)/dt_*(time-t1);
  }

  return 0.;
}

double Disk3D_BB::transmission1date(double nu, double dsem,
			       double* fake,
			       double co[8]) const{
  GYOTO_DEBUG << endl;
  double * temperature = const_cast<double*>(getEmissquant());

  double dist_unit = GYOTO_G_OVER_C_SQUARE*gg_->getMass();
  
  double risco;
  switch (gg_->getCoordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    risco = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getRms();
    break;
  default:
    throwError("Disk3D_BB::emission1date(): bad COORDKIND"
	       ", should be BL corrdinates");
  }

  double rcur=co[1];
  double th=co[2];

  //cout << "rcur, risco= " << rcur << " " << risco << endl;

  if (rcur > rout() || rcur < risco) return 0.;

  size_t i[4]; // {i_nu, i_phi, i_z, i_r}
  getIndices(i,co,nu);
  size_t naxes[4];
  getEmissquantNaxes(naxes);
  size_t nnu=naxes[0], nphi=naxes[1], nz=naxes[2];
  double TT = temperature[i[3]*nphi*nz*nnu+i[2]*nphi*nnu+i[1]*nnu+i[0]];
  //This is local temperature in K
  
  spectrumBB_->setTemperature(TT);
  double BnuT=(*spectrumBB_)(nu); //Planck function
  double jnu=emission1date(nu,dsem,fake,co); // Emission coef
  double alphanu=0.; //absorption coef.
  if (BnuT==0.){
    /*
      BnuT can be 0 in the region close to ISCO where density
      decreases very fast. Then jnu should be 0 too (density~0).
      If both are 0, then nothing happens (no absorption, no emission,
      it's free space). Thus leave alphanu=0.
      If jnu!=0 then alphanu is not defined, this should not happen.
     */
    if (jnu!=0.)
      throwError("In Disk3D_BB::transmission1date absorption coef. undefined!");
  }else{
    alphanu=jnu/BnuT;
  }
  //cout << "at end transmission1date nu alphanu argexp= " << nu << " " << TT << " " << jnu << " " << BnuT << " " << exp(-alphanu*dsem*dist_unit) << endl;
  //Thermal emission assumed, use Kirchhoff alphanu=jnu/Bnu
  return exp(-alphanu*dsem*dist_unit); 
}

double Disk3D_BB::transmission(double nuem, double dsem, double* co) const {

  GYOTO_DEBUG << endl;
  double time = co[0], tcomp=tinit_;
  int ifits=1;
  while(time>tcomp && ifits<nb_times_){
    tcomp+=dt_;
    ifits++;
  }
  double* fake;
  if (ifits==1 || ifits==nb_times_){
    const_cast<Disk3D_BB*>(this)->copyQuantities(ifits); //awful trick to avoid problems with constness of function transmission -> to improve
    return transmission1date(nuem,dsem,fake,co);
  }else{
    double I1, I2;
    const_cast<Disk3D_BB*>(this)->copyQuantities(ifits-1);
    I1=transmission1date(nuem,dsem,fake,co);
    const_cast<Disk3D_BB*>(this)->copyQuantities(ifits);
    I2=transmission1date(nuem,dsem,fake,co);
    double t1 = tinit_+(ifits-2)*dt_;
    return I1+(I2-I1)/dt_*(time-t1);
  }

  return double(flag_radtransf_);
}


void Disk3D_BB::setMetric(SmartPointer<Metric::Generic> gg) {
  //Metric must be KerrBL (see emission function)
  string kind = gg->getKind();
  if (kind != "KerrBL")
    throwError
      ("Disk3D_BB::setMetric(): metric must be KerrBL");
  Disk3D::setMetric(gg);
}

int Disk3D_BB::setParameter(std::string name, std::string content) {
  if (name == "File") {
    dirname_ = new char[strlen(content.c_str())+1];
    strcpy(dirname_,content.c_str());
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dirname_)) == NULL) {
      throwError("In Disk3D_BB.C constructor : bad dirname_");
    }
    
    nb_times_=0;
    while ((dirp = readdir(dp)) != NULL) {
      nb_times_++;
    }
    nb_times_-=2; //for directories . and .. 
    
    /*
      NB: ***Caution***, here it is assumed that dirname_ 
      contains ONLY the FITS files, nothing else.
    */
    closedir(dp);
    
    GYOTO_DEBUG << "FITS directory, number of FITS files= " << 
      dirname_ << " " << nb_times_ << endl;
    
    if (nb_times_<1) 
      throwError("In Disk3D_BB.C: bad nb_times_ value");
    
    temperature_array_ = new double*[nb_times_] ;
    velocity_array_ = new double*[nb_times_] ;

    double nu0b, zminb, zmaxb, rinb, routb;
    size_t nnub, nphib, nzb, nrb;
    
    for (int i=1; i<=nb_times_; i++) {
      ostringstream stream_name ;
      stream_name << dirname_ << "pseudoN3D" 
		  << setw(4) << setfill('0') 
		  << i << ".fits.gz" ;
      
      string filename = stream_name.str();
      GYOTO_DEBUG << "Reading FITS file: " << filename << endl ;
      fitsRead(filename);
      size_t naxes[4];
      getEmissquantNaxes(naxes);
      size_t nnu=naxes[0], nphi=naxes[1], 
	nz=naxes[2], nr=naxes[3];
      size_t nel1=nnu*nphi*nz*nr, nel2=3*nr*nz*nphi;
      //save temperature
      if (getEmissquant()){
	double * emtemp = const_cast<double*>(getEmissquant());
	temperature_array_[i-1] = new double[nel1];
	for (size_t j=0;j<nel1;j++)
	  temperature_array_[i-1][j]=emtemp[j];
      }else throwError("In Disk3D_BB::setParameter: Temperature must be supplied");
      //save velocity
      if (getVelocity()){
	double * veltemp = const_cast<double*>(getVelocity());
	velocity_array_[i-1] = new double[nel2];
	for (size_t j=0;j<nel2;j++)
	  velocity_array_[i-1][j]=veltemp[j];
      }else throwError("In DynmicalDisk::setParameter: Velocity must be supplied");
      
      //check grid is constant
      if (i==1){
	nu0b=nu0();nnub=nnu;
	nphib=nphi;
	zminb=zmin();zmaxb=zmax();nzb=nz;
	rinb=rin();routb=rout();nrb=nr;
      }
      
      if (
	  nu0()!=nu0b || nnu!=nnub
	  || nphi!=nphib
	  || zmin()!=zminb || zmax()!=zmaxb || nz!=nzb
	  || rin()!=rinb || rout()!=routb || nr!=nrb
	  ) throwError("Disk3D_BB::setParameter Grid is not constant!");
    }
      
  }
  else if (name=="tinit") tinit_=atof(content.c_str());
  else if (name=="dt") dt_=atof(content.c_str());
  else return Disk3D::setParameter(name, content);
  return 0;
}
      
#ifdef GYOTO_USE_XERCES
void Disk3D_BB::fillElement(FactoryMessenger *fmp) const {
  if (tinit_) fmp->setParameter("tinit", tinit_);
  if (dt_) fmp->setParameter("dt", dt_);
  Disk3D::fillElement(fmp);
}

#endif
