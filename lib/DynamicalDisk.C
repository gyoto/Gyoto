/*
    Copyright 2011, 2018 Frederic Vincent, Thibaut Paumard

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
#include "GyotoDynamicalDisk.h"
#include "GyotoProperty.h"
#include "GyotoUtils.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"
#include "GyotoKerrKS.h"

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

/// Properties

GYOTO_PROPERTY_START(DynamicalDisk)
GYOTO_PROPERTY_DOUBLE(DynamicalDisk, tinit, tinit)
GYOTO_PROPERTY_DOUBLE(DynamicalDisk, dt, dt)
GYOTO_PROPERTY_END(DynamicalDisk, PatternDiskBB::properties)

///


DynamicalDisk::DynamicalDisk() :
  PatternDiskBB(),
  dirname_(NULL),
  tinit_(0.),
  dt_(1.),
  nb_times_(0),
  nnu_(0), nphi_(0), nr_(0),
  emission_array_(NULL),
  velocity_array_(NULL),
  radius_array_(NULL)
{
  kind_="DynamicalDisk";
  GYOTO_DEBUG << "DynamicalDisk Construction" << endl;
  std::string str = "";
  dirname_ = new char[str.length() + 1];
  strcpy(dirname_, str.c_str());
}

DynamicalDisk::DynamicalDisk(const DynamicalDisk& o) :
  PatternDiskBB(o),
  dirname_(NULL),
  tinit_(o.tinit_),
  dt_(o.dt_),
  nb_times_(0),
  nnu_(o.nnu_), nphi_(o.nphi_), nr_(o.nr_),
  emission_array_(NULL),
  velocity_array_(NULL),
  radius_array_(NULL)
{
  GYOTO_DEBUG << "DynamicalDisk Copy" << endl;
#ifdef GYOTO_USE_CFITSIO
  if (o.dirname_) {
    dirname_ = new char[strlen(o.dirname_)+1];
    strcpy(dirname_,o.dirname_);
  }
  if (!nb_times_) return;
  emission_array_ = new double*[nb_times_] ;
  velocity_array_ = new double*[nb_times_] ;
  radius_array_   = new double*[nb_times_] ;
  for (int i=1; i<=nb_times_; i++) {
    size_t nel1=nnu_*nphi_*nr_, nel2=2*nr_*nphi_;
    emission_array_[i-1] = new double[nel1];
    velocity_array_[i-1] = new double[nel2];
    radius_array_  [i-1] = new double[nr_ ];
    memcpy(emission_array_[i-1], o.emission_array_[i-1], nel1*sizeof(double));
    memcpy(velocity_array_[i-1], o.velocity_array_[i-1], nel2*sizeof(double));
    memcpy(radius_array_  [i-1], o.radius_array_  [i-1], nr_ *sizeof(double));
  }
#endif
}
DynamicalDisk* DynamicalDisk::clone() const
{ return new DynamicalDisk(*this); }

DynamicalDisk::~DynamicalDisk() {
  GYOTO_DEBUG << "DynamicalDisk Destruction" << endl;
  for (int i=1; i<=nb_times_; i++) {
    if (emission_array_) delete [] emission_array_[i-1];
    if (velocity_array_) delete [] velocity_array_[i-1];
    if (radius_array_)   delete [] radius_array_[i-1];
  }
  if (emission_array_) delete [] emission_array_;
  if (velocity_array_) delete [] velocity_array_;
  if (radius_array_)   delete [] radius_array_  ;
  emission_array_ = NULL;
  velocity_array_ = NULL;
  radius_array_   = NULL;
  nb_times_ = 0;
  if (dirname_) delete dirname_;
}

double const * DynamicalDisk::getVelocity() const { return PatternDiskBB::getVelocity(); }

void DynamicalDisk::copyQuantities(int iq) {
  if (iq<1 || iq>nb_times_)
    GYOTO_ERROR("In DynamicalDisk::copyQuantities: incoherent value of iq");

  setEmission(emission_array_[iq-1]);
  setVelocity(velocity_array_[iq-1]);
  radius(radius_array_[iq-1]);
}

void DynamicalDisk::nullifyQuantities() {
  // This function is necessary to remove a malloc error
  // in PatternDisk destructor; I don't understand why (May 2017)
  setEmission(NULL);
  setVelocity(NULL);
  radius(NULL);
}

void DynamicalDisk::getVelocity(double const pos[4], double vel[4]) {
  double time = pos[0], tcomp=tinit_;
  int ifits=1;
  //cout << "time in dyna velo: " << time << " " << tinit_ << endl;
  while(time>tcomp && ifits<nb_times_){
    tcomp+=dt_;
    ifits++;
  }
  //cout << "ifits= " << ifits << endl;
  if (ifits==1 || ifits==nb_times_){
    copyQuantities(ifits);
    PatternDiskBB::getVelocity(pos,vel);
    //cout << "in dyna velo no interpo: " << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3] << endl;
  }else{
    double vel1[4], vel2[4];
    copyQuantities(ifits-1);
    //cerr << "in dyna go to get vel of pattern" << endl;
    PatternDiskBB::getVelocity(pos,vel1);
    //cerr << "in dyna after get vel of pattern" << endl;
    copyQuantities(ifits);
    PatternDiskBB::getVelocity(pos,vel2);
    double t1 = tinit_+(ifits-2)*dt_; // t1 <= time <= t1+dt_    
    for (int ii=0;ii<4;ii++){ // 1st order interpol
      vel[ii]=vel1[ii]+(vel2[ii]-vel1[ii])/dt_*(time-t1);
    }
    //cout << "in dyna velo dph/dtau_1, dph/dtau_2, dph/dtau, time, t1, dt= " << vel1[3] << " " << vel2[3] << " " << vel[3] << " " << time << " " << t1 << " " << dt_ << endl;
  }
  nullifyQuantities();
}

double DynamicalDisk::emission(double nu, double dsem,
			       state_t const &cp,
			       double const co[8]) const {
  GYOTO_DEBUG << endl;
  double time = co[0], tcomp=tinit_;
  //cout << "time in dyna emission: " << time << " " << tinit_ << endl;
  int ifits=1;
  while(time>tcomp && ifits<nb_times_){
    tcomp+=dt_;
    ifits++;
  }
  //cout << "ifits= " << ifits << endl;
  //if (ifits==1 || time>tinit_+(nb_times_-1)*dt_){
  if (ifits==1 || ifits==nb_times_){
    const_cast<DynamicalDisk*>(this)->copyQuantities(ifits); //awful trick to avoid problems with constness of function emission -> to improve
    double Iem=PatternDiskBB::emission(nu,dsem,cp,co);
    const_cast<DynamicalDisk*>(this)->nullifyQuantities();
    //cout << "in dyna emis no interpo: " << Iem << endl;
    return Iem;
  }else{
    double I1, I2;
    const_cast<DynamicalDisk*>(this)->copyQuantities(ifits-1);
    I1=PatternDiskBB::emission(nu,dsem,cp,co);
    const_cast<DynamicalDisk*>(this)->copyQuantities(ifits);
    I2=PatternDiskBB::emission(nu,dsem,cp,co);
    double t1 = tinit_+(ifits-2)*dt_;
    const_cast<DynamicalDisk*>(this)->nullifyQuantities();
    //cout << "in dyna I1, I2, time, t1, dt= " << I1 << " " << I2 << " " << time << " " << t1 << " " << dt_ << endl;
    //cout << "in dyna interpo time= " << I1+(I2-I1)/dt_*(time-t1) << endl;
    return I1+(I2-I1)/dt_*(time-t1);
  }
  return 0.;
}

std::string DynamicalDisk::file() const {return dirname_?dirname_:"";}
void DynamicalDisk::file(std::string const &fname) {
#ifdef GYOTO_USE_CFITSIO
    if (nb_times_) {
      // first free current arrays, if any
      for (int i=1; i<=nb_times_; i++) {
	if (emission_array_) delete [] emission_array_[i-1];
	if (velocity_array_) delete [] velocity_array_[i-1];
	if (radius_array_)   delete [] radius_array_  [i-1];
      }
      if (emission_array_) delete [] emission_array_;
      if (velocity_array_) delete [] velocity_array_;
      if (radius_array_)   delete [] radius_array_  ;
      emission_array_ = NULL;
      velocity_array_ = NULL;
      radius_array_   = NULL;
      nb_times_ = 0;
    }

    if (dirname_) delete dirname_;
    dirname_ = new char[strlen(fname.c_str())+1];
    strcpy(dirname_,fname.c_str());
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dirname_)) == NULL) {
      GYOTO_ERROR("In DynamicalDisk.C constructor : bad dirname_");
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
      GYOTO_ERROR("In DynamicalDisk.C: bad nb_times_ value");
    
    emission_array_ = new double*[nb_times_] ;
    velocity_array_ = new double*[nb_times_] ;
    radius_array_ = new double*[nb_times_] ;
    
    for (int i=1; i<=nb_times_; i++) {
      ostringstream stream_name ;
      stream_name << dirname_ << "data" 
		  << setw(4) << setfill('0') 
		  << i << ".fits.gz" ;
      
      string filename = stream_name.str();
      GYOTO_DEBUG << "Reading FITS file: " << filename << endl ;
      fitsRead(filename);
      size_t naxes[3];
      getIntensityNaxes(naxes);
      if (i==1){
	nnu_=naxes[0],nphi_=naxes[1],nr_=naxes[2];
      }else{
	if (nnu_!=naxes[0] || nphi_!=naxes[1] || nr_!=naxes[2]){
	  GYOTO_ERROR("In DynDisk::file: grid dimensions changing!");
	}
      }
      //	nel = (nnu=naxes[0])*(nphi=naxes[1])*(nr=naxes[2]);
      size_t nel1=nnu_*nphi_*nr_, nel2=2*nr_*nphi_;
      //save emission
      if (getIntensity()){
	double * emtemp = const_cast<double*>(getIntensity());
	emission_array_[i-1] = new double[nel1];
	for (size_t j=0;j<nel1;++j)
	  emission_array_[i-1][j]=emtemp[j];
      }else GYOTO_ERROR("In DynmicalDisk::file: Emission must be supplied");
      //save velocity
      if (getVelocity()){
	double * veltemp = const_cast<double*>(getVelocity());
	velocity_array_[i-1] = new double[nel2];
	for (size_t j=0;j<nel2;++j)
	  velocity_array_[i-1][j]=veltemp[j];
      }else GYOTO_ERROR("In DynmicalDisk::file: Velocity must be supplied");
      //save radius
      if (getGridRadius()){
      double * radtemp = const_cast<double*>(getGridRadius());
      radius_array_[i-1] = new double[nr_];
      for (size_t j=0;j<nr_;++j)
	radius_array_[i-1][j]=radtemp[j];
      }else GYOTO_ERROR("In DynmicalDisk::file: Radius must be supplied");
      //save other quantities
    }
#else
    GYOTO_ERROR("This Gyoto has no FITS i/o");
#endif
}

void DynamicalDisk::tinit(double t) {tinit_=t;}
double DynamicalDisk::tinit()const{return tinit_;}

void DynamicalDisk::dt(double t) {dt_=t;}
double DynamicalDisk::dt()const{return dt_;}

void DynamicalDisk::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const{
  if (p.name == "File"){
    fmp->setParameter("File", dirname_ );
  }
  else ThinDisk::fillProperty(fmp, p);
}