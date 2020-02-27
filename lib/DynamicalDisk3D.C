/*
    Copyright 2013-2016, 2018-2020 Frederic Vincent, Thibaut Paumard

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
#include "GyotoDynamicalDisk3D.h"
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

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(DynamicalDisk3D)
GYOTO_PROPERTY_FILENAME(DynamicalDisk3D, File, file)
GYOTO_PROPERTY_DOUBLE(DynamicalDisk3D, tinit, tinit)
GYOTO_PROPERTY_DOUBLE(DynamicalDisk3D, dt, dt)
GYOTO_PROPERTY_BOOL(DynamicalDisk3D,
		    TemperatureGrid, IntensityGrid, temperature)
GYOTO_PROPERTY_DOUBLE(DynamicalDisk3D, PLindex, PLindex)
GYOTO_PROPERTY_BOOL(DynamicalDisk3D,
		    WithVelocity, NoVelocity, withVelocity)
GYOTO_PROPERTY_DOUBLE(DynamicalDisk3D, FloorTemperature, floorTemperature)
GYOTO_PROPERTY_END(DynamicalDisk3D, Disk3D::properties)

DynamicalDisk3D::DynamicalDisk3D() :
  Disk3D(),
  spectrumBB_(NULL),
  temperature_(1),
  dirname_(NULL),
  tinit_(0.),
  dt_(1.),
  nb_times_(1),
  PLindex_(3),
  novel_(0),
  floortemperature_(0),
  emission_array_(NULL),
  velocity_array_(NULL),
  absorption_array_(NULL)
{
  GYOTO_DEBUG << "DynamicalDisk3D Construction" << endl;
  spectrumBB_ = new Spectrum::BlackBody(); 
}

DynamicalDisk3D::DynamicalDisk3D(const DynamicalDisk3D& o) :
  Disk3D(o),
  spectrumBB_(NULL),
  temperature_(o.temperature_),
  dirname_(NULL),
  tinit_(o.tinit_),
  dt_(o.dt_),
  nb_times_(o.nb_times_),
  PLindex_(o.PLindex_),
  novel_(o.novel_),
  floortemperature_(o.floortemperature_),
  emission_array_(NULL),
  velocity_array_(NULL),
  absorption_array_(NULL)
{
  GYOTO_DEBUG << "DynamicalDisk3D Copy" << endl;
  if (o.spectrumBB_()) spectrumBB_=o.spectrumBB_->clone();

  if (o.dirname_){
    //dirname_ copy
    size_t length = strlen(o.dirname_)+1;
    dirname_ = new char[length];
    memcpy(dirname_, o.dirname_, length);
  }
  if (o.emission_array_ && o.velocity_array_){
    // emission_array_ and velocity_array_ copy
    size_t naxes[4];
    getEmissquantNaxes(naxes);
    size_t nnu=naxes[0], nphi=naxes[1], nz=naxes[2], nr=naxes[3];
    // Allocate
    emission_array_   = new double*[nb_times_];
    velocity_array_      = new double*[nb_times_];
    
    // Copy
    size_t nel1=nnu*nphi*nz*nr,
      nel2=3*nphi*nz*nr,
      szt=nel1*sizeof(double),
      szv=nel2*sizeof(double);
    for (int i=1; i<=nb_times_; i++) {
      emission_array_[i-1] = new double[nel1];
      velocity_array_[i-1] = new double[nel2];
      memcpy(emission_array_[i-1], o.emission_array_[i-1], szt);
      memcpy(velocity_array_[i-1], o.velocity_array_[i-1], szv);
    }

    // If absorption is given, copy
    if (o.absorption_array_){
      absorption_array_   = new double*[nb_times_];
      for (int i=1; i<=nb_times_; i++) {
	absorption_array_[i-1] = new double[nel1];
	memcpy(absorption_array_[i-1], o.absorption_array_[i-1], szt);
      }
    }
  }
  
  
}
DynamicalDisk3D* DynamicalDisk3D::clone() const
{ return new DynamicalDisk3D(*this); }

bool DynamicalDisk3D::isThreadSafe() const {
  // spectrumBB_ is not handled by a property.
  return Disk3D::isThreadSafe()
    && (!spectrumBB_ || spectrumBB_ -> isThreadSafe());
}

DynamicalDisk3D::~DynamicalDisk3D() {
  GYOTO_DEBUG << "DynamicalDisk3D Destruction" << endl;
  if (emission_array_) delete [] emission_array_;
  if (absorption_array_) delete [] absorption_array_;
  if (velocity_array_) delete [] velocity_array_;
}

double const * DynamicalDisk3D::getVelocity() const { return Disk3D::getVelocity(); }

void DynamicalDisk3D::copyQuantities(int iq) {
  if (iq<1 || iq>nb_times_)
    GYOTO_ERROR("In DynamicalDisk3D::copyQuantities: incoherent value of iq");
  setEmissquant(emission_array_[iq-1]);
  if (absorption_array_) opacity(absorption_array_[iq-1]);
  setVelocity(velocity_array_[iq-1]);
}

void DynamicalDisk3D::getVelocity(double const pos[4], double vel[4]) {
  if (novel_){
    // Velocity of emitted particle is not provided (only bulk velocity
    // is known). Then put velocity to default, redshift factor will
    // be constant.
    vel[0]=1.;vel[1]=0.;vel[2]=0.;vel[3]=0.;
  }else{
    double rcur=pos[1];
    double risco;
    switch (gg_->coordKind()) {
    case GYOTO_COORDKIND_SPHERICAL:
      {
      string kin = gg_->kind();
      if (kin == "KerrBL")
	risco = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getRms();
      else if (kin == "Minkowski")
	risco = 6.;
      else
	risco=0., GYOTO_ERROR("In DynamicalDisk3D::getVelocity: bad metric");
      break;
      }
    default:
      GYOTO_ERROR("DynamicalDisk3D::getVelocity: bad COORDKIND");
      risco=0.;
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
}

double DynamicalDisk3D::emission1date(double nu, double dsem,
			       state_t const &,
			       double const co[8]) const{
  GYOTO_DEBUG << endl;
  
  double * emiss = const_cast<double*>(getEmissquant());

  double risco;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    {
      string kin = gg_->kind();
      if (kin == "KerrBL")
	risco = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getRms();
      else if (kin == "Minkowski")
	risco = 6.;
      else
	risco=0., GYOTO_ERROR("In DynamicalDisk3D::getVelocity: bad metric");
    break;
    }
  default:
    GYOTO_ERROR("DynamicalDisk3D::emission1date(): bad COORDKIND"
	       ", should be BL corrdinates");
    risco=0.;
  }

  double rcur=co[1];

  if (rcur*fabs(sin(co[2])) > rout() || rcur < risco) return 0.;

  size_t i[4]; // {i_nu, i_phi, i_z, i_r}
  getIndices(i,co,nu);
  size_t naxes[4];
  getEmissquantNaxes(naxes);
  size_t nnu=naxes[0], nphi=naxes[1], nz=naxes[2];
  double emissq = emiss[i[3]*nphi*nz*nnu+i[2]*nphi*nnu+i[1]*nnu+i[0]];
  //This is a quantity from which emiss. coef. j_nu is known
  //this is temperature in K if temperature_=1
  //or directly j_nu, up to the nu dependence, if temperature_=0

  double Ires=-1.; // emitted specific intensity

  if (!flag_radtransf_){ // optically thick case
    
    if (temperature_){
      spectrumBB_->temperature(emissq);
      Ires=(*spectrumBB_)(nu);
      //cout << "in emis: " << emissq << " " << Ires << endl;
    }else{
      Ires=emissq;
    }

  }else{//                 // optically thin case

    if (temperature_){
      // cout << "in DD3 Ires= " << rcur << " " << th  << " " << emissq << " " << floortemperature_ << " ";
      if (emissq<floortemperature_){
	//cout << "return 0 " << emissq << " " << floortemperature_ << endl;
	//throwwError("test dynad3d");
	Ires=0.;
      }else{
	// BB radiation
	spectrumBB_->temperature(emissq);
	Ires=(*spectrumBB_)(nu);
	//cout << "return  " << emissq << " " << Ires << endl;
	// BELOW: BREMS computation for 2012 RWI paper
	// //SI value of cylindrical r coordinate:
	// double dist_unit = gg_->unitLength();
	// /*
	//   Following fact (in SI units) is defined by:
	//   fact=RR/(Mm*kappa*gamma)
	//   with:
	//   RR=8.3144621; // Perfect gas constant in SI
	//   Mm=6e-4; //This is N_avogadro*M_atomicmassunit/gamma in kg/mol
	//   kappa=3e10; // p = kappa*density^gamma
	//   gamma=1.66667;
	//   [see DynaDisk3D.i]
	// */
	// double fact=2.77149e-07;
	// //See DynaDisk3D.i, or paper, for relation between emissq and density
	// double density=pow(fact*emissq,1.5);// 1.5 is 1/(gamma-1)
	// //density is in SI units, kg/m^3
	
	// /*** Computing emission coef: ***/
	
	// //Emission coef jnu for thermal bremsstrahlung 
	// // (see RybickiLightman 5.14a)
	
	// /* 
	//    fact2=1/4pi * 2^5*pi*e^6/(3*me*c^3) * sqrt(2pi/(3*kB*me)) * 1/mu^2 
	//    in SI, with: me=electron mass, e=electron charge, kB=boltzman,
	//    mu=atomic mass unit
	//    Anyway this factor has no importance, 
	//    we are interested in relative values
	// */
	// double fact2=7.83315e-12;
	// double hok=4.79924e-11; //planck cst / boltzman cst
	// double jnu = fact2 * 1./sqrt(emissq) * density*density
	//   * exp(-hok*nu/emissq);
	
	// //Elementary intensity added by current dsem segment of worldline
	// //in SI units:
	
	// Ires=jnu*dsem*dist_unit; // usd e.g. for 3D RWI computation
	// //cout << "stuff: " << density << " " << fact2 << " " << emissq << " " << dsem << " " << jnu << endl;
      }
      // cout << Ires << endl; // TEST
    }else{
      //Ires=Iem; // used e.g. for GC blob computation
      double dist_unit = gg_->unitLength()*100.; // unit length in cgs
      double jnu = emissq*pow(nu,-(PLindex_-1.)/2.);
      Ires=jnu*dsem*dist_unit;
    }

  }

  return Ires;

}

double DynamicalDisk3D::emission(double nu, double dsem,
			       state_t const &cph,
			       double const co[8]) const {
  GYOTO_DEBUG << endl;
  double time = co[0], tcomp=tinit_;
  int ifits=1;
  while(time>tcomp && ifits<nb_times_){
    tcomp+=dt_;
    ifits++;
  }

  if (ifits==1 || ifits==nb_times_){
    const_cast<DynamicalDisk3D*>(this)->copyQuantities(ifits); //awful trick to avoid problems with constness of function emission -> to improve
    return emission1date(nu,dsem,cph,co);
  }else{
    double I1, I2;
    const_cast<DynamicalDisk3D*>(this)->copyQuantities(ifits-1);
    I1=emission1date(nu,dsem,cph,co);
    const_cast<DynamicalDisk3D*>(this)->copyQuantities(ifits);
    I2=emission1date(nu,dsem,cph,co);
    double t1 = tinit_+(ifits-2)*dt_;
    return I1+(I2-I1)/dt_*(time-t1);
  }

  return 0.;
}

double DynamicalDisk3D::transmission1date(double nu, double dsem,
			       state_t const &,
			       double const co[8]) const{
  GYOTO_DEBUG << endl;
  if (!flag_radtransf_) return 0.;
  
  double risco;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    {
      string kin = gg_->kind();
      if (kin == "KerrBL")
	risco = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> getRms();
      else if (kin == "Minkowski")
	risco = 6.;
      else
	risco=0., GYOTO_ERROR("In DynamicalDisk3D::getVelocity: bad metric");
    break;
    }
  default:
    GYOTO_ERROR("DynamicalDisk3D::emission1date(): bad COORDKIND"
	       ", should be BL corrdinates");
    risco=0.;
  }

  double rcur=co[1];
 
  if (rcur*fabs(sin(co[2])) > rout() || rcur < risco) return 0.;

  size_t i[4]; // {i_nu, i_phi, i_z, i_r}
  getIndices(i,co,nu);
  size_t naxes[4];
  getEmissquantNaxes(naxes);
  size_t nnu=naxes[0], nphi=naxes[1], nz=naxes[2];

  if (temperature_){
    double * emiss = const_cast<double*>(getEmissquant());
    double emissq = emiss[i[3]*nphi*nz*nnu+i[2]*nphi*nnu+i[1]*nnu+i[0]];
    //emissq is local temperature in K

    //cout << "in trans DD3: " << rcur << " " << th << " " << emissq << " " << floortemperature_ << endl;
    if (emissq<floortemperature_) return 1.;
    else return 0.;
    
    // BELOW: absorption for RWI 2012 paper
    // spectrumBB_->temperature(emissq);
    // double BnuT=(*spectrumBB_)(nu); //Planck function
    // double jnu=emission1date(nu,dsem,NULL,co); // Emission coef * ds
    // double alphanu=0.; //absorption coef.
    // if (BnuT==0.){
    //   /*
    // 	BnuT can be 0 in the region close to ISCO where density
    // 	decreases very fast. Then jnu should be 0 too (density~0).
    // 	If both are 0, then nothing happens (no absorption, no emission,
    // 	it's free space). Thus leave alphanu=0.
    // 	If jnu!=0 then alphanu is not defined, this should not happen.
    //   */
    //   if (jnu!=0.){
    // 	cout << "r= " << rcur << " " << emissq << " " << jnu << " " << BnuT << endl;
    // 	GYOTO_ERROR("In DynamicalDisk3D::"
    // 		   "transmission1date absorption coef. undefined!");
    //   }
    // }else{
    //   alphanu=jnu/BnuT;
    // }
    // //Thermal emission assumed, use Kirchhoff alphanu=jnu/Bnu
    // return exp(-alphanu); // the dsem factor is already included
    //                       //in alphanu via jnu=emission1date(...,dsem,...)
  }else{
    if (absorption_array_){
      double * abs = const_cast<double*>(opacity());
      double absq = abs[i[3]*nphi*nz*nnu+i[2]*nphi*nnu+i[1]*nnu+i[0]];
      double dist_unit = gg_->unitLength()*100.; //dist unit in cgs
      double alphanu=absq*pow(nu,-(PLindex_+4.)/2.);
      return exp(-alphanu*dsem*dist_unit);
    }else{
      GYOTO_ERROR("In DynamicalDisk3D: in non-BB optically thin case, "
		 "opacity should be provided");
    }
  }
  GYOTO_ERROR("BUG: should not reach this point!");
  return 0.; // avoid pedantic warning
}

double DynamicalDisk3D::transmission(double nuem, double dsem, state_t const &cp, double const *co) const {

  GYOTO_DEBUG << endl;
  double time = co[0], tcomp=tinit_;
  int ifits=1;
  while(time>tcomp && ifits<nb_times_){
    tcomp+=dt_;
    ifits++;
  }

  if (ifits==1 || ifits==nb_times_){
    const_cast<DynamicalDisk3D*>(this)->copyQuantities(ifits); //awful trick to avoid problems with constness of function transmission -> to improve
    return transmission1date(nuem,dsem,cp,co);
  }else{
    double I1, I2;
    const_cast<DynamicalDisk3D*>(this)->copyQuantities(ifits-1);
    I1=transmission1date(nuem,dsem,cp,co);
    const_cast<DynamicalDisk3D*>(this)->copyQuantities(ifits);
    I2=transmission1date(nuem,dsem,cp,co);
    double t1 = tinit_+(ifits-2)*dt_;
    return I1+(I2-I1)/dt_*(time-t1);
  }

  return double(flag_radtransf_);
}


void DynamicalDisk3D::metric(SmartPointer<Metric::Generic> gg) {
  //Metric must be KerrBL (see emission function)
  string kin = gg->kind();
  if (kin != "KerrBL" && kin != "Minkowski")
    GYOTO_ERROR
      ("DynamicalDisk3D::metric(): metric must be KerrBL");
  Disk3D::metric(gg);
}

void DynamicalDisk3D::file(std::string const &content) {
#ifdef GYOTO_USE_CFITSIO
    int withopacity=0;

    dirname_ = new char[strlen(content.c_str())+1];
    strcpy(dirname_,content.c_str());
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dirname_)) == NULL) {
      GYOTO_ERROR("In DynamicalDisk3D.C constructor : bad dirname_");
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
      GYOTO_ERROR("In DynamicalDisk3D.C: bad nb_times_ value");

    //check whether absorption is provided
    {
      // make declarations local (avoid warning)
      ostringstream stream_name ;
      stream_name << dirname_ << "data3D0001.fits.gz"; 
      string filename = stream_name.str();
      fitsRead(filename);
    }
    if (opacity()) {
      //cout << "WITH OPACITY" << endl;
      withopacity=1;
    }

    //initialize emission, absorption, velocity arrays
    emission_array_ = new double*[nb_times_] ;
    if (withopacity) absorption_array_ = new double*[nb_times_] ;
    else absorption_array_=NULL;
    velocity_array_ = new double*[nb_times_] ;

    double nu0b=0., zminb=0., zmaxb=0., rinb=0., routb=0.;
    size_t nnub=0, nphib=0, nzb=0, nrb=0;
    
    //fill in the arrays
    for (int i=1; i<=nb_times_; i++) {
      ostringstream stream_name ;
      stream_name << dirname_ << "data3D" 
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

      //save emission
      if (getEmissquant()){
	double * emtemp = const_cast<double*>(getEmissquant());
	emission_array_[i-1] = new double[nel1];
	for (size_t j=0;j<nel1;j++)
	  emission_array_[i-1][j]=emtemp[j];
      }else {
	GYOTO_ERROR("In DynamicalDisk3D::file(fname): "
		   "Emission must be supplied");
      }

      //save absorption (if any)
      if (withopacity){
	if (opacity()){
	  double * abstemp = const_cast<double*>(opacity());
	  absorption_array_[i-1] = new double[nel1];
	  for (size_t j=0;j<nel1;j++)
	    absorption_array_[i-1][j]=abstemp[j];
	  //cout << "SAVING ABS ARRAY" << endl;
	}else{
	  GYOTO_ERROR("In DynamicalDisk3D::file(fname): "
		     "Absorption should be supplied here");
	}
      }

      //save velocity
      if (getVelocity()){
	double * veltemp = const_cast<double*>(getVelocity());
	velocity_array_[i-1] = new double[nel2];
	for (size_t j=0;j<nel2;j++)
	  velocity_array_[i-1][j]=veltemp[j];
      }else{
	GYOTO_ERROR("In DynmicalDisk::file(fname): "
		   "Velocity must be supplied");
      }
      
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
	  ) GYOTO_ERROR("DynamicalDisk3D::file(fname) Grid is not constant!");
    }
#else
    GYOTO_ERROR("This Gyoto has no FITS i/o"); 
#endif     
}
std::string DynamicalDisk3D::file() const {return dirname_;}

void DynamicalDisk3D::tinit(double t) {tinit_=t;}
double DynamicalDisk3D::tinit()const{return tinit_;}

void DynamicalDisk3D::dt(double t) {dt_=t;}
double DynamicalDisk3D::dt()const{return dt_;}

void DynamicalDisk3D::PLindex(double t) {PLindex_=t;}
double DynamicalDisk3D::PLindex()const{return PLindex_;}

void DynamicalDisk3D::floorTemperature(double t) {floortemperature_=t;}
double DynamicalDisk3D::floorTemperature()const{return floortemperature_;}

void DynamicalDisk3D::temperature(bool t) {temperature_=t;}
bool DynamicalDisk3D::temperature() const {return temperature_;}

void DynamicalDisk3D::withVelocity(bool t) {novel_=!t;}
bool DynamicalDisk3D::withVelocity() const {return !novel_;}
