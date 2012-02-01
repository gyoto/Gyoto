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
#include "GyotoDisk3D.h"
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

Disk3D::Disk3D() :
  Generic("Disk3D"), filename_(""),
  emissquant_(NULL), velocity_(NULL),
  dnu_(1.), nu0_(0), nnu_(0),
  dphi_(0.), nphi_(0), repeat_phi_(1),
  dz_(0.), zmin_(-DBL_MAX), nz_(0), zmax_(DBL_MAX),
  dr_(0.), rin_(-DBL_MAX), nr_(0), rout_(DBL_MAX)
{
  GYOTO_DEBUG << "Disk3D Construction" << endl;
}

Disk3D::Disk3D(const Disk3D& o) :
  Generic(o), filename_(o.filename_),
  emissquant_(NULL), velocity_(NULL),
  dnu_(o.dnu_), nu0_(o.nu0_), nnu_(o.nnu_),
  dphi_(o.dphi_), nphi_(o.nphi_), repeat_phi_(o.repeat_phi_),
  dz_(o.dz_), zmin_(o.zmin_), nz_(o.nz_), zmax_(o.zmax_),
  dr_(o.dr_), rin_(o.rin_), nr_(o.nr_), rout_(o.rout_)
{
  GYOTO_DEBUG << "Disk3D Copy" << endl;
  size_t ncells = 0;
  if (o.emissquant_) {
    emissquant_ = new double[ncells = nnu_ * nphi_ * nz_ * nr_];
    memcpy(emissquant_, o.emissquant_, ncells * sizeof(double));
  }
  if (o.velocity_) {
    velocity_ = new double[ncells = 3 * nphi_ * nz_ * nr_];
    memcpy(velocity_, o.velocity_, ncells * sizeof(double));
  }
}
Disk3D* Disk3D::clone() const
{ return new Disk3D(*this); }

Disk3D::~Disk3D() {
  GYOTO_DEBUG << "Disk3D Destruction" << endl;
  if (emissquant_) delete [] emissquant_;
  if (velocity_) delete [] velocity_;
}

void Disk3D::setEmissquant(double * pattern) {
  emissquant_ = pattern;
}

void Disk3D::setVelocity(double * pattern) {
  velocity_ = pattern;
}

void Disk3D::copyEmissquant(double const *const pattern, size_t const naxes[4]) {
  GYOTO_DEBUG << endl;
  if (emissquant_) {
    GYOTO_DEBUG << "delete [] emissquant_;" << endl;
    delete [] emissquant_; emissquant_ = NULL;
  }
  if (pattern) {
    size_t nel;
    if (nphi_ != naxes[1]) {
      GYOTO_DEBUG <<"nphi_ changed, freeing velocity_" << endl;
      if (velocity_) { delete [] velocity_; velocity_= NULL; }
    }
    if (nz_ != naxes[2]) {
      GYOTO_DEBUG <<"nz_ changed, freeing velocity_" << endl;
      if (velocity_) { delete [] velocity_; velocity_= NULL; }
    }
    if (nr_ != naxes[3]) {
      GYOTO_DEBUG <<"nr_ changed, freeing velocity_" << endl;
      if (velocity_) { delete [] velocity_; velocity_= NULL; }
    }
    if (!(nel=(nnu_ = naxes[0]) * (nphi_=naxes[1]) * (nz_=naxes[2]) * (nr_=naxes[3])))
      throwError( "dimensions can't be null");
    dr_ = (rout_ - rin_) / nr_;
    dz_ = (zmax_ - zmin_) / nz_;
    dphi_ = 2.*M_PI/double(nphi_*repeat_phi_);
    GYOTO_DEBUG << "allocate emissquant_;" << endl;
    emissquant_ = new double[nel];
    GYOTO_DEBUG << "pattern >> emissquant_" << endl;
    memcpy(emissquant_, pattern, nel*sizeof(double));
  }
}

double const * const Disk3D::getEmissquant() const { return emissquant_; }
void Disk3D::getEmissquantNaxes( size_t naxes[3] ) const
{ 
  naxes[0] = nnu_; naxes[1] = nphi_; naxes[2] = nz_; 
  naxes[3] = nr_;
}

void Disk3D::copyVelocity(double const *const velocity, size_t const naxes[3]) {
  GYOTO_DEBUG << endl;
  if (velocity_) {
    GYOTO_DEBUG << "delete [] velocity_;\n";
    delete [] velocity_; velocity_ = NULL;
  }
  if (velocity) {
    if (!emissquant_) throwError("Please use copyEmissquant() before copyVelocity()");
    if (nphi_ != naxes[0] || nz_ != naxes[1] || nr_ != naxes[2])
      throwError("emissquant_ and velocity_ have inconsistent dimensions");
    GYOTO_DEBUG << "allocate velocity_;" << endl;
    velocity_ = new double[3*nphi_*nz_*nr_];
    GYOTO_DEBUG << "velocity >> velocity_" << endl;
    memcpy(velocity_, velocity, 3*nphi_*nz_*nr_*sizeof(double));
  }
}
double const * const Disk3D::getVelocity() const { return velocity_; }

void Disk3D::repeatPhi(size_t n) {
  repeat_phi_ = n;
  dphi_=2.*M_PI/double(nphi_*repeat_phi_);
}
size_t Disk3D::repeatPhi() const { return repeat_phi_; }

void Disk3D::nu0(double freq) { nu0_ = freq; }
double Disk3D::nu0() const { return nu0_; }

void Disk3D::dnu(double dfreq) { dnu_ = dfreq; }
double Disk3D::dnu() const { return dnu_; }

void Disk3D::rin(double rrin) {
  rin_ = rrin;
  if (nr_) dr_ = (rout_-rin_) / nr_;
}
double Disk3D::rin() {return rin_;}

void Disk3D::dr(double ddr) {dr_ = ddr;}
double Disk3D::dr() {return dr_;}

void Disk3D::rout(double rrout) {
  rout_ = rrout;
  if (nr_) dr_ = (rout_-rin_) / nr_;
}
double Disk3D::rout() {return rout_;}

void Disk3D::zmin(double zzmin) {
  zmin_ = zzmin;
  if (nz_) dz_ = (zmax_-zmin_) / nz_;
}
double Disk3D::zmin() {return zmin_;}

void Disk3D::dz(double ddz) {dz_ = ddz;}
double Disk3D::dz() {return dz_;}

void Disk3D::zmax(double zzmax) {
  zmax_ = zzmax;
  if (nz_) dz_ = (zmax_-zmin_) / nz_;
}
double Disk3D::zmax() {return zmax_;}


void Disk3D::fitsRead(string filename) {
  GYOTO_MSG << "Disk3D reading FITS file: " << filename << endl;

  filename_ = filename;
  char*     pixfile   = const_cast<char*>(filename_.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  int       anynul    = 0;
  long      tmpl;
  double    tmpd;
  long      naxes []  = {1, 1, 1, 1};
  long      fpixel[]  = {1,1,1,1};
  long      inc   []  = {1,1,1,1};
  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  GYOTO_DEBUG << "Disk3D::readFile(): opening file" << endl;
  if (fits_open_file(&fptr, pixfile, 0, &status)) throwCfitsioError(status) ;

  ////// READ FITS KEYWORDS COMMON TO ALL TABLES ///////

  GYOTO_DEBUG << "Disk3D::readFile(): read RepeatPhi_" << endl;
  fits_read_key(fptr, TLONG, "GYOTO Disk3D RepeatPhi", &tmpl,
		NULL, &status);
  if (status) {
    if (status == KEY_NO_EXIST) status = 0; // not fatal
    else throwCfitsioError(status) ;
  } else repeat_phi_ = size_t(tmpl); // RepeatPhi found

  GYOTO_DEBUG << "Disk3D::readFile(): read Rin_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO Disk3D Rin", &tmpd,
		NULL, &status);
  if (status) {
    //if (status == KEY_NO_EXIST) status = 0; // not fatal -> now fatal in 3D
    throwCfitsioError(status) ;
  } else {
    rin_ = tmpd; // InnerRadius found
  }
  GYOTO_DEBUG << "Disk3D::readFile(): read Rout_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO Disk3D Rout", &tmpd,
		NULL, &status);
  if (status) {
    //if (status == KEY_NO_EXIST) status = 0; // not fatal -> now fatal in 3D
    throwCfitsioError(status) ;
  } else {
    rout_ = tmpd; // OuterRadius found
  }

  GYOTO_DEBUG << "Disk3D::readFile(): read Zmin_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO Disk3D Zmin", &tmpd,
		NULL, &status);
  if (status) {
    throwCfitsioError(status) ;
  } else {
    zmin_ = tmpd; // Zmin found
  }
  GYOTO_DEBUG << "Disk3D::readFile(): read Zmax_" << endl;
  fits_read_key(fptr, TDOUBLE, "GYOTO Disk3D Zmax", &tmpd,
		NULL, &status);
  if (status) {
    throwCfitsioError(status) ;
  } else {
    zmax_ = tmpd; // Zmax found
  }

  ////// FIND MANDATORY EMISSION HDU, READ KWDS & DATA ///////
  GYOTO_DEBUG << "Disk3D::readFile(): search emissquant HDU" << endl;
  if (fits_movnam_hdu(fptr, ANY_HDU,
		      const_cast<char*>("GYOTO Disk3D emissquant"),
		      0, &status))
    throwCfitsioError(status) ;
  GYOTO_DEBUG << "Disk3D::readFile(): get image size" << endl;
  if (fits_get_img_size(fptr, 4, naxes, &status)) throwCfitsioError(status) ;

  //update nu0_, nnu_, dnu_;
  nnu_ = naxes[0]; 
  double CRPIX1;
  GYOTO_DEBUG << "Disk3D::readFile(): read CRPIX1, CRVAL1, CDELT1"
		    << endl;
  fits_read_key(fptr, TDOUBLE, "CRVAL1", &nu0_, NULL, &status);
  fits_read_key(fptr, TDOUBLE, "CDELT1", &dnu_, NULL, &status);
  fits_read_key(fptr, TDOUBLE, "CRPIX1", &CRPIX1, NULL, &status);
  if (status) throwCfitsioError(status) ;
  if (CRPIX1 != 1) nu0_ -= dnu_*(CRPIX1 - 1.);

  // update nphi_, dphi_
  nphi_ = naxes[1];
  dphi_ = 2.*M_PI/double(nphi_*repeat_phi_);

  // update nz_, nr_, dz_, dr_
  nz_ = naxes[2];
  nr_ = naxes[3];
  dr_ = (rout_ - rin_) / nr_;
  dz_ = (zmax_ - zmin_) / nz_;

  if (emissquant_) { delete [] emissquant_; emissquant_ = NULL; }
  emissquant_ = new double[nnu_ * nphi_ * nz_ * nr_];
  if (debug())
    cerr << "Disk3D::readFile(): read emission: "
	 << "nnu_=" << nnu_ << ", nphi_="<<nphi_ << ", nz_="<<nz_ << ", nr_="<<nr_ << "...";
  if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc,
		       0, emissquant_,&anynul,&status)) {
    GYOTO_DEBUG << " error, trying to free pointer" << endl;
    delete [] emissquant_; emissquant_=NULL;
    throwCfitsioError(status) ;
  }
  GYOTO_DEBUG << " done." << endl;

  ////// FIND MANDATORY VELOCITY HDU ///////

  fits_movnam_hdu(fptr, ANY_HDU,
		  const_cast<char*>("GYOTO Disk3D velocity"),
		  0, &status);
  if (status) {
    if (status == BAD_HDU_NUM)
      throwCfitsioError(status) ;
  } else {
    if (fits_get_img_size(fptr, 4, naxes, &status)) throwCfitsioError(status) ;
    if (   size_t(naxes[0]) != size_t(3)
	   || size_t(naxes[1]) != nphi_
	   || size_t(naxes[2]) != nz_
	   || size_t(naxes[3]) != nr_)
      throwError("Disk3D::readFile(): velocity array not conformable");
    if (velocity_) { delete [] velocity_; velocity_ = NULL; }
    velocity_ = new double[3 * nphi_ * nz_ * nr_];
    if (fits_read_subset(fptr, TDOUBLE, fpixel, naxes, inc, 
			 0, velocity_,&anynul,&status)) {
      delete [] velocity_; velocity_=NULL;
      throwCfitsioError(status) ;
    }
  }

  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}

void Disk3D::fitsWrite(string filename) {
  if (!emissquant_) throwError("Disk3D::fitsWrite(filename): nothing to save!");
  filename_ = filename;
  char*     pixfile   = const_cast<char*>(filename_.c_str());
  fitsfile* fptr      = NULL;
  int       status    = 0;
  long      naxes []  = {nnu_, nphi_, nz_, nr_};
  long      fpixel[]  = {1,1,1,1};
  char * CNULL=NULL;

  char      ermsg[31] = ""; // ermsg is used in throwCfitsioError()

  ////// CREATE FILE
  GYOTO_DEBUG << "creating file" << endl;
  fits_create_file(&fptr, pixfile, &status);
  fits_create_img(fptr, DOUBLE_IMG, 4, naxes, &status);
  if (status) throwCfitsioError(status) ;

  ////// WRITE FITS KEYWORDS COMMON TO ALL TABLES ///////
  if (repeat_phi_!=1)
    fits_write_key(fptr, TLONG,
		   const_cast<char*>("GYOTO Disk3D RepeatPhi"),
		   &repeat_phi_, CNULL, &status);

  if (rin_ > -DBL_MAX){
    fits_write_key(fptr, TDOUBLE,
		   const_cast<char*>("GYOTO Disk3D Rin"),
		   &rin_, CNULL, &status);
  }else{
    cout << "Disk3D::fitsRead Error rin_ not set!" << endl;
    status=1;
    throwCfitsioError(status)
  }

  if (rout_ < DBL_MAX){
    fits_write_key(fptr, TDOUBLE,
		   const_cast<char*>("GYOTO Disk3D Rout"),
		   &rout_, CNULL, &status);
  }else{
    cout << "Disk3D::fitsRead Error rout_ not set!" << endl;
    status=1;
    throwCfitsioError(status)
  }

  if (zmin_ > -DBL_MAX){
    fits_write_key(fptr, TDOUBLE,
		   const_cast<char*>("GYOTO Disk3D Zmin"),
		   &zmin_, CNULL, &status);
  }else{
    cout << "Disk3D::fitsRead Error zmin_ not set!" << endl;
    status=1;
    throwCfitsioError(status)
  }

  if (zmax_ < DBL_MAX){
    fits_write_key(fptr, TDOUBLE,
		   const_cast<char*>("GYOTO Disk3D Zmax"),
		   &zmax_, CNULL, &status);
  }else{
    cout << "Disk3D::fitsRead Error zmax_ not set!" << endl;
    status=1;
    throwCfitsioError(status)    
  }

  ////// SAVE EMISSION IN PRIMARY HDU ///////
  GYOTO_DEBUG << "saving emissquant_\n";
  fits_write_key(fptr, TSTRING,
		 const_cast<char*>("EXTNAME"),
		 const_cast<char*>("GYOTO Disk3D emissquant"),
		 CNULL, &status);
  fits_write_key(fptr, TDOUBLE,
		 const_cast<char*>("CRVAL1"),
		 &nu0_, CNULL, &status);
  fits_write_key(fptr, TDOUBLE,
		 const_cast<char*>("CDELT1"),
		 &dnu_, CNULL, &status);
  double CRPIX1 = 1.;
  fits_write_key(fptr, TDOUBLE,
		 const_cast<char*>("CRPIX1"),
		 &CRPIX1, CNULL, &status);
  fits_write_pix(fptr, TDOUBLE, fpixel, nnu_*nphi_*nz_*nr_, emissquant_, &status);
  if (status) throwCfitsioError(status) ;

  ////// SAVE MANDATORY VELOCITY HDU ///////
  if (velocity_) {
    GYOTO_DEBUG << "saving velocity_\n";
    naxes[0]=3;
    fits_create_img(fptr, DOUBLE_IMG, 4, naxes, &status);
    fits_write_key(fptr, TSTRING, const_cast<char*>("EXTNAME"),
		   const_cast<char*>("GYOTO Disk3D velocity"),
		   CNULL, &status);
    fits_write_pix(fptr, TDOUBLE, fpixel, 3*nphi_*nz_*nr_, velocity_, &status);
    if (status) throwCfitsioError(status) ;
  }

  ////// CLOSING FILE ///////
  GYOTO_DEBUG << "close FITS file\n";
  if (fits_close_file(fptr, &status)) throwCfitsioError(status) ;
  fptr = NULL;
}

void Disk3D::getIndices(size_t i[4], double const co[4], double nu) const {
  GYOTO_DEBUG << "dnu_="<<dnu_<<", dphi_="<<dphi_
	      <<", dz_="<<dz_<<", dr_="<<dr_<<endl;
  if (nu <= nu0_) i[0] = 0;
  else {
    i[0] = size_t((nu-nu0_)/dnu_);
    if (i[0] >= nnu_) i[0] = nnu_-1;
  }
  
  double rr,zz,phi, tt=co[0]; //cylindrical coord
  switch (gg_ -> getCoordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    double rs=co[1];
    zz= rs*cos(co[2]);
    rr= sqrt(rs*rs-zz*zz);
    phi= co[3];
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    zz= co[3];
    double xx=co[1], yy=co[2];
    rr= sqrt(xx*xx+yy*yy);
    phi= atan2(yy, xx);
    break;
  default:
    throwError("Disk3D::getIndices(): unknown COORDKIND");
  }

  //Phi indice
  while (phi<0) phi += 2.*M_PI;
  i[1] = size_t(phi/dphi_) % nphi_;

  //z indice
  if (zz<0. && zmin_>=0.) zz*=-1.; //if zmin>=0, assume disk is symmetric
  i[2] = size_t((zz-zmin_)/dz_);
  if (i[2] == nz_) i[2] = nz_ - 1;
  else if (i[2] > nz_) throwError("In Disk3D::getIndices() impossible indice value for z");

  //r indice
  i[3] = size_t((rr-rin_)/dr_);
  if (i[3] == nr_) i[3] = nr_ - 1;
  else if (i[3] > nr_) throwError("In Disk3D::getIndices() impossible indice value for r");

}

void Disk3D::getVelocity(double const pos[4], double vel[4]) {
    if (velocity_) {
      /*if (dir_ != 1)
      throwError("Disk3D::getVelocity(): "
      "dir_ should be 1 if velocity_ is provided");*/
    size_t i[4]; // {i_nu, i_phi, i_z, i_r}
    getIndices(i, pos);
    double phiprime=velocity_[i[3]*3*nphi_*nz_+i[2]*3*nphi_+i[1]*3+0];
    double zprime=velocity_[i[3]*3*nphi_*nz_+i[2]*3*nphi_+i[1]*3+1];
    double rprime=velocity_[i[3]*3*nphi_*nz_+i[2]*3*nphi_+i[1]*3+2];
    switch (gg_->getCoordKind()) {
    case GYOTO_COORDKIND_SPHERICAL:
      {
	//Formula from derivation of rsph^2=rcyl^2+zz^2
	//and rsph*cos(th)=zz where rsph and rcyl
	//are spherical and cylindrical radii
	double rsph=pos[1], th=pos[2], ph=pos[3],
	  zz= rsph*cos(th), rcyl=sqrt(rsph*rsph-zz*zz);
	vel[1] = (rcyl*rprime+zz*zprime)/rsph;
	vel[2] = (vel[1]*cos(th)-zprime)/(rsph*sin(th));
	vel[3] = phiprime;
	vel[0] = gg_->SysPrimeToTdot(pos, vel+1);
	vel[1] *= vel[0];
	vel[2] *= vel[0];
	vel[3] *= vel[0];
      }
      break;
    case GYOTO_COORDKIND_CARTESIAN:
      throwError("Disk3D::getVelocity(): metric must be in "
		 "spherical coordinates");
      break;
    default:
      throwError("Disk3D::getVelocity(): unknown COORDKIND");
    } 
    }else throwError("In Disk3D::getVelocity(): velocity_==NULL!");
}

int Disk3D::Impact(Photon *ph, size_t index,
			       Astrobj::Properties *data) {
  /*  size_t indices[4];
  getIndices()
  size_t iphi=71, iz=71, ir=191;
  double phiprime=velocity_[ir*3*nphi_*nz_+iz*3*nphi_+iphi*3+0];
  double zprime=velocity_[ir*3*nphi_*nz_+iz*3*nphi_+iphi*3+1];
  double rprime=velocity_[ir*3*nphi_*nz_+iz*3*nphi_+iphi*3+2];
  cout << "v term= " << phiprime << " " << zprime << " " << rprime << endl;
  double emiss=emissquant_[ir*nphi_*nz_+iz*nphi_+iphi];
  cout << "em term= " << emiss << endl;
  double rr=rin_+ir*dr_;
  double zz=zmin_+iz*dz_;
  double phi=iphi*dphi_;
  cout << "pos= " << phi << " " << zz << " " << rr << endl;
  return 0;*/

  if (debug()) 
    cerr << "DEBUG: in DiskFromFile3D::Impact()" << endl;
  double coord_ph_hit[8], coord_obj_hit[8];
  double rcross, phicross;
  double tmp, frac;
  double coord1[8], coord2[8];
 
  ph->getCoord(index, coord1);
  ph->getCoord(index+1, coord2);

  /* Transforms theta and phi in coord1,2 so that 
     theta is in [0,pi] and phi in [0,2pi] */
  double thetatmp=coord1[2], phitmp=coord1[3];
  while (thetatmp>M_PI) thetatmp-=2.*M_PI;
  while (thetatmp<-M_PI) thetatmp+=2.*M_PI;//then theta in [-pi,pi]
  if (thetatmp<0.) {
    thetatmp=-thetatmp;//then theta in [0,pi]
    phitmp+=M_PI;//thus, same point x,y,z
  }
  while (phitmp>2.*M_PI) phitmp-=2.*M_PI;
  while (phitmp<0.) phitmp+=2.*M_PI;//then phi in [0,2pi]
  coord1[2]=thetatmp;
  coord1[3]=phitmp;
  thetatmp=coord2[2], phitmp=coord2[3];
  while (thetatmp>M_PI) thetatmp-=2.*M_PI;
  while (thetatmp<-M_PI) thetatmp+=2.*M_PI;//then theta in [-pi,pi]
  if (thetatmp<0.) {
    thetatmp=-thetatmp;//then theta in [0,pi]
    phitmp+=M_PI;//thus, same point x,y,z
  }
  while (phitmp>2.*M_PI) phitmp-=2.*M_PI;
  while (phitmp<0.) phitmp+=2.*M_PI;//then phi in [0,2pi]
  coord2[2]=thetatmp;
  coord2[3]=phitmp;

  double theta1=coord1[2], theta2=coord2[2];

  // HEURISTIC TESTS TO PREVENT TOO MANY INTEGRATION STEPS
  // Speeds up a lot!
  // Idea: no test if r1,r2 > factr*rdiskmax_ AND z1,z2 have same sign
  double r1=coord1[1], r2=coord2[1],
    z1=r1*cos(theta1), z2=r2*cos(theta2);
  double factr=2.;
  double rtol=factr*rout_;
  if (coord1[1]>rtol && coord2[1]>rtol && z1*z2>0.)
    return 0;

  int iicross=-1;

  double t1=coord1[0], t2=coord2[0];

  double deltatmin=0.1, deltat12=fabs(t2-t1)*0.1;
  //Break the worldline in pieces of "size" deltat:
  double deltat= deltat12 < deltatmin ? deltat12 : deltatmin;
  double tcur=t2, myrcur=coord2[1], thetacur=coord2[2], phicur=coord2[3], 
    zcur=myrcur*cos(thetacur),rcur=sqrt(myrcur*myrcur-zcur*zcur),
    tdotcur, rdotcur, thetadotcur, phidotcur;
  //NB: myrcur is r in spherical coord, rcur is r in cylindrical coord
  int hit2=0, hit1=0;
  int i1=0,i1bef,i2,i2bef,i3,i3bef,i4,i5,i6,i7,i8;
  double phitab, ztab, rtab;

  //Flags that will be 1 if phicur<phimin, phicur>phimax, zcur<zmin:
  int philittle=0,phibig=0,zlittle=0;
  
  double tin, tout;//time of entry, exit from the disk
  
  //cout << "t1,t2, r1,r2= " << t1 << " " << t2 << " " << coord1[1] << " " << coord2[1] << endl;
  
  //Following while loop determines (if any) the entry and out points
  //along the geodesic that goes throw the disk
  
  //cout << "-----> Find entry: " << endl;
  // FIND ENTRY POINT (if any)
  while (tcur>t1+deltat 
	 && 
	 (
	  ((zmin_<0. && zcur<zmin_) || (zmin_>=0. && zcur<-zmax_)) 
	  || zcur>zmax_ || rcur>rout_ || rcur<rin_)
	 ){
    //NB: condition on zmin assumes disk is symmetric in z if zmin>=0
    //NBB: tcur must stay >t1+deltat coz getCoord can be called for:
    //t1<=tcur<=t2
    
    tcur-=deltat;
    ph -> getCoord( &tcur, 1, &myrcur, &thetacur, &phicur,
		    &tdotcur, &rdotcur, &thetadotcur, &phidotcur);
    //Cylindrical z and r coordinates
    zcur=myrcur*cos(thetacur);
    rcur=sqrt(myrcur*myrcur-zcur*zcur);
    //cout << "t1,tcur,t2=" << t1 << " " << tcur << " " << t2 << endl;
    
  }
  
  if (tcur>t1+deltat) {
    //HIT
    tin=tcur;
  }else{
    //NO HIT
    return 0;
  }
  
  // INTEGRATION BETWEEN ENTRY AND EXIT
  int indisk=1;
  //cout << "-----> Integrate!" << endl;
  while (indisk && tcur>t1+deltat){//Continue integration until either 
                                   //we get out of the disk, or time
                                   //becomes smaller than t1 
    tcur-=deltat;
    //cout << "New point at t1,tcur,t2 -- rcur= " << t1 << " " << tcur << " " << t2 << " -- " << rcur << endl;
    ph -> getCoord( &tcur, 1, &myrcur, &thetacur, &phicur,
		    &tdotcur, &rdotcur, &thetadotcur, &phidotcur);
    //Cylindrical z and r coordinates
    zcur=myrcur*cos(thetacur);
    rcur=sqrt(myrcur*myrcur-zcur*zcur);
    if (
	((zmin_<0. && zcur<zmin_) || (zmin_>=0. && zcur<-zmax_))
	|| zcur>zmax_ || rcur>rout_ || rcur<rin_
	){
      indisk=0;
    }else{ //Compute emission loop
      //cout << "Still in disk" << endl;
      //Insure phi is in 0-2pi
      thetatmp=thetacur, phitmp=phicur;
      while (thetatmp>M_PI) thetatmp-=2.*M_PI;
      while (thetatmp<-M_PI) thetatmp+=2.*M_PI;//then theta in [-pi,pi]
      if (thetatmp<0.) {
      thetatmp=-thetatmp;//then theta in [0,pi]
      phitmp+=M_PI;//thus, same point x,y,z
      }
      while (phitmp>2.*M_PI) phitmp-=2.*M_PI;
      while (phitmp<0.) phitmp+=2.*M_PI;//then phi in [0,2pi]
      thetacur=thetatmp;
      phicur=phitmp;

      coord_ph_hit[0]=tcur;
      coord_ph_hit[1]=myrcur;//spherical radius
      coord_ph_hit[2]=thetacur;
      coord_ph_hit[3]=phicur;
      coord_ph_hit[4]=tdotcur;
      coord_ph_hit[5]=rdotcur;
      coord_ph_hit[6]=thetadotcur;
      coord_ph_hit[7]=phidotcur;

      for (int ii=0;ii<4;ii++) coord_obj_hit[ii]=coord_ph_hit[ii];
      getVelocity(coord_obj_hit, coord_obj_hit+4);
      
      processHitQuantities(ph, coord_ph_hit, coord_obj_hit, deltat, data);
      if (!flag_radtransf_) indisk=0;//not to go on integrating
      
    } //end compute emission else loop
  }

  return 1;

}

double Disk3D::emission(double nu, double dsem,
				    double *,
				    double co[8]) const{
  return 1.;

}

double Disk3D::transmission(double nu, double dsem, double*co) const {
  return 0.;
}

int Disk3D::setParameter(std::string name, std::string content) {
  if      (name == "File")          fitsRead( content );
  else return Generic::setParameter(name, content);
  return 0;
}

#ifdef GYOTO_USE_XERCES
void Disk3D::fillElement(FactoryMessenger *fmp) const {
  fmp->setParameter("File", (filename_.compare(0,1,"!") ?
			     filename_ :
			     filename_.substr(1)));
  Generic::fillElement(fmp);
}

void Disk3D::setParameters(FactoryMessenger* fmp) {
  string name, content;
  setMetric(fmp->getMetric());
  while (fmp->getNextParameter(&name, &content)) {
    if  (name == "File") setParameter(name, fmp -> fullPath(content));
    else setParameter(name, content);
  }
}
#endif
