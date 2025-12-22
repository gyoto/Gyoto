/*
    Copyright 2024 Aimar Nicolas, Irene Urso

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

#include "GyotoSimBridge.h"
#include "GyotoFactoryMessenger.h"
#include "GyotoKerrBL.h"

#ifdef GYOTO_USE_CFITSIO
#include <fitsio.h>
#define throwCfitsioError(status) \
    { fits_get_errstatus(status, ermsg); GYOTO_ERROR(ermsg); }
#endif

using namespace std;
using namespace Gyoto;
using namespace Gyoto::Astrobj;

#include "GyotoProperty.h"
GYOTO_PROPERTY_START(SimBridge)
GYOTO_PROPERTY_STRING(SimBridge, Directory, directory)
GYOTO_PROPERTY_STRING(SimBridge, Filename, filename)
GYOTO_PROPERTY_DOUBLE(SimBridge, GammaMin, gammaMin)
GYOTO_PROPERTY_DOUBLE(SimBridge, GammaMax, gammaMax)
GYOTO_PROPERTY_BOOL(SimBridge, CircularMotion, NoCircularMotion, circularMotion)
GYOTO_PROPERTY_BOOL(SimBridge, CunninghamVel, NoCunninghamVel, cunninghamVel)
GYOTO_PROPERTY_DOUBLE(SimBridge, PLindex, PLindex)
GYOTO_PROPERTY_DOUBLE(SimBridge, FloorTemperature, floorTemperature)
GYOTO_PROPERTY_STRING(SimBridge, EmissionType, emissionType)
GYOTO_PROPERTY_STRING(SimBridge, BoundaryConditions, boundaryConditions)
GYOTO_PROPERTY_DOUBLE(SimBridge, Magnetization, magnetization)
GYOTO_PROPERTY_END(SimBridge, Standard::properties)

SimBridge::SimBridge() :
  Standard("SimBridge"),
  FitsRW(),
  spectrumBB_(NULL),
  spectrumKappaSynch_(NULL),
  spectrumPLSynch_(NULL),
  spectrumThermalSynch_(NULL),
  dirname_("None"),
  fname_("data"),
  temperature_(true),
  circularmotion_(false),
  cunninghamvel_(false),
  BinFile_(true),
  emission_("None"),
  PLindex_(1.),
  gammaMin_(1.),
  gammaMax_(1.),
  magnetizationParameter_(1.),
  floortemperature_(0.),
  time_array_(NULL),
  x1_array_(NULL),
  x2_array_(NULL),
  x3_array_(NULL),
  nu_array_(NULL),
  ntime_(1),
  nx1_(1),
  nx2_(1),
  nx3_(1),
  nnu_(0),
  boundCond_(NULL)
{
  boundCond_ = new string[5];
  boundCond_[0]=boundCond_[1]=boundCond_[2]=boundCond_[3]=boundCond_[4]="None";

  spectrumKappaSynch_ = new Spectrum::KappaDistributionSynchrotron();
  spectrumPLSynch_ = new Spectrum::PowerLawSynchrotron();
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
  spectrumBB_ = new Spectrum::BlackBody();
}

SimBridge::SimBridge(const SimBridge& orig) :
  Standard(orig),
  FitsRW(orig),
  spectrumBB_(NULL),
  spectrumKappaSynch_(NULL),
  spectrumPLSynch_(NULL),
  spectrumThermalSynch_(NULL),
  dirname_(orig.dirname_),
  fname_(orig.fname_),
  temperature_(orig.temperature_),
  circularmotion_(orig.circularmotion_),
  cunninghamvel_(orig.cunninghamvel_),
  BinFile_(orig.BinFile_),
  emission_(orig.emission_),
  PLindex_(orig.PLindex_),
  gammaMin_(orig.gammaMin_),
  gammaMax_(orig.gammaMax_),
  magnetizationParameter_(orig.magnetizationParameter_),
  floortemperature_(orig.floortemperature_),
  time_array_(NULL),
  x1_array_(NULL),
  x2_array_(NULL),
  x3_array_(NULL),
  nu_array_(NULL),
  ntime_(orig.ntime_),
  nx1_(orig.nx1_),
  nx2_(orig.nx2_),
  nx3_(orig.nx3_),
  nnu_(orig.nnu_),
  boundCond_(NULL)
{
  if (orig.boundCond_){
    boundCond_  = new string[5];
    for (int ii=0; ii<5; ii++){
      boundCond_[ii] = orig.boundCond_[ii];
    }
  }
  if (orig.time_array_){
    time_array_ = new double[ntime_];
    memcpy(time_array_, orig.time_array_, ntime_*sizeof(double));
  }

  if (orig.x1_array_){
    x1_array_ = new double[nx1_];
    memcpy(x1_array_, orig.x1_array_, nx1_*sizeof(double));
  }

  if (orig.x2_array_){
    x2_array_ = new double[nx2_];
    memcpy(x2_array_, orig.x2_array_, nx2_*sizeof(double));
  }

  if (orig.x3_array_){
    x3_array_ = new double[nx3_];
    memcpy(x3_array_, orig.x3_array_, nx3_*sizeof(double));
  }
  
  if (orig.nu_array_){
    nu_array_ = new double[nnu_];
    memcpy(nu_array_, orig.nu_array_, nnu_*sizeof(double));
  }

  if (orig.spectrumKappaSynch_()) spectrumKappaSynch_=orig.spectrumKappaSynch_->clone();
  if (orig.spectrumPLSynch_()) spectrumPLSynch_=orig.spectrumPLSynch_->clone();
  if (orig.spectrumThermalSynch_()) spectrumThermalSynch_=orig.spectrumThermalSynch_->clone();
  if (orig.spectrumBB_()) spectrumBB_=orig.spectrumBB_->clone();

  for (int ii=0; ii<4; ii++){
    emisInFile_[ii]=orig.emisInFile_[ii];
    absInFile_[ii]=orig.absInFile_[ii];
    if (ii!=3)
      rotInFile_[ii]=orig.rotInFile_[ii];
  }
}

SimBridge* SimBridge::clone() const { return new SimBridge(*this); }

SimBridge::~SimBridge() {
  if (debug()) cout << "DEBUG: SimBridge::~SimBridge()\n";
  if (boundCond_) delete [] boundCond_;

  if (time_array_) delete[] time_array_;
  if (x1_array_) delete[] x1_array_;
  if (x2_array_) delete[] x2_array_;
  if (x3_array_) delete[] x3_array_;
  if (nu_array_) delete[] nu_array_;

} 

string SimBridge::className() const { return  string("SimBridge"); }
string SimBridge::className_l() const { return  string("simbridge"); }

void SimBridge::directory(std::string const &d){
  dirname_=d;
}
std::string SimBridge::directory() const{
  return dirname_;
}

void SimBridge::filename(std::string const &f){
  fname_=f;
  if (dirname_=="None")
    GYOTO_ERROR("Please set the directory before the filenames");
  
  // Reading and save dimensions of arrays in header
  ostringstream stream_name ;
  stream_name << dirname_ << fname_ << setw(4) << setfill('0') << 0 << ".fits" ;
      
  string filename = stream_name.str();
  GYOTO_DEBUG << "Reading FITS file: " << filename << endl ;
  
  fitsfile* fptr = NULL;
  
  fptr = FitsRW::fitsOpen(filename);

  //  double* tmp;

  ntime_        = FitsRW::fitsReadKey(fptr, "NB_X0");
  time_array_   = new double[ntime_];
  FitsRW::fitsReadHDUData(fptr, "X0", time_array_, ntime_);
  GYOTO_DEBUG << "ntime_" << ntime_ << endl;
  GYOTO_DEBUG << "time array: " << time_array_[0] << " - " << time_array_[ntime_-1] << endl;
  
  nx1_          = FitsRW::fitsReadKey(fptr, "NB_X1");
  x1_array_     = new double[nx1_];
  FitsRW::fitsReadHDUData(fptr, "X1", x1_array_, nx1_);
  GYOTO_DEBUG << "nx1_" << nx1_ << endl;
  GYOTO_DEBUG << "x1 array: " << x1_array_[0] << " - " << x1_array_[nx1_-1] << endl;

  nx2_          = FitsRW::fitsReadKey(fptr, "NB_X2");
  x2_array_     = new double[nx2_];
  FitsRW::fitsReadHDUData(fptr, "X2", x2_array_, nx2_);
  GYOTO_DEBUG << "nx2_" << nx2_ << endl;
  GYOTO_DEBUG << "x2 array: " << x2_array_[0] << " - " << x2_array_[nx2_-1] << endl;

  nx3_         = FitsRW::fitsReadKey(fptr, "NB_X3");
  x3_array_    = new double[nx3_];
  FitsRW::fitsReadHDUData(fptr, "X3", x3_array_, nx3_);
  GYOTO_DEBUG << "nx3_" << nx3_ << endl;
  GYOTO_DEBUG << "x3 array: " << x3_array_[0] << " - " << x3_array_[nx3_-1] << endl;


  int status = 0;
  string key = "NB_FREQ";
  double tmpd;
  int tmpi;
  fits_movabs_hdu(fptr, 1, &tmpi, &status);
  fits_read_key(fptr, TDOUBLE, const_cast<char*>(key.c_str()), &tmpd, NULL, &status);
  if(status==0){
    nnu_       = FitsRW::fitsReadKey(fptr, "NB_FREQ");
    nu_array_  = new double[nnu_];
    FitsRW::fitsReadHDUData(fptr, "FREQ", nu_array_, nnu_);
    GYOTO_DEBUG << "nnu_" << nnu_ << endl;
    GYOTO_DEBUG << "nu array: " << nu_array_[0] << " - " << nu_array_[nnu_-1] << endl;
  }
  

  BinFile_ = FitsRW::fitsReadKey(fptr, "BINFILE");

  // Filling the flags
  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>("J_I"), 0, &status);
  if (status==0){
    emisInFile_[0]=1;
    temperature_=false;
  } else {
    emisInFile_[0]=0;
  }

  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>("J_Q"), 0, &status);
  emisInFile_[1]=status?0:1;

  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>("J_U"), 0, &status);
  emisInFile_[2]=status?0:1;

  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>("J_V"), 0, &status);
  emisInFile_[3]=status?0:1;

  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>("ALPHA_I"), 0, &status);
  absInFile_[0]=status?0:1;
  
  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>("ALPHA_Q"), 0, &status);
  absInFile_[1]=status?0:1;

  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>("ALPHA_U"), 0, &status);
  absInFile_[2]=status?0:1;

  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>("ALPHA_V"), 0, &status);
  absInFile_[3]=status?0:1;

  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>("R_Q"), 0, &status);
  rotInFile_[0]=status?0:1;
  
  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>("R_U"), 0, &status);
  rotInFile_[1]=status?0:1;

  fits_movnam_hdu(fptr, ANY_HDU, const_cast<char*>("R_V"), 0, &status);
  rotInFile_[2]=status?0:1;

  FitsRW::fitsClose(fptr);
}

std::string SimBridge::filename() const{
  return fname_;
}

void SimBridge::PLindex(double pl){
  PLindex_=pl;
}
double SimBridge::PLindex()const{
  return PLindex_;
}

void SimBridge::gammaMin(double gmin){
  gammaMin_=gmin;
}
double SimBridge::gammaMin() const{
  return gammaMin_;
}

void SimBridge::gammaMax(double gmax){
  gammaMax_=gmax;
}
double SimBridge::gammaMax() const{
  return gammaMax_;
}

void SimBridge::circularMotion(bool t){
  circularmotion_=t;
}
bool SimBridge::circularMotion() const{
  return circularmotion_;
}

void SimBridge::cunninghamVel(bool t){
  cunninghamvel_=t;
}
bool SimBridge::cunninghamVel() const{
  return cunninghamvel_;
}

void SimBridge::floorTemperature(double t){
  floortemperature_=t;
}
double SimBridge::floorTemperature()const{
  return floortemperature_;
}

void SimBridge::magnetization(double ss){
  magnetizationParameter_=ss;
}

double SimBridge::magnetization() const{
  return magnetizationParameter_;
}

void SimBridge::boundaryConditions(string const &sbc){
  std::string token;
  std::istringstream tokenStream(sbc);
  int wordCount = 0;

  while (std::getline(tokenStream, token, ' ') && wordCount < 5) {
    std::istringstream commaStream(token);
    while (std::getline(commaStream, token, ',') && wordCount < 5) {
      if (!token.empty()) {
        switch (wordCount) {
          case 0: boundCond_[0] = token; break;
          case 1: boundCond_[1] = token; break;
          case 2: boundCond_[2] = token; break;
          case 3: boundCond_[3] = token; break;
          case 4: boundCond_[4] = token; break;
        }
        ++wordCount;
      }
    }
  }
}
std::string SimBridge::boundaryConditions() const{
  std::string list;
  for (int ii=0;ii<5;ii++)
    list += boundCond_[ii] + " ";
  return list;
}

void SimBridge::emissionType(std::string const &kind){
  if(kind == "Thermal")
    emission_ = "Thermal";
  else if(kind == "Kappa")
    emission_ = "Kappa";
  else if (kind == "PL")
    emission_ = "PL";
  else if (kind == "BB")
    emission_ = "BlackBody";
  else
    GYOTO_ERROR("unknown electron distribution!");
}

std::string SimBridge::emissionType() const{
  return emission_;
}

void SimBridge::radiativeQ(double *Inu, double *Qnu, double *Unu, double *Vnu,
       Eigen::Matrix4d *Onu, double const *nuem , size_t nbnu, double dsem,
       state_t const &coord_ph, double const *coord_obj) const {

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  if (dirname_=="None")
      GYOTO_ERROR("In SimBridge RadiativeQ : dirname_ not defined, please use directory(string)");
  
  double tcur=coord_ph[0]; // in M units # TBC
  int nfile=ntime_==1?1:2; // 1 if ntime=1, 2 otherwise
  //cout << "tcur: " << tcur << endl;
  int index=getIndex(tcur); //index of files to be loaded (index and index+1)
  //cout << "index: " << index << endl;
  
  // Creating arrays
  long ncells = nx1_*nx2_*nx3_; // Number of cells for each time
  double* density_array=NULL;
  double* temperature_array=NULL;
  double** magneticfield_array=NULL;
  double** emission_array=NULL;
  double** absorption_array=NULL;
  double** rotation_array=NULL;
  if (temperature_){
    density_array = new double[nfile*ncells];
    temperature_array = new double[nfile*ncells];
    
    if (BinFile_){
      magneticfield_array = new double*[4];
      for (int ii=0; ii<4; ii++){
        magneticfield_array[ii] = new double[nfile*ncells];
      }
    }
  } else {
    emission_array = new double*[4];
    absorption_array = new double*[4];
    rotation_array = new double*[3];
    for (int ii=0; ii<4; ii++){
      if (emisInFile_[ii]){
        emission_array[ii] = new double[nfile*ncells];
      }else{
        emission_array[ii] = nullptr;
      }

      if (absInFile_[ii]) {absorption_array[ii] = new double[nfile*ncells];}
      else {absorption_array[ii] = nullptr;}
    }
    for (int ii=0; ii<3; ii++){
      if (rotInFile_[ii]){
        rotation_array[ii] = new double[nfile*ncells];
      }else{
        rotation_array[ii] = nullptr;
      }
    }
  }

  // Opening and reading Files
  //  double* tmp;
  double time_interpo[nfile];
  for (int ii=0; ii<nfile; ii++){
    ostringstream stream_name ;
    stream_name << dirname_ << fname_ << setw(4) << setfill('0') << index+ii << ".fits" ;
        
    string filename = stream_name.str();
    GYOTO_DEBUG << "Reading FITS file: " << filename << endl ;
    //cout << "Reading FITS file: " << filename << endl ;
    fitsfile* fptr = NULL;
    
    fptr = FitsRW::fitsOpen(filename);

    int ntime = FitsRW::fitsReadKey(fptr, "NB_X0");
    int nx1   = FitsRW::fitsReadKey(fptr, "NB_X1");
    int nx2   = FitsRW::fitsReadKey(fptr, "NB_X2");
    int nx3   = FitsRW::fitsReadKey(fptr, "NB_X3");

    if (ntime!=ntime_ || nx1!=nx1_ || nx2!=nx2_ || nx3!=nx3_)
      GYOTO_ERROR("In SimBridge RadiativeQ : size of arrays in FITS file different from initial file");

    time_interpo[ii] = FitsRW::fitsReadKey(fptr, "TIME");

    double* emplacement;

    if (temperature_){
      // The arrays stored in the FITS files are density and temperature + eventually magnetic field
      // Read NUMBERDENSITY
      emplacement = density_array + ii * ncells;
      FitsRW::fitsReadHDUData(fptr, "NUMBERDENSITY", emplacement, ncells);

      // Read TEMPERATURE
      emplacement = temperature_array + ii * ncells;
      FitsRW::fitsReadHDUData(fptr, "TEMPERATURE", emplacement, ncells);
      
      if (BinFile_) {
        for (int jj = 0; jj < 4; jj++) {
          std::ostringstream hdu_name;
          hdu_name << "B" << jj;
          emplacement = magneticfield_array[jj] + ii * ncells;
          FitsRW::fitsReadHDUData(fptr, hdu_name.str().c_str(), emplacement, ncells);
        }
      }
    }else{
      // The arrays stored in the FITS files are (at least) emission(s).
      if (emission_array==NULL or absorption_array==NULL or rotation_array==NULL)
        GYOTO_ERROR("In SimBridge RadiativeQ : emission, absorption and/or rotation arrays are not properly set.");
      if (nu_array_==NULL)
        GYOTO_ERROR("In SimBridge RadiativeQ : frequency array is not set");
      
      if (emisInFile_[0]){
        emplacement = emission_array[0] + ii * ncells;
        FitsRW::fitsReadHDUData(fptr, "J_I", emplacement, ncells);
      }
      if (emisInFile_[1]){
        emplacement = emission_array[1] + ii * ncells;
        FitsRW::fitsReadHDUData(fptr, "J_Q", emplacement, ncells);
      }
      if (emisInFile_[2]){
        emplacement = emission_array[2] + ii * ncells;
        FitsRW::fitsReadHDUData(fptr, "J_U", emplacement, ncells);
      }
      if (emisInFile_[3]){
        emplacement = emission_array[3] + ii * ncells;
        FitsRW::fitsReadHDUData(fptr, "J_V", emplacement, ncells);
      }
      if (absInFile_[0]){
        emplacement = absorption_array[0] + ii * ncells;
        FitsRW::fitsReadHDUData(fptr, "ALPHA_I", emplacement, ncells);
      }
      if (absInFile_[1]){
        emplacement = absorption_array[1] + ii * ncells;
        FitsRW::fitsReadHDUData(fptr, "ALPHA_Q", emplacement, ncells);
      }
      if (absInFile_[2]){
        emplacement = absorption_array[2] + ii * ncells;
        FitsRW::fitsReadHDUData(fptr, "ALPHA_U", emplacement, ncells);
      }
      if (absInFile_[3]){
        emplacement = absorption_array[3] + ii * ncells;
        FitsRW::fitsReadHDUData(fptr, "ALPHA_V", emplacement, ncells);
      }
      if (rotInFile_[0]){
        emplacement = rotation_array[0] + ii * ncells;
        FitsRW::fitsReadHDUData(fptr, "R_Q", emplacement, ncells);
      }
      if (rotInFile_[1]){
        emplacement = rotation_array[1] + ii * ncells;
        FitsRW::fitsReadHDUData(fptr, "R_U", emplacement, ncells);
      }
      if (rotInFile_[2]){
        emplacement = rotation_array[2] + ii * ncells;
        FitsRW::fitsReadHDUData(fptr, "R_V", emplacement, ncells);
      }
    }
    FitsRW::fitsClose(fptr);
    
    //cout << "Emission array" << endl;
    //for (int ii=0; ii<ncells; ii++){
    //  cout << emission_array[0][ii] << " ";
    //}
    //cout << endl;
  } // End of reading FITS files. At this point, the relevent arrays should be filled.

  double x1 = coord_ph[1], x2 = coord_ph[2], x3 = coord_ph[3];
  double pos[4] = {tcur, x1, x2, x3};
  double vel[4];
  //cout << "pos: " << tcur << " " << x1 << " " << x2*180/3.141592 << " " << x3*180/3.141592 << endl;
  const_cast<SimBridge*>(this)->getVelocity(pos, vel);
  //cout << "vel: " << vel[0] << " " << vel[1] << " " << vel[2] << " " << vel[3] << endl;

  // Defining emission, absoprtion and rotation coefficients for the transmission matrix
  GYOTO_DEBUG << " Initializing radiative coefficients arrays." << endl;
  double jInu[nbnu], jQnu[nbnu], jUnu[nbnu], jVnu[nbnu];
  double aInu[nbnu], aQnu[nbnu], aUnu[nbnu], aVnu[nbnu];
  double rotQnu[nbnu], rotUnu[nbnu], rotVnu[nbnu];
  for (size_t ii=0; ii<nbnu; ++ii){
    // Initialze them to -1 to create error if not updated
    jInu[ii]=-1.;
    jQnu[ii]=-1.;
    jUnu[ii]=-1.;
    jVnu[ii]=-1.;
    aInu[ii]=-1.;
    aQnu[ii]=-1.;
    aUnu[ii]=-1.;
    aVnu[ii]=-1.;
    rotQnu[ii]=-1.;
    rotUnu[ii]=-1.;
    rotVnu[ii]=-1.;
  }

  double Chi=0.;
  if (temperature_){
    // Interpolate number density and temperature
    double Xq[4] = {tcur, x1, x2, x3};
    int X_params[4] = {nfile, nx1_, nx2_, nx3_};
    double** X;
    X = new double*[4];
    X[0] = time_interpo;
    X[1] = x1_array_;
    X[2] = x2_array_;
    X[3] = x3_array_;

    //cout << "Xq : " << Xq[0] << ", " << Xq[1] << ", " << Xq[2] << ", " << Xq[3] << endl;
    double number_density = max(interpolate(4, density_array, Xq, X, X_params, boundCond_), 0.); // prevent negative values
    double temperature    = max(interpolate(4, temperature_array, Xq, X, X_params, boundCond_),floortemperature_);
    //cout << "ne, Te at (t,r,theta, phi) : "  << number_density << ", " << temperature << ", (" << tcur << "," << x1 << "," << x2 << "," << x3  << ")" << endl;

    int avg=0; // flag for magnetic field average for synchrotron
    if (!BinFile_ && magneticConfig_=="None")
      avg = 1;
    
    double B4vect[4]={0.,0.,0.,0.};
    double theta_mag, BB;
    if (!avg){
      if (BinFile_){
        B4vect[0] = interpolate(4, magneticfield_array[0], Xq, X, X_params, boundCond_);
        B4vect[1] = interpolate(4, magneticfield_array[1], Xq, X, X_params, boundCond_);
        B4vect[2] = interpolate(4, magneticfield_array[2], Xq, X, X_params, boundCond_);
        B4vect[3] = interpolate(4, magneticfield_array[3], Xq, X, X_params, boundCond_);
        
        // compute norm of 3D magnetic field
        double B4vectproj[4];
        for (int ii=0;ii<4;ii++)
          B4vectproj[ii] = B4vect[ii];
        gg_->projectFourVect(&coord_ph[0],B4vectproj,vel);
        BB = gg_->norm(&coord_ph[0],B4vectproj);
      }else {
        computeB4vect(B4vect, magneticConfig_, coord_obj, coord_ph);
        BB = sqrt(4.*M_PI*magnetizationParameter_*GYOTO_PROTON_MASS_CGS*GYOTO_C_CGS*GYOTO_C_CGS*number_density);
      }
      theta_mag = get_theta_mag(B4vect, coord_ph, vel);
      Chi=getChi(B4vect, coord_ph, vel); // this is EVPA
    }else{
      theta_mag = 0.; // Average so we don't care
      BB = sqrt(4.*M_PI*magnetizationParameter_*GYOTO_PROTON_MASS_CGS*GYOTO_C_CGS*GYOTO_C_CGS*number_density);
    }
    X[0] = NULL;
    X[1] = NULL;
    X[2] = NULL;
    X[3] = NULL;
    delete[] X;
  
    double nu0 = GYOTO_ELEMENTARY_CHARGE_CGS*BB
      /(2.*M_PI*GYOTO_ELECTRON_MASS_CGS*GYOTO_C_CGS); // cyclotron freq
    double Theta = GYOTO_BOLTZMANN_CGS*temperature
      /(GYOTO_ELECTRON_MASS_CGS*GYOTO_C2_CGS);
    
    //cout << "BB, nu0, Theta : " << BB << ", " << nu0 << ", " << Theta << endl;
    if (number_density!=0. && temperature!=0. && BB!=0.){
      if (emission_=="Thermal"){
        double besselK2 = bessk(2, 1./Theta);
        spectrumThermalSynch_->temperature(temperature);
        spectrumThermalSynch_->numberdensityCGS(number_density);
        spectrumThermalSynch_->angle_averaged(avg); 
        spectrumThermalSynch_->angle_B_pem(theta_mag);
        spectrumThermalSynch_->cyclotron_freq(nu0);
        spectrumThermalSynch_->besselK2(besselK2);
        spectrumThermalSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu, aInu, aQnu, aUnu, aVnu, rotQnu, rotUnu, rotVnu, nuem, nbnu);
      }else if (emission_=="Kappa"){
        double hypergeom = Gyoto::hypergeom(PLindex_+1., Theta);
        spectrumKappaSynch_->kappaindex(PLindex_+1.);
        spectrumKappaSynch_->numberdensityCGS(number_density);
        spectrumKappaSynch_->angle_averaged(avg);
        spectrumKappaSynch_->angle_B_pem(theta_mag);
        spectrumKappaSynch_->cyclotron_freq(nu0);
        spectrumKappaSynch_->thetae(Theta);
        spectrumKappaSynch_->hypergeometric(hypergeom);
        spectrumKappaSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu, aInu, aQnu, aUnu, aVnu, rotQnu, rotUnu, rotVnu, nuem, nbnu);
      }else if (emission_ == "PL"){
        spectrumPLSynch_->numberdensityCGS(number_density);
        spectrumPLSynch_->angle_averaged(avg);
        spectrumPLSynch_->angle_B_pem(theta_mag);
        spectrumPLSynch_->cyclotron_freq(nu0);
        spectrumPLSynch_->PLindex(PLindex_);
        spectrumPLSynch_->gamma_min(gammaMin_);
        spectrumPLSynch_->gamma_max(gammaMax_);
        spectrumPLSynch_->radiativeQ(jInu, jQnu, jUnu, jVnu, aInu, aQnu, aUnu, aVnu, rotQnu, rotUnu, rotVnu, nuem, nbnu);   
      }else if (emission_=="BlackBody"){
        spectrumBB_->temperature(temperature);
        for (size_t ii=0; ii<nbnu; ++ii){
          jInu[ii]=(*spectrumBB_)(nuem[ii]);
          jQnu[ii]=0.;
          jUnu[ii]=0.;
          jVnu[ii]=0.;
          aInu[ii]=0.;
          aQnu[ii]=0.;
          aUnu[ii]=0.;
          aVnu[ii]=0.;
          rotQnu[ii]=0.;
          rotUnu[ii]=0.;
          rotVnu[ii]=0.;
        } 
      }else{
        GYOTO_ERROR("Unknown emission type");
      }
    } else {
      for (size_t ii=0; ii<nbnu; ++ii){
        // Set coefficients to zero to avoid NaN and errors
        jInu[ii]=0.;
        jQnu[ii]=0.;
        jUnu[ii]=0.;
        jVnu[ii]=0.;
        aInu[ii]=0.;
        aQnu[ii]=0.;
        aUnu[ii]=0.;
        aVnu[ii]=0.;
        rotQnu[ii]=0.;
        rotUnu[ii]=0.;
        rotVnu[ii]=0.;
      }
    }
  }else{
    for (size_t ii=0; ii<nbnu; ++ii){
      GYOTO_DEBUG << "Setup arguments for the interpolation of radiative coefficients." << endl;
      double Xq[5] = {tcur, x1, x2, x3, nuem[ii]};
      int X_params[5] = {nfile, nx1_, nx2_, nx3_, nnu_};
      //cout << "coord: " << Xq[0] << " " << Xq[1] << " " << Xq[2] << " " << Xq[3] << " " << Xq[4] << endl;
      double** X;
      X = new double*[5];
      X[0] = time_interpo;
      X[1] = x1_array_;
      X[2] = x2_array_;
      X[3] = x3_array_;
      X[4] = nu_array_;
      std::string cond_limits[5] = {boundCond_[0], boundCond_[1], boundCond_[2], boundCond_[3], boundCond_[4]};
      //cout << "Boundary conditions: " << cond_limits[0] << ", " << cond_limits[1] << ", " << cond_limits[2] << ", " << cond_limits[3] << ", " << cond_limits[4] << endl;
      
      GYOTO_DEBUG << "Perform the interpolation of radiative coefficients." << endl;
      jInu[ii]=interpolate(5, emission_array[0], Xq, X, X_params, cond_limits);
      GYOTO_DEBUG << "Interpolating polarized emission coefficients." << endl;
      jQnu[ii]=emisInFile_[1]?interpolate(5, emission_array[1],   Xq, X, X_params, cond_limits):0.;
      jUnu[ii]=emisInFile_[2]?interpolate(5, emission_array[2],   Xq, X, X_params, cond_limits):0.;
      jVnu[ii]=emisInFile_[3]?interpolate(5, emission_array[3],   Xq, X, X_params, cond_limits):0.;
      GYOTO_DEBUG << "Interpolating polarized absorption coefficients." << endl;
      aInu[ii]=absInFile_[0]?interpolate(5, absorption_array[0], Xq, X, X_params, cond_limits):0.;
      aQnu[ii]=absInFile_[1]?interpolate(5, absorption_array[1], Xq, X, X_params, cond_limits):0.;
      aUnu[ii]=absInFile_[2]?interpolate(5, absorption_array[2], Xq, X, X_params, cond_limits):0.;
      aVnu[ii]=absInFile_[3]?interpolate(5, absorption_array[3], Xq, X, X_params, cond_limits):0.;
      GYOTO_DEBUG << "Interpolating polarized rotation coefficients." << endl;
      rotQnu[ii]=rotInFile_[0]?interpolate(5, rotation_array[0], Xq, X, X_params, cond_limits):0.;
      rotUnu[ii]=rotInFile_[1]?interpolate(5, rotation_array[1], Xq, X, X_params, cond_limits):0.;
      rotVnu[ii]=rotInFile_[2]?interpolate(5, rotation_array[2], Xq, X, X_params, cond_limits):0.;

      GYOTO_DEBUG << "Freeing interpolation arrays." << endl;
      X[0] = NULL;
      X[1] = NULL;
      X[2] = NULL;
      X[3] = NULL;
      X[4] = NULL;
      delete[] X;
    }
  }
  
  GYOTO_DEBUG << "SimBridge::RadiativeQ: Compute radiative increment" << endl;
  for (size_t ii=0; ii<nbnu; ++ii) {
    //cout << "In SimBridge: jInu, jQnu, jUnu, jVnu: " << jInu[ii] << ", " << jQnu[ii] << ", " << jUnu[ii] << ", " << jVnu[ii] << endl;
    //cout << "In SimBridge: aInu, aQnu, aUnu, aVnu: " << aInu[ii] << ", " << aQnu[ii] << ", " << aUnu[ii] << ", " << aVnu[ii] << endl;
    //cout << "In SimBridge: rQnu, rUnu, rVnu: " << rotQnu[ii] << ", " << rotUnu[ii] << ", " << rotVnu[ii] << endl;
    Eigen::Vector4d Jstokes=rotateJs(jInu[ii], jQnu[ii], jUnu[ii], jVnu[ii], Chi);
    Eigen::Matrix4d Omat = Omatrix(aInu[ii], aQnu[ii], aUnu[ii], aVnu[ii], rotQnu[ii], rotUnu[ii], rotVnu[ii], Chi, dsem);
    Eigen::Matrix4d Pmat = Pmatrix(aInu[ii], aQnu[ii], aUnu[ii], aVnu[ii], rotQnu[ii], rotUnu[ii], rotVnu[ii], sin(2.*Chi), cos(2.*Chi), dsem);

    // Computing the increment of the Stokes parameters. Equivalent to dInu=exp(-anu*dsem)*jnu*dsem in the non-polarised case.
    Eigen::Vector4d Stokes=Pmat*Jstokes;

    if (Stokes(0) <=0.){
      Inu[ii] = Qnu[ii] =  Unu[ii] = Vnu[ii] = 0.;
      Onu[ii] = Eigen::Matrix4d::Identity();
    } else {
      Inu[ii] = Stokes(0);
      Qnu[ii] = Stokes(1);
      Unu[ii] = Stokes(2);
      Vnu[ii] = Stokes(3);
      Onu[ii] = Omat;
    }
    
    if (Inu[ii]<0.)
      GYOTO_ERROR("In SimBridge::radiativeQ(): Inu<0");
    if (isnan(Inu[ii]) or isnan(Qnu[ii]) or isnan(Unu[ii]) or isnan(Vnu[ii]) or isnan(Onu[ii](0,0)))
      GYOTO_ERROR("In SimBridge::radiativeQ(): Snu or Taunu is nan");
    if (isinf(Inu[ii]) or isinf(Qnu[ii]) or isinf(Unu[ii]) or isinf(Vnu[ii]) or isinf(Onu[ii](0,0)))
      GYOTO_ERROR("In Simbridge::radiativeQ(): Snu or Taunu is infinite");
  }

  GYOTO_DEBUG << "Freeing interpolation arrays." << endl;
  // Freeing arrays
  if(density_array) delete[] density_array;
  if(temperature_array) delete[] temperature_array;

  if (magneticfield_array) {
    for (int ii = 0; ii < 4; ++ii) {
      if (magneticfield_array[ii]) {
        delete[] magneticfield_array[ii];
      }
    }
    delete[] magneticfield_array;
  }

  if (emission_array) {
    for (int ii = 0; ii < 4; ++ii) {
      if (emission_array[ii]) {
        delete[] emission_array[ii];
      }
    }
    delete[] emission_array;
  }

  if (absorption_array) {
    for (int ii = 0; ii < 4; ++ii) {
      if (absorption_array[ii]) {
        delete[] absorption_array[ii];
      }
    }
    delete[] absorption_array;
  }

  if (rotation_array) {
    for (int ii = 0; ii < 3; ++ii) {
      if (rotation_array[ii]) {
        delete[] rotation_array[ii];
      }
    }
    delete[] rotation_array;
  }
  GYOTO_DEBUG << "Done." << endl;
}

void SimBridge::radiativeQ(double Inu[], double Taunu[], double const nu_em[], size_t nbnu,
			  double dsem, state_t const &coord_ph, double const coord_obj[8]) const{

  double * Qnu           = new double[nbnu];
  double * Unu           = new double[nbnu];
  double * Vnu           = new double[nbnu];
  Eigen::Matrix4d * Onu  = new Eigen::Matrix4d[nbnu];

  const_cast<SimBridge*>(this)->radiativeQ(Inu, Qnu, Unu, Vnu, Onu, nu_em, nbnu, dsem, coord_ph, coord_obj);
  for (size_t ii=0; ii<nbnu; ii++){
    Taunu[ii] = Onu[ii](0,0);
  }
  delete [] Qnu;
  delete [] Unu;
  delete [] Vnu;
  delete [] Onu;
}

void SimBridge::getVelocity(double const pos[4], double vel[4]){
  if (circularmotion_){
    GYOTO_DEBUG << "SimBridge: computing keplerian velocity." << endl;
    double risco = 0.; // Warning: v>c below 1M in Minkowski
    if (gg_->kind()!="Minkowski" && gg_->kind()!="Hayward")
      risco=gg_->getRms();
    // let risco=0 if metric is Minko; then ISCO not defined
    // let also risco=0 for Hayward as we would need to compute it numerically and give it in xml Metric field; not implemented so far
    double rr=pos[1]; // radius
    //cout << "rr: " << rr << ", risco: " << risco << endl;
    if (rr > risco){
      // Keplerian velocity above ISCO
      gg_ -> circularVelocity(pos, vel, 1);
    }else{
      if (cunninghamvel_ && gg_->kind()=="KerrBL"){
        // See formulas in Gralla, Lupsasca & Marrone 2020, Eqs B8-B14
	// initally from Cunnigham 1975
        double SPIN = static_cast<SmartPointer<Metric::KerrBL> >(gg_) -> spin();
        double lambda_ms = (risco*risco - 2.*SPIN*sqrt(risco) + SPIN*SPIN)/(pow(risco,1.5) - 2.*sqrt(risco) + SPIN),
        gamma_ms = sqrt(1.-2./(3.*risco)),
	delta = rr*rr - 2.*rr + SPIN*SPIN,
	hh = (2.*rr - SPIN*lambda_ms)/delta;
	
	vel[0] = gamma_ms*(1.+2./rr*(1.+hh)); // this is: -Ems*g^{tt} + Lms*g^{tp}
	vel[1] = -sqrt(2./(3.*risco))*pow(risco/rr-1.,1.5); // this is: -sqrt{(-1 - g_{tt}*u^t - g_{pp}*u^p - 2*g_{tp}*u^t*u^p)/grr}
	vel[2] = 0.;
        vel[3] = gamma_ms/(rr*rr)*(lambda_ms+SPIN*hh);
      }else{
        cout << "WARNING: Radial motion below the radius of the ISCO for the keplerian rotation! Cunningham prescription only implemented for the Kerr metric." << endl;
        double grr = gg_->gmunu(pos,1,1),
               guptt = gg_->gmunu_up(pos,0,0),
               guptp = gg_->gmunu_up(pos,0,3);
        vel[0] = -guptt;
        vel[1] = -sqrt((-1.-guptt)/grr);
        vel[2] = 0;
        vel[3] = -guptp;
      }
    }
  }else{
    if (dirname_=="None")
        GYOTO_ERROR("In SimBridge RadiativeQ : dirname_ not defined, please use directory(string)");
  
    double tcur=pos[0]; // in M units # TBC
  
    int nfile=ntime_==1?1:2; // 1 if ntime=1, 2 otherwise
    int index=getIndex(tcur); //index of files to be loaded (index and index+1)

    // Creating the velocity arrays
    long ncells = nx1_*nx2_*nx3_; // Number of cells for each time
    double** velocity_array = new double*[3];
    for (int ii=0; ii<3; ii++){
      velocity_array[ii] = new double[nfile*ncells];
    }

    // Reading FITS File
    // double* tmp;
    double time_interpo[nfile];
    for (int ii=0; ii<nfile; ii++){
      ostringstream stream_name ;
      stream_name << dirname_ << fname_ << setw(4) << setfill('0') << index+ii << ".fits" ;
        
      string filename = stream_name.str();
      GYOTO_DEBUG << "Reading FITS file: " << filename << endl ;
    
      fitsfile* fptr = NULL;
    
      fptr = FitsRW::fitsOpen(filename);

      int ntime = FitsRW::fitsReadKey(fptr, "NB_X0");
      int nx1   = FitsRW::fitsReadKey(fptr, "NB_X1");
      int nx2   = FitsRW::fitsReadKey(fptr, "NB_X2");
      int nx3   = FitsRW::fitsReadKey(fptr, "NB_X3");

      if (ntime!=ntime_ || nx1!=nx1_ || nx2!=nx2_ || nx3!=nx3_)
        GYOTO_ERROR("In SimBridge RadiativeQ : size of arrays in FITS file different from initial file");
    
      time_interpo[ii] = FitsRW::fitsReadKey(fptr, "TIME");
    
      double* emplacement = velocity_array[0] + ii * ncells;
      FitsRW::fitsReadHDUData(fptr, "VELOCITY1", emplacement, ncells);
      emplacement = velocity_array[1] + ii * ncells;
      FitsRW::fitsReadHDUData(fptr, "VELOCITY2", emplacement, ncells);
      emplacement = velocity_array[2] + ii * ncells;
      FitsRW::fitsReadHDUData(fptr, "VELOCITY3", emplacement, ncells);

      FitsRW::fitsClose(fptr);
    }

    double Xq[4] = {pos[0], pos[1], pos[2], pos[3]};
    int X_params[4] = {nfile, nx1_, nx2_, nx3_};
    double** X;
    X = new double*[4];
    X[0] = time_interpo;
    X[1] = x1_array_;
    X[2] = x2_array_;
    X[3] = x3_array_;

    vel[0] = 1.;
    vel[1] = interpolate(4, velocity_array[0], Xq, X, X_params, boundCond_);
    vel[2] = interpolate(4, velocity_array[1], Xq, X, X_params, boundCond_);
    vel[3] = interpolate(4, velocity_array[2], Xq, X, X_params, boundCond_);
    //cout << "vel : " << vel[0] << ", "  << vel[1] << ", "  << vel[2] << ", "  << vel[3] << endl;

    X[0] = NULL;
    X[1] = NULL;
    X[2] = NULL;
    X[3] = NULL;
    delete[] X;

    gg_->normalizeFourVel(pos, vel);

    for (int ii=0; ii<3; ii++){
      if (velocity_array[ii]){
        delete [] velocity_array[ii];
        velocity_array[ii] = NULL;
      }
    }
  }
  GYOTO_DEBUG << "SimBridge: Done." <<endl;
}
  

double SimBridge::operator()(double const coord[4]) {
  double rr=0.; // radius
  double distance;
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    {
    rr = coord[1];
    break;
    }
  case GYOTO_COORDKIND_CARTESIAN:
    {
    rr = pow(coord[1]*coord[1]+coord[2]*coord[2]+coord[3]*coord[3], 0.5);
    break;
    }
  default:
    {
    GYOTO_ERROR("In SimBridge::operator(): Unknown coordinate system kind");
    }
  }
  distance = rr - rmax_;
  return distance;
}

int SimBridge::getIndex(double const tcur) const {
  if (ntime_==1){
    return 0;
  }else{
    if (tcur<time_array_[0])
      return 0; // Set index at 0, the output of the interpolation will be defined by the boundary conditions
    else if (tcur>time_array_[ntime_-1])
       return ntime_-2; // Set index at ntime_-1, the output of the interpolation will be defined by the boundary conditions
    else{
      int ii=0;
      while (tcur > time_array_[ii] and ii<ntime_-1){
        ii+=1;
      }
      return ii - 1;
    }
  }
}
