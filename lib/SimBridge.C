/*
    Copyright 2024 Aimar Nicolas

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
GYOTO_PROPERTY_DOUBLE(SimBridge, GammaMin, gammaMin)
GYOTO_PROPERTY_DOUBLE(SimBridge, GammaMax, gammaMax)
GYOTO_PROPERTY_BOOL(SimBridge, TemperatureGrid, IntensityGrid, temperature)
GYOTO_PROPERTY_DOUBLE(SimBridge, PLindex, PLindex)
GYOTO_PROPERTY_DOUBLE(SimBridge, FloorTemperature, floorTemperature)
GYOTO_PROPERTY_STRING(SimBridge, EmissionType, emissionType)
GYOTO_PROPERTY_END(SimBridge, Standard::properties)


SimBridge::SimBridge() :
  Standard("SimBridge"),
  FitsRW(),
  dirname_("None"),
  fname_("data"),
  magneticConfig_("None"),
  magnetizationParameter_(1.),
  emission_("None"),
  gammaMin_(1.),
  gammaMax_(1.),
  PLindex_(1.),
  floortemperature_(0.),
  temperature_(true),
  openingAngle_(M_PI/2),
  ntime_(1),
  nx1_(1),
  nx2_(1),
  nx3_(1),
  nnu_(0),
  time_array_(NULL),
  x1_array_(NULL),
  x2_array_(NULL),
  x3_array_(NULL),
  nu_array_(NULL),
  boundCond_(NULL),
  spectrumKappaSynch_(NULL),
  spectrumPLSynch_(NULL),
  spectrumThermalSynch_(NULL)
{
  boundCond_ = new string[5];
  boundCond_[0]=boundCond_[1]=boundCond_[2]=boundCond_[3]=boundCond_[4]="None";

  spectrumKappaSynch_ = new Spectrum::KappaDistributionSynchrotron();
  spectrumPLSynch_ = new Spectrum::PowerLawSynchrotron();
  spectrumThermalSynch_ = new Spectrum::ThermalSynchrotron();
}

SimBridge::SimBridge(const SimBridge& orig) :
  Standard(orig),
  FitsRW(orig),
  dirname_(orig.dirname_),
  fname_(orig.fname_),
  magneticConfig_(orig.magneticConfig_),
  magnetizationParameter_(orig.magnetizationParameter_),
  emission_(orig.emission_),
  gammaMin_(orig.gammaMin_),
  gammaMax_(orig.gammaMax_),
  PLindex_(orig.PLindex_),
  floortemperature_(orig.floortemperature_),
  temperature_(orig.temperature_),
  openingAngle_(orig.openingAngle_),
  ntime_(orig.ntime_),
  nx1_(orig.nx1_),
  nx2_(orig.nx2_),
  nx3_(orig.nx3_),
  nnu_(orig.nnu_),
  time_array_(NULL),
  x1_array_(NULL),
  x2_array_(NULL),
  x3_array_(NULL),
  nu_array_(NULL),
  boundCond_(NULL),
  spectrumKappaSynch_(NULL),
  spectrumPLSynch_(NULL),
  spectrumThermalSynch_(NULL)
{
  if (orig.boundCond_){
    boundCond_  = new string[5];
    for (int ii=0; ii<5; ii++){
      boundCond_[ii] = orig.boundCond_[ii];
    }
  }
  if (orig.time_array_){
    time_array_ = new double[ntime_];
    memcpy(orig.time_array_, time_array_, ntime_*sizeof(double));
  }

  if (orig.x1_array_){
    x1_array_ = new double[nx1_];
    memcpy(orig.x1_array_, x1_array_, nx1_*sizeof(double));
  }

  if (orig.x2_array_){
    x2_array_ = new double[nx2_];
    memcpy(orig.x2_array_, x2_array_, nx2_*sizeof(double));
  }

  if (orig.x3_array_){
    x3_array_ = new double[nx3_];
    memcpy(orig.x3_array_, x3_array_, nx3_*sizeof(double));
  }
  
  if (orig.nu_array_){
    nu_array_ = new double[nnu_];
    memcpy(orig.nu_array_, nu_array_, nnu_*sizeof(double));
  }

  if (orig.spectrumKappaSynch_()) spectrumKappaSynch_=orig.spectrumKappaSynch_->clone();
  if (orig.spectrumPLSynch_()) spectrumPLSynch_=orig.spectrumPLSynch_->clone();
  if (orig.spectrumThermalSynch_()) spectrumThermalSynch_=orig.spectrumThermalSynch_->clone();
}

SimBridge* SimBridge::clone() const { return new SimBridge(*this); }

SimBridge::~SimBridge() {
  if (debug()) cerr << "DEBUG: SimBridge::~SimBridge()\n";
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

void SimBridge::filename(std::string f){
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

  ntime_        = FitsRW::fitsReadKey(fptr, "NB_X0");
  time_array_   = new double[ntime_];
  time_array_   = FitsRW::fitsReadHDUData(fptr, "X0");

  nx1_          = FitsRW::fitsReadKey(fptr, "NB_X1");
  x1_array_     = new double[nx1_];
  x1_array_     = FitsRW::fitsReadHDUData(fptr, "X1");

  nx2_          = FitsRW::fitsReadKey(fptr, "NB_X2");
  x2_array_     = new double[nx2_];
  x2_array_     = FitsRW::fitsReadHDUData(fptr, "X2");

  nx3_         = FitsRW::fitsReadKey(fptr, "NB_X3");
  x3_array_    = new double[nx3_];
  x3_array_    = FitsRW::fitsReadHDUData(fptr, "X3");

  try{
    nnu_       = FitsRW::fitsReadKey(fptr, "NB_FREQ");
    nu_array_  = new double[nnu_];
    nu_array_  = FitsRW::fitsReadHDUData(fptr, "FREQ");
  }catch(...){
    GYOTO_DEBUG << "No frequency founded in FITS file." << endl;
  }
  FitsRW::fitsClose(fptr);

  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    {
    critical_value_ = x1_array_[nx1_-1];
    break;
    }
  case GYOTO_COORDKIND_CARTESIAN:
    {
    double xmax = x1_array_[nx1_-1];
    double ymax = x2_array_[nx2_-1];
    double zmax = x3_array_[nx2_-1];
    critical_value_ = 50.; //pow(xmax*xmax+ymax*ymax+zmax*zmax, 0.5);
    break;
    }
  default:
    {
    GYOTO_ERROR("In SimBridge::filename : Unknown coordinate system kind");
    }
  }
  safety_value_ = critical_value_*1.1;
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

void SimBridge::temperature(bool t){
  temperature_=t;
}
bool SimBridge::temperature() const{
  return temperature_;
}

void SimBridge::floorTemperature(double t){
  floortemperature_=t;
}
double SimBridge::floorTemperature()const{
  return floortemperature_;
}

void SimBridge::magneticConfiguration(string config){
  magneticConfig_=config;
}
string SimBridge::magneticConfiguration() const{
  return magneticConfig_;
}

void SimBridge::magnetization(double ss){
  magnetizationParameter_=ss;
}

double SimBridge::magnetization() const{
  return magnetizationParameter_;
}

void SimBridge::boundaryConditions(std::string x0BC, std::string x1BC, std::string x2BC, std::string x3BC, std::string freqBC){
  boundCond_[0] = x0BC;
  boundCond_[1] = x1BC;
  boundCond_[2] = x2BC;
  boundCond_[3] = x3BC;
  boundCond_[4] = freqBC;
}
std::string * SimBridge::boundaryConditions() const{
  return boundCond_;
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
    throwError("unknown electron distribution!");
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
  
  int index; //index of files to be loaded (index and index+1)
  if (ntime_==1){
    index=0;
  }else{
    if (tcur<time_array_[0])
      index=0; // Set index at 0, the output of the interpolation will be defined by the boundary conditions
    if (tcur>time_array_[ntime_-1])
      index=ntime_-2; // Set index at ntime_-1, the output of the interpolation will be defined by the boundary conditions
  }

  // Creating arrays
  long ncells = nx1_*nx2_*nx3_; // Number of cells for each time
  double* density_array       = new double[2*ncells];
  double* temperature_array   = new double[2*ncells];
  
  double** magneticfield_array = new double*[4];
  for (int ii=0; ii<4; ii++){
    magneticfield_array[ii] = new double[2*ncells];
  }

  double** emission_array      = new double*[4];
  double** absorption_array    = new double*[4];
  double** rotation_array      = new double*[3];
  for (int ii=0; ii<4; ii++){
    emission_array[ii]   = new double[2*ncells];
    absorption_array[ii] = new double[2*ncells];
  }
  for (int ii=0; ii<3; ii++){
    rotation_array[ii] = new double[2*ncells];
  }

  // Opening and reading Files
  double* tmp;
  double time_interpo[2];
  for (int ii=0; ii<2; ii++){
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

    if (temperature_){
      // The arrays stored in the FITS files are density and temperature + eventually magnetic field
      // Read NUMBERDENSITY
      tmp = FitsRW::fitsReadHDUData(fptr, "NUMBERDENSITY");
      std::memcpy(density_array + ii * ncells, tmp, ncells * sizeof(double));
      delete[] tmp;

      // Read TEMPERATURE
      tmp = FitsRW::fitsReadHDUData(fptr, "TEMPERATURE");
      std::memcpy(temperature_array + ii * ncells, tmp, ncells * sizeof(double));
      delete[] tmp;
      try {
        for (int jj = 0; jj < 4; jj++) {
          std::ostringstream hdu_name;
          hdu_name << "B" << jj;
          tmp = FitsRW::fitsReadHDUData(fptr, hdu_name.str().c_str());
          std::memcpy(magneticfield_array[jj] + ii * ncells, tmp, ncells * sizeof(double));
          delete[] tmp;
        }
      } catch (...) {
        // Free magnetic field arrays if not found in FITS files and if ii == 0
        if (magneticfield_array && ii == 0) {
          for (int jj = 0; jj < 4; jj++) {
            if (magneticfield_array[jj]) {
              delete[] magneticfield_array[jj];
              magneticfield_array[jj] = nullptr;
            }
          }
          delete[] magneticfield_array;
          magneticfield_array = nullptr;
        }
      }
    }else{
      // The arrays stored in the FITS files are (at least) emission(s).
      if (emission_array==NULL or absorption_array==NULL or rotation_array==NULL)
        GYOTO_ERROR("In SimBridge RadiativeQ : emission, absorption and/or rotation arrays are not properly set.");
      if (nu_array_==NULL)
        GYOTO_ERROR("In SimBridge RadiativeQ : frequency array is not set");
      
      tmp = FitsRW::fitsReadHDUData(fptr, "J_I");
      memcpy(emission_array[0]+ii*ncells, tmp, ncells*sizeof(double));
      delete[] tmp;
      try {
        tmp = FitsRW::fitsReadHDUData(fptr, "J_Q");
        memcpy(emission_array[1]+ii*ncells, tmp, ncells*sizeof(double));
        delete[] tmp;
        tmp = FitsRW::fitsReadHDUData(fptr, "J_U");
        memcpy(emission_array[2]+ii*ncells, tmp, ncells*sizeof(double));
        delete[] tmp;
        tmp = FitsRW::fitsReadHDUData(fptr, "J_V");
        memcpy(emission_array[3]+ii*ncells, tmp, ncells*sizeof(double));
        delete[] tmp;
      } catch (...) {
        memset(emission_array[1], 0, 2*ncells*sizeof(double));
        memset(emission_array[2], 0, 2*ncells*sizeof(double));
        memset(emission_array[3], 0, 2*ncells*sizeof(double));
      }
      try{
        tmp = FitsRW::fitsReadHDUData(fptr, "ALPHA_I");
        memcpy(absorption_array[0]+ii*ncells, tmp, ncells*sizeof(double));
        delete[] tmp;
      } catch (...) {
        memset(absorption_array[0], 0, 2*ncells*sizeof(double));
      }
      try{
        tmp = FitsRW::fitsReadHDUData(fptr, "ALPHA_Q");
        memcpy(absorption_array[1]+ii*ncells, tmp, ncells*sizeof(double));
        delete[] tmp;
        tmp = FitsRW::fitsReadHDUData(fptr, "ALPHA_U");
        memcpy(absorption_array[2]+ii*ncells, tmp, ncells*sizeof(double));
        delete[] tmp;
        tmp = FitsRW::fitsReadHDUData(fptr, "ALPHA_V");
        memcpy(absorption_array[3]+ii*ncells, tmp, ncells*sizeof(double));
        delete[] tmp;
        tmp = FitsRW::fitsReadHDUData(fptr, "R_Q");
        memcpy(rotation_array[1]+ii*ncells, tmp, ncells*sizeof(double));
        delete[] tmp;
        tmp = FitsRW::fitsReadHDUData(fptr, "R_U");
        memcpy(rotation_array[2]+ii*ncells, tmp, ncells*sizeof(double));
        delete[] tmp;
        tmp = FitsRW::fitsReadHDUData(fptr, "R_V");
        memcpy(rotation_array[3]+ii*ncells, tmp, ncells*sizeof(double));
        delete[] tmp;
      } catch (...) {
        memset(absorption_array[1], 0, 2*ncells*sizeof(double));
        memset(absorption_array[2], 0, 2*ncells*sizeof(double));
        memset(absorption_array[3], 0, 2*ncells*sizeof(double));
        memset(rotation_array[0],   0, 2*ncells*sizeof(double));
        memset(rotation_array[1],   0, 2*ncells*sizeof(double));
        memset(rotation_array[2],   0, 2*ncells*sizeof(double));
      }      
    }
    FitsRW::fitsClose(fptr);
  } // End of reading FITS files. At this point, the relevent arrays should be filled.

  double x1 = coord_ph[1], x2 = coord_ph[2], x3 = coord_ph[3];
  double pos[4] = {tcur, x1, x2, x3};
  double vel[4];
  const_cast<SimBridge*>(this)->getVelocity(pos, vel);

  // Defining emission, absoprtion and rotation coefficients for the transmission matrix
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
    int X_params[4] = {2, nx1_, nx2_, nx3_};
    double** X;
    X = new double*[4];
    X[0] = time_interpo;
    X[1] = x1_array_;
    X[2] = x2_array_;
    X[3] = x3_array_;

    double number_density = interpolate(4, density_array, Xq, X, X_params, boundCond_);
    double temperature    = max(interpolate(4, temperature_array, Xq, X, X_params, boundCond_),floortemperature_);
    //cout << "ne, Te at (t,r,theta, phi) : "  << number_density << ", " << temperature << ", (" << tcur << "," << x1 << "," << x2 << "," << x3  << ")" << endl;
    delete[] X;

    int avg=0; // flag for magnetic field average for synchrotron
    if (magneticfield_array==NULL && magneticConfig_=="None")
      avg = 1;
    
    double B4vect[4]={0.,0.,0.,0.};
    double theta_mag, BB;
    if (!avg){
      if (magneticfield_array){
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
      double Xq[5] = {tcur, x1, x2, x3, nuem[ii]};
      int X_params[5] = {2, nx1_, nx2_, nx3_, nnu_};
      double** X;
      X = new double*[5];
      X[0] = time_interpo;
      X[1] = x1_array_;
      X[2] = x2_array_;
      X[3] = x3_array_;
      X[4] = nu_array_;
      std::string cond_limits[5] = {boundCond_[0], boundCond_[1], boundCond_[2], boundCond_[3], "None"};

      jInu[ii]=interpolate(5, emission_array[0],   Xq, X, X_params, cond_limits);
      jQnu[ii]=interpolate(5, emission_array[1],   Xq, X, X_params, cond_limits);
      jUnu[ii]=interpolate(5, emission_array[2],   Xq, X, X_params, cond_limits);
      jVnu[ii]=interpolate(5, emission_array[3],   Xq, X, X_params, cond_limits);
      aInu[ii]=interpolate(5, absorption_array[0], Xq, X, X_params, cond_limits);
      aQnu[ii]=interpolate(5, absorption_array[1], Xq, X, X_params, cond_limits);
      aUnu[ii]=interpolate(5, absorption_array[2], Xq, X, X_params, cond_limits);
      aVnu[ii]=interpolate(5, absorption_array[3], Xq, X, X_params, cond_limits);
      rotQnu[ii]=interpolate(5, rotation_array[0], Xq, X, X_params, cond_limits);
      rotUnu[ii]=interpolate(5, rotation_array[1], Xq, X, X_params, cond_limits);
      rotVnu[ii]=interpolate(5, rotation_array[2], Xq, X, X_params, cond_limits);

      delete[] X;
    }
  }
  
  for (size_t ii=0; ii<nbnu; ++ii) {
    //cout << "In SimBridge: jInu, jQnu, jUnu, jVnu: " << jInu[ii] << ", " << jQnu[ii] << ", " << jUnu[ii] << ", " << jVnu[ii] << endl;
    //cout << "In SimBridge: aInu, aQnu, aUnu, aVnu: " << aInu[ii] << ", " << aQnu[ii] << ", " << aUnu[ii] << ", " << aVnu[ii] << endl;
    //cout << "In SimBridge: rQnu, rUnu, rVnu: " << rotQnu[ii] << ", " << rotUnu[ii] << ", " << rotVnu[ii] << endl;
    Eigen::Vector4d Jstokes=rotateJs(jInu[ii], jQnu[ii], jUnu[ii], jVnu[ii], Chi)*dsem*gg_->unitLength();
    Eigen::Matrix4d Omat = Omatrix(aInu[ii], aQnu[ii], aUnu[ii], aVnu[ii], rotQnu[ii], rotUnu[ii], rotVnu[ii], Chi, dsem);
    Eigen::Vector4d Stokes=Omat*Jstokes;
    Inu[ii] = Stokes(0);
    Qnu[ii] = Stokes(1);
    Unu[ii] = Stokes(2);
    Vnu[ii] = Stokes(3);
    Onu[ii] = Omat;

    if (Inu[ii]<0.)
      GYOTO_ERROR("In Plasmoid::radiativeQ(): Inu<0");
    if (Inu[ii]!=Inu[ii] or Onu[ii](0,0)!=Onu[ii](0,0))
      GYOTO_ERROR("In Plasmoid::radiativeQ(): Inu or Taunu is nan");
    if (Inu[ii]==Inu[ii]+1. or Onu[ii](0,0)==Onu[ii](0,0)+1.)
      GYOTO_ERROR("In Plasmoid::radiativeQ(): Inu or Taunu is infinite");
  }


  // Freeing arrays
  delete[] density_array;
  delete[] temperature_array;

  for (int ii = 0; ii < 4; ++ii) {
    if (magneticfield_array) {
      delete[] magneticfield_array[ii];
    }
    delete[] emission_array[ii];
    delete[] absorption_array[ii];
  }
  if (magneticfield_array) {
    delete[] magneticfield_array;
  }
  delete[] emission_array;
  delete[] absorption_array;
  for (int ii = 0; ii < 3; ++ii) {
    delete[] rotation_array[ii];
  }
  delete[] rotation_array;
}

void SimBridge::radiativeQ(double Inu[], double Taunu[], double const nu_em[], size_t nbnu,
			  double dsem, state_t const &coord_ph, double const coord_obj[8]) const{

  double * Qnu           = new double[nbnu];
  double * Unu           = new double[nbnu];
  double * Vnu           = new double[nbnu];
  Eigen::Matrix4d * Onu  = new Eigen::Matrix4d[nbnu];

  const_cast<SimBridge*>(this)->radiativeQ(Inu, Qnu, Unu, Vnu, Onu, nu_em, nbnu, dsem, coord_ph, coord_obj);
  for (int ii=0; ii<nbnu; ii++){
    Taunu[ii] = Onu[ii](0,0);
  }
  delete [] Qnu;
  delete [] Unu;
  delete [] Vnu;
  delete [] Onu;
}

void SimBridge::getVelocity(double const pos[4], double vel[4]){
  if (dirname_=="None")
      GYOTO_ERROR("In SimBridge RadiativeQ : dirname_ not defined, please use directory(string)");
  
  double tcur=pos[0]; // in M units # TBC
  
  int index; //index of files to be loaded (index and index+1)
  if (ntime_==1){
    index=0;
  }else{
    if (tcur<time_array_[0])
      index=0; // Set index at 0, the output of the interpolation will be defined by the boundary conditions
    if (tcur>time_array_[ntime_-1])
      index=ntime_-2; // Set index at ntime_-1, the output of the interpolation will be defined by the boundary conditions
  }

  // Creating the velocity arrays
  long ncells = nx1_*nx2_*nx3_; // Number of cells for each time
  double** velocity_array = new double*[4];
  for (int ii=0; ii<4; ii++){
    velocity_array[ii] = new double[2*ncells];
  }

  // Reading FITS File
  double* tmp;
  bool vel4=true;
  double time_interpo[2];
  for (int ii=0; ii<2; ii++){
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
    try{
      tmp = FitsRW::fitsReadHDUData(fptr, "VELOCITY0");
      memcpy(velocity_array[0]+ii*ncells, tmp, ncells*sizeof(double));
      delete[] tmp;
    } catch (...) {
      vel4 = false;
    }
    tmp = FitsRW::fitsReadHDUData(fptr, "VELOCITY1");
    memcpy(velocity_array[1]+ii*ncells, tmp, ncells*sizeof(double));
    delete[] tmp;
    tmp = FitsRW::fitsReadHDUData(fptr, "VELOCITY2");
    memcpy(velocity_array[2]+ii*ncells, tmp, ncells*sizeof(double));
    delete[] tmp;
    tmp = FitsRW::fitsReadHDUData(fptr, "VELOCITY3");
    memcpy(velocity_array[3]+ii*ncells, tmp, ncells*sizeof(double));
    delete[] tmp;

    FitsRW::fitsClose(fptr);
  }

  double Xq[4] = {pos[0], pos[1], pos[2], pos[3]};
  int X_params[4] = {2, nx1_, nx2_, nx3_};
  double** X;
  X = new double*[4];
  X[0] = time_interpo;
  X[1] = x1_array_;
  X[2] = x2_array_;
  X[3] = x3_array_;

  if (vel4)
    vel[0] = interpolate(4, velocity_array[0], Xq, X, X_params, boundCond_);
  else
    vel[0] = 1;
  vel[1] = interpolate(4, velocity_array[1], Xq, X, X_params, boundCond_);
  vel[2] = interpolate(4, velocity_array[2], Xq, X, X_params, boundCond_);
  vel[3] = interpolate(4, velocity_array[3], Xq, X, X_params, boundCond_);

  delete[] X;

  if (!vel4)
    gg_->normalizeFourVel(pos, vel);

  for (int ii=0; ii<4; ii++){
    if (velocity_array[ii]){
      delete [] velocity_array[ii];
      velocity_array[ii] = NULL;
    }
  }
}

double SimBridge::operator()(double const coord[4]) {
  /*
  double rr=0., rmax=0.; // radius
  switch (gg_->coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    {
    rr = coord[1];
    rmax = x1_array_[nx1_-1];
    break;
    }
  case GYOTO_COORDKIND_CARTESIAN:
    {
    rr = pow(coord[1]*coord[1]+coord[2]*coord[2]+coord[3]*coord[3], 0.5);
    double xmax = x1_array_[nx1_-1];
    double ymax = x2_array_[nx2_-1];
    double zmax = x3_array_[nx2_-1];
    rmax = pow(xmax*xmax+ymax*ymax+zmax*zmax, 0.5);
    break;
    }
  default:
    {
    GYOTO_ERROR("In SimBridge::operator(): Unknown coordinate system kind");
    }
  }
  
  double distance = rr - rmax;
  return distance;*/

  // zpos: modulus of altitude above equatorial plane
  // rproj: radius projected in the equatorial plane
  double zpos=0., rproj=0.;
  
  switch (gg_ -> coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    rproj  = coord[1]*sin(coord[2]);
    zpos  = fabs(coord[1]*cos(coord[2]));
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    zpos  = fabs(coord[3]);
    rproj  = sqrt(coord[1]*coord[1]+coord[2]*coord[2]);
    break;
  default:
    GYOTO_ERROR("SimBridge::operator(): unknown COORDKIND");
  }
  double zdisk = 0.;   // zdisk is fixed at zero rproj <= rinner,
  // then the distance to the disk is always positive
  double rproj_lim=x1_array_[0];
  if (rproj > rproj_lim) // usual linear surface above rproj_lim
    zdisk = (rproj - rproj_lim) * tan(M_PI/2. - openingAngle_) ; 
  return zpos - zdisk; // >0 outside, <0 inside flared disk
}

void SimBridge::openingAngle(double angle){
  openingAngle_ = angle;
}

double SimBridge::openingAngle() const{
  return openingAngle_;
}