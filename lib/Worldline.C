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

#include <GyotoWorldline.h>
#include <GyotoUtils.h>
#include <GyotoProperty.h>
#include <GyotoFactoryMessenger.h>
#include <iostream>
#include <sstream>
#include <cmath>
#include <cfloat>
#include <cstdlib>
#include <cstring>
#include <iomanip>
using namespace std;
using namespace Gyoto;

GYOTO_PROPERTY_BOOL(Worldline,
		    HighOrderImages, PrimaryOnly, secondary,
		    Object::properties);
GYOTO_PROPERTY_DOUBLE(Worldline, RelTol, relTol, &HighOrderImages);
GYOTO_PROPERTY_DOUBLE(Worldline, AbsTol, absTol, &RelTol);
GYOTO_PROPERTY_DOUBLE(Worldline, DeltaMaxOverR, deltaMaxOverR, &AbsTol);
GYOTO_PROPERTY_DOUBLE(Worldline, DeltaMax, deltaMax, &DeltaMaxOverR);
GYOTO_PROPERTY_DOUBLE(Worldline, DeltaMin, deltaMin, &DeltaMax);
GYOTO_PROPERTY_STRING(Worldline, Integrator, integrator, &DeltaMin);
GYOTO_PROPERTY_SIZE_T(Worldline, MaxIter, maxiter, &Integrator);
GYOTO_PROPERTY_BOOL(Worldline, Adaptive, NonAdaptive, adaptive, &MaxIter);
GYOTO_PROPERTY_DOUBLE_UNIT(Worldline, Delta, delta, &Adaptive);
GYOTO_PROPERTY_VECTOR_DOUBLE(Worldline, InitCoord, initCoord, &Delta);
GYOTO_PROPERTY_METRIC(Worldline, Metric, metric, &InitCoord);
GYOTO_PROPERTY_FINALIZE(Worldline, &::Metric);

Worldline::Worldline() : stopcond(0), metric_(NULL),
                         imin_(1), i0_(0), imax_(0), adaptive_(1),
			 secondary_(1),
			 delta_(GYOTO_DEFAULT_DELTA),
			 tmin_(-DBL_MAX), cst_(NULL), cst_n_(0),
			 wait_pos_(0), init_vel_(NULL),
			 maxiter_(GYOTO_DEFAULT_MAXITER),
			 delta_min_(GYOTO_DEFAULT_DELTA_MIN),
			 delta_max_(GYOTO_DEFAULT_DELTA_MAX),
			 delta_max_over_r_(GYOTO_DEFAULT_DELTA_MAX_OVER_R),
			 abstol_(GYOTO_DEFAULT_ABSTOL),
			 reltol_(GYOTO_DEFAULT_RELTOL)
{ 
  xAllocate();

#ifdef HAVE_BOOST
  state_ = new Worldline::IntegState::Boost(this, "runge_kutta_fehlberg78");
#else
  state_ = new Worldline::IntegState::Legacy(this);
#endif
  state_ -> init();
}

Worldline::Worldline(const Worldline& orig) :
  metric_(NULL),
  x_size_(orig.x_size_), imin_(orig.imin_), i0_(orig.i0_), imax_(orig.imax_),
  adaptive_(orig.adaptive_), secondary_(orig.secondary_),
  delta_(orig.delta_), tmin_(orig.tmin_), cst_(NULL), cst_n_(orig.cst_n_),
  wait_pos_(orig.wait_pos_), init_vel_(NULL),
  maxiter_(orig.maxiter_),
  delta_min_(orig.delta_min_),
  delta_max_(orig.delta_max_),
  delta_max_over_r_(orig.delta_max_over_r_),
  abstol_(orig.abstol_),
  reltol_(orig.reltol_),
  state_(NULL)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  if (orig.metric_()) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "cloning metric\n";
#   endif
    metric_=orig.metric_->clone();
  }

  state_ = orig.state_->clone(this);

  xAllocate(x_size_);
  size_t sz = get_nelements()*sizeof(double);
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "sz="<<sz<<", imin_="<<imin_<<endl;
# endif
  memcpy(x0_+imin_, orig.x0_+imin_, sz);
  memcpy(x1_+imin_, orig.x1_+imin_, sz);
  memcpy(x2_+imin_, orig.x2_+imin_, sz);
  memcpy(x3_+imin_, orig.x3_+imin_, sz);
  memcpy(x0dot_+imin_, orig.x0dot_+imin_, sz);
  memcpy(x1dot_+imin_, orig.x1dot_+imin_, sz);
  memcpy(x2dot_+imin_, orig.x2dot_+imin_, sz);
  memcpy(x3dot_+imin_, orig.x3dot_+imin_, sz);
  if (orig.cst_ && cst_n_) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "cloning constants of motion\n";
#   endif
    cst_ = new double [cst_n_];
    memcpy(cst_, orig.cst_, cst_n_*sizeof(double));
  }
  if (orig.init_vel_) {
    init_vel_ = new double [3];
    memcpy(init_vel_, orig.init_vel_, 3*sizeof(double));
  }
  state_ -> init();
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done\n";
# endif
}

Worldline::Worldline(Worldline *orig, size_t i0, int dir, double step_max) :
  metric_(orig->metric_),
//  x_size_(orig.x_size_), imin_(orig.imin_), i0_(orig.i0_), imax_(orig.imax_),
  adaptive_(orig->adaptive_), secondary_(orig->secondary_),
  delta_(orig->delta_), tmin_(orig->tmin_), cst_n_(orig->cst_n_),
  wait_pos_(orig->wait_pos_), init_vel_(NULL),
  maxiter_(orig->maxiter_),
  delta_min_(orig->delta_min_),
  delta_max_(orig->delta_max_),
  delta_max_over_r_(orig->delta_max_over_r_),
  abstol_(orig->abstol_),
  reltol_(orig->reltol_)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  state_ = orig->state_->clone(this);
  double d1 = orig->x0_[i0], d2 = orig->x0_[i0+dir];
  x_size_= size_t(fabs(d1-d2)/step_max)+2;
  double step = (d2-d1)/double(x_size_-1);
  xAllocate(x_size_);
  imin_=0; imax_=x_size_-1; i0_=(dir==1?imin_:imax_);
  x0_[i0_]=d1;
  size_t i=i0_;
  for (i=i0_+dir; i>imin_ && i<imax_; i+=dir) x0_[i] = x0_[i-dir]+step;
  x0_[i]=d2;

  orig->getCoord(x0_, x_size_, x1_, x2_, x3_, x0dot_, x1dot_, x2dot_, x3dot_);

# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
    {
      GYOTO_DEBUG << "(Worldline*, "<<i0<<", "<<dir<<", "<<step_max<<")"<<endl;
      GYOTO_DEBUG << "d1="<<d1<<", d2="<<d2<<endl;
      GYOTO_DEBUG_ARRAY(x0_, x_size_);
      GYOTO_DEBUG_ARRAY(x1_, x_size_);
      GYOTO_DEBUG_ARRAY(x2_, x_size_);
      GYOTO_DEBUG_ARRAY(x3_, x_size_);
      GYOTO_DEBUG_ARRAY(x0dot_, x_size_);
      GYOTO_DEBUG_ARRAY(x1dot_, x_size_);
      GYOTO_DEBUG_ARRAY(x2dot_, x_size_);
      GYOTO_DEBUG_ARRAY(x3dot_, x_size_);
    }
  GYOTO_ENDIF_DEBUG
# endif

  if (orig->cst_ && cst_n_) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG << "cloning constants of motion\n";
#   endif
    cst_ = new double [cst_n_];
    memcpy(cst_, orig->cst_, cst_n_*sizeof(double));
  }
  /*
  if (orig->init_vel_) {
    init_vel_ = new double [3];
    memcpy(init_vel_, orig->init_vel_, 3*sizeof(double));
  }
  */
  state_ -> init();
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done\n";
# endif
}

Worldline::~Worldline(){
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  if (metric_) metric_ -> unhook(this);
  delete[] x0_;
  delete[] x1_;
  delete[] x2_;
  delete[] x3_;
  delete[] x0dot_;
  delete[] x1dot_;
  delete[] x2dot_;
  delete[] x3dot_;
  if (cst_) delete [] cst_;
  if (init_vel_) delete[] init_vel_;
  state_=NULL;
}
void Worldline::xAllocate() {xAllocate(GYOTO_DEFAULT_X_SIZE);}

void Worldline::xAllocate(size_t sz)
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(sz);
# endif
  x_size_ = sz ;
  x0_ = new double[x_size_];
  x1_ = new double[x_size_];
  x2_ = new double[x_size_];
  x3_ = new double[x_size_];
  x0dot_ = new double[x_size_];
  x1dot_ = new double[x_size_];
  x2dot_ = new double[x_size_];
  x3dot_ = new double[x_size_];
}

void Worldline::xExpand(double* &x, int dir) {
  double * old;
  size_t offset=(dir==1)?0:x_size_;
  size_t i;
  size_t nsize=2*x_size_;

  old=x;
  x=new double[nsize];
  for (i=imin_;i<=imax_;++i) x[i+offset]=old[i];
  delete [] old;
}

size_t Worldline::xExpand(int dir) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(dir);
# endif

  xExpand(x0_, dir);
  xExpand(x1_, dir);
  xExpand(x2_, dir);
  xExpand(x3_, dir);
  xExpand(x0dot_, dir);
  xExpand(x1dot_, dir);
  xExpand(x2dot_, dir);
  xExpand(x3dot_, dir);

  size_t retval=(dir==1)?(x_size_-1):x_size_;
  size_t offset=(dir==1)?0:x_size_;
  x_size_*=2;
  
  imin_+=offset;
  i0_+=offset;
  imax_+=offset;

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG<< ", xsize_=" << x_size_
		    << ", imin_=" << imin_
		    << ", i0_=" << i0_
		    << ", imax_=" << imax_;
# endif

  return retval;
}

void Worldline::metric(SmartPointer<Metric::Generic> gg) {
  // Unhook from previous metric
  if (metric_) metric_ -> unhook(this);

  // Set the Metric
  metric_=gg;

  // Hook to new metric
  if (metric_) metric_ -> hook(this);

  // Reinit integrator: some cache the metric
  state_->init();

  // Reset integration
  reInit();

}

void Worldline::tell(Gyoto::Hook::Teller* msg) {
  if (msg != metric_) {
    throwError("Worldline::tell(): wrong Teller");
  }
  reInit();
}

#ifdef GYOTO_USE_XERCES

void Worldline::fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const {
  if (p.name == "InitCoord") {
    if (imin_ <= imax_) {
      double coord[8];
      getInitialCoord(coord);
      if (getMass()) {
	// For massive particule, express initial condition with 3-velocity
	double vel[3] = {coord[5]/coord[4], coord[6]/coord[4], coord[7]/coord[4]};
	fmp -> setParameter ("Position", coord, 4);
	fmp -> setParameter ("Velocity", vel, 3);
      } else {
	// For massless particle, only 4-velocity is meaningfull
	fmp -> setParameter("InitCoord", coord, 8);
      }
    }
    Property const * const * parent = p.parents;
    if (parent) {
      for ( ; *parent; ++parent) {
	fillProperty(fmp, **parent);
      } 
    }
    return;
  }
  Object::fillProperty(fmp, p);
}

void Worldline::setParameters(FactoryMessenger* fmp) {
  wait_pos_ = 1;
  metric(fmp->metric());
  Object::setParameters(fmp);
  wait_pos_ = 0;
  if (init_vel_) {
    delete[] init_vel_; init_vel_=NULL;
    throwError("Worldline::setParameters(): "
	       "Velocity was found but not Position");
  }
}
#endif

int Worldline::setParameter(std::string name,
			    std::string content,
			    std::string unit) {
  double coord[8];
  char* tc = const_cast<char*>(content.c_str());
  if (name=="InitialCoordinate") {
    name=="InitCoord";
    return Object::setParameter(name, content, unit);
  } else if (name=="Position") {
    if (FactoryMessenger::parseArray(content, coord, 4) != 4)
      throwError("Worldline \"Position\" requires exactly 4 tokens");
    if (init_vel_) {
      setInitCoord(coord, init_vel_);
      delete[] init_vel_; init_vel_=NULL;
    } else setPosition(coord);
    wait_pos_ = 0;
  } else if (name=="Velocity") {
    if (FactoryMessenger::parseArray(content, coord, 3) != 3)
      throwError("Worldline \"Velocity\" requires exactly 3 tokens");
    if (wait_pos_) {
      if (init_vel_) delete [] init_vel_;
      init_vel_ = new double[3];
      memcpy(init_vel_, coord, 3*sizeof(double));
    } else setVelocity(coord);
  }
  else return Object::setParameter(name, content, unit);
  return 0;
}

void Worldline::integrator(std::string const &type) {
  if (type=="Legacy") state_ = new IntegState::Legacy(this);
#ifdef HAVE_BOOST
  else state_ = new IntegState::Boost(this, type);
#else
  else throwError("unrecognized integrator (recompile with boost?)");
#endif
}

std::string Worldline::integrator() const {
  return state_->kind();
}

SmartPointer<Metric::Generic> Worldline::metric() const { return metric_; }

string Worldline::className() const { return  string("Worldline"); }
string Worldline::className_l() const { return  string("worldline"); }

void Worldline::initCoord(std::vector<double> const &v) {
  if (v.size() != 8)
    throwError("Worldline::initCoord() requires an 8-element vector");
  double c[8];
  for (size_t i=0; i<8; ++i) c[i]=v[i];
  setInitCoord(c);
}

std::vector<double> Worldline::initCoord() const {
  std::vector<double> coord(8, 0.);
  coord[0] = x0_[i0_];
  coord[1] = x1_[i0_];
  coord[2] = x2_[i0_];
  coord[3] = x3_[i0_];
  coord[4] = x0dot_[i0_];
  coord[5] = x1dot_[i0_];
  coord[6] = x2dot_[i0_];
  coord[7] = x3dot_[i0_];
  return coord;
}

void Worldline::setInitCoord(const double coord[8], int dir) {
  GYOTO_DEBUG_ARRAY(coord, 8);
  if (dir==0) dir = getMass() ? 1 : -1;
  imin_=imax_=i0_=(dir==1?0:x_size_-1);
  x0_[i0_]=coord[0];
  x1_[i0_]=coord[1];
  x2_[i0_]=coord[2];
  x3_[i0_]=coord[3];
  x0dot_[i0_]=coord[4];
  x1dot_[i0_]=coord[5];
  x2dot_[i0_]=coord[6];
  x3dot_[i0_]=coord[7];
  reInit();
}

void Worldline::setInitCoord(double pos[4], double v[3], int dir) {
  if (!getMass())
    throwError("Worldline::setInitCoord(pos, vel) "
	       "only makes sense for massive particles");
  if (!metric_)
    throwError("Please set metric before calling "
	       "Worldline::setInitCoord(double pos[4], double vel[3])");
  double tdot0=metric_->SysPrimeToTdot(pos, v);
  GYOTO_DEBUG_EXPR(tdot0);
  double coord[8]={pos[0], pos[1], pos[2], pos[3],
		   tdot0, v[0]*tdot0, v[1]*tdot0, v[2]*tdot0};
  GYOTO_DEBUG_ARRAY(coord, 8);
  setInitCoord(coord, dir);
}

void Worldline::setInitialCondition(SmartPointer<Metric::Generic> met,
				    const double coord[8],
				    const int dir)
{
  metric(met);
  setInitCoord(coord, dir);
}

void Worldline::setPosition(double pos[4]) {
  double vel[] = {0., 0., 0.};
  setInitCoord(pos, vel);
}

void Worldline::setVelocity(double vel[3]) {
  double coord[8];
  getInitialCoord(coord);
  setInitCoord(coord, vel);
}



void Worldline::reset() { if (imin_<=imax_) imin_=imax_=i0_; }
void Worldline::reInit() {
  if (imin_ <= imax_) {
    reset();
    double coord[8];
    getInitialCoord(coord);
    GYOTO_DEBUG_ARRAY(coord, 8);
    if (metric_) {
      if ( metric_() -> coordKind() == GYOTO_COORDKIND_SPHERICAL 
	   && x2_[i0_]==0. ) {
	if (verbose() >= GYOTO_SEVERE_VERBOSITY)
	  cerr << "SEVERE: Worldline::reInit(): Kicking particle off z axis\n";
	x2_[i0_]=coord[2]=1e-10;
      }
      metric_ -> setParticleProperties(this,coord);
    }
  }
}


void Worldline::xStore(size_t ind, double coord[8])
{
  x0_[ind] = coord[0];
  x1_[ind] = coord[1];
  x2_[ind] = coord[2];
  x3_[ind] = coord[3];
  x0dot_[ind] = coord[4];
  x1dot_[ind] = coord[5];
  x2dot_[ind] = coord[6];
  x3dot_[ind] = coord[7];
}

void Worldline::xFill(double tlim) {

  //GYOTO_DEBUG<< "In xFill" << endl;
  int dir;
  stopcond=0;
  size_t ind;

  // Check whether anything needs to be done,
  // Determine direction,
  // Allocate memory.
  if (tlim > x0_[imax_]) {
    // integrate forward
    dir = 1;
    ind = (imax_==x_size_-1)?xExpand(1):imax_;
  } else if (tlim < x0_[imin_]) {
    // integrate backward
    dir = -1;
    ind = (imin_==0)?xExpand(-1):imin_;
  } else return ; // nothing to do

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG<< "Integrating worldline " ;
# endif
  
  // Set up integration
  double MassPart=getMass();
  if (MassPart==1.) {
#   if GYOTO_DEBUG_ENABLED
    GYOTO_IF_DEBUG
    cerr << "of massive particule ....." << endl;
    GYOTO_ENDIF_DEBUG
#   endif
  }else if(MassPart==0.){
#   if GYOTO_DEBUG_ENABLED
    GYOTO_IF_DEBUG
    cerr << "of 0-mass particule ....." << endl;
    GYOTO_ENDIF_DEBUG
#   endif
  }else{
    throwError("In Worldline.C Unrecognized mass.");
    //GYOTO_DEBUG<< "of unrecognized mass (!!) particule ....." << endl;
    //equations of geodesics written for a mass=1 star
  }
 
  double coord[8]={x0_[ind], x1_[ind], x2_[ind], x3_[ind],
		   x0dot_[ind], x1dot_[ind], x2dot_[ind], x3dot_[ind]};

  GYOTO_DEBUG << "IntegState initialization" << endl;

  state_->init(this, coord, dir*delta_);
    //delta_ = initial integration step (defaults to 0.01)

  GYOTO_DEBUG << "IntegState initialized" << endl;

  size_t mycount=0;// to prevent infinite integration

  while (!stopcond) {
    mycount++;
    ind+=dir;

    stopcond= state_ -> nextStep(coord);

    //if (stopcond && debug()) cout << "stopcond from integrator" << endl;
    if (mycount==maxiter_) {
      stopcond=1;
      Error ( "***WARNING STOP: in Worldline.C unexplained stop !!!" );
    }
    // store particle's trajectory for later use
    xStore(ind, coord);

    // Check stop condition and whether we need to expand the arrays
    if (dir==1) {
      if (coord[0]>tlim) stopcond=1;
      if ((!stopcond) & (ind==x_size_-1)) {
	imax_=x_size_-1;
	ind=xExpand(1);
      }
    } else {
      if (coord[0]<tlim) {
	stopcond=1;
      }
      if ((!stopcond) & (ind==0)) {
	imin_=0;
	ind=xExpand(-1);
      }
    }
  }
  if (dir==1) imax_=ind; // tell when we actually stopped integrating
  else imin_=ind;

}

size_t Worldline::get_nelements() const { return imax_-imin_+1; }

size_t Worldline::getImin() const {return imin_;}
size_t Worldline::getImax() const {return imax_;}
size_t Worldline::getI0() const {return i0_;}

void Worldline::get_t(double *dest) const
{ memcpy(dest, x0_+imin_, sizeof(double)*(imax_-imin_+1)); }

void Worldline::get_xyz(double *x, double *y, double *z) const {
  size_t n;
  int coordkind = metric_ -> coordKind();
  switch(coordkind) {
  case GYOTO_COORDKIND_SPHERICAL: 
    for (n=imin_;n<=imax_;++n) {
      x[n-imin_]=x1_[n]*sin(x2_[n])*cos(x3_[n]);
      y[n-imin_]=x1_[n]*sin(x2_[n])*sin(x3_[n]);
      z[n-imin_]=x1_[n]*cos(x2_[n]);
    }
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    for (n=imin_;n<=imax_;++n) {
      x[n-imin_]=x1_[n];
      y[n-imin_]=x2_[n];
      z[n-imin_]=x3_[n];
    }
    break;
  default: Gyoto::throwError("in Worldline::get_xyz: Incompatible coordinate kind");
  }

}
void Worldline::getCartesian(double const * const dates, size_t const n_dates,
			     double * const x, double * const y,
			     double * const z, double * const xprime,
			     double * const yprime,  double * const zprime)
{
  double *x1, *x2, *x3, *x0dot, *x1dot, *x2dot, *x3dot, tauprime;
  int coordkind = metric_ -> coordKind();
  size_t di;
  double rprime, thetaprime, phiprime, costheta, sintheta, cosphi, sinphi, r;

  x0dot = new double[n_dates];
  x1dot = new double[n_dates];
  x2dot = new double[n_dates];
  x3dot = new double[n_dates];

  switch(coordkind) {
  case GYOTO_COORDKIND_SPHERICAL: 
    x1    = new double[n_dates];
    x2    = new double[n_dates];
    x3    = new double[n_dates];
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    x1=x; x2=y; x3=z;
    break;
  default: 
    Gyoto::throwError("in Worldline::get_xyz: unknown coordinate kind");
    x1=x2=x3=NULL; // fix warning
  }

  getCoord(dates, n_dates, x1, x2, x3, x0dot, x1dot, x2dot, x3dot);
  
  switch(coordkind) {
  case GYOTO_COORDKIND_SPHERICAL: 
    for (di=0; di<n_dates; ++di) {
      r = x1[di];
      costheta = cos(x2[di]); sintheta = sin(x2[di]); 
      cosphi   = cos(x3[di]); sinphi   = sin(x3[di]);
      // x, y, z
      x[di] = r * sintheta * cosphi;
      y[di] = r * sintheta * sinphi;
      z[di] = r * costheta;

      if (xprime || yprime || zprime) {
	// dx/dt, dy/dt, dz/dt
	tauprime   = 1./x0dot[di];
	rprime     = x1dot[di]*tauprime;
	thetaprime = x2dot[di]*tauprime;
	phiprime   = x3dot[di]*tauprime;
	if (xprime)
	  xprime[di] = rprime * sintheta * cosphi
	    + r * thetaprime * costheta * cosphi
	    - r * phiprime * sintheta * sinphi;
	if (yprime)
	  yprime[di] = rprime * sintheta * sinphi
	    + r * thetaprime * costheta * sinphi
	    + r * phiprime * cosphi;
	if (zprime)
	  zprime[di] = rprime * costheta
	    - r * thetaprime * sintheta
	    ;
      }
    }
    delete [] x1;
    delete [] x2;
    delete [] x3;
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    if (xprime || yprime || zprime) {
      for (di=0; di<n_dates; ++di) {
	tauprime = 1./x0dot[di];
	if (xprime) xprime[di] = x1dot[di]*tauprime;
	if (yprime) yprime[di] = x2dot[di]*tauprime;
	if (zprime) zprime[di] = x3dot[di]*tauprime;
      }
    }
    break;
  default: Gyoto::throwError("in Worldline::get_xyz: unknown coordinate kind");
  }

  delete [] x0dot;
  delete [] x1dot;
  delete [] x2dot;
  delete [] x3dot;
}

void Worldline::getCoord(double const * const dates, size_t const n_dates,
			 double * const x1,
			 double * const x2,    double * const x3,
			 double * const x0dot, double * const x1dot,
			 double * const x2dot, double * const x3dot)
{

  size_t curl=imin_, curm, curh=imax_;
       // Current indices for the binary search. 
       // To avoid overflows, compute curm as curl+(curh-curl)/2
  size_t di=0; // current date index
  double date; // current date

  // For the interpolation
  double bestl[8], besth[8], resl[8], resh[8]; // i/o for myrk4
  double factl, facth;
  double tausecond, dtaul, dtauh, dtl, dth, Dt, Dtm1, tauprimel, tauprimeh;
  double second, primel, primeh, pos[4], vel[3], tdot;
  int i;
  stringstream ss;


  for (di=0; di<n_dates; ++di) {
    date = dates[di];
    if (date == x0_[imax_]) {
      if (x1)       x1[di] =    x1_[imax_];
      if (x2)       x2[di] =    x2_[imax_];
      if (x3)       x3[di] =    x3_[imax_];
      if (x0dot) x0dot[di] = x0dot_[imax_];
      if (x1dot) x1dot[di] = x1dot_[imax_];
      if (x2dot) x2dot[di] = x2dot_[imax_];
      if (x3dot) x3dot[di] = x3dot_[imax_];
      if (metric_->coordKind() == GYOTO_COORDKIND_SPHERICAL){
	double pos2[8]={0.,0.,x2[di],x3[di],0.,0.,x2dot[di],0.};
	checkPhiTheta(pos2);
	x2[di]=pos2[2];x3[di]=pos2[3];x2dot[di]=pos2[6];
      }
      continue;
    } else if (date > x0_[imax_]) {
      curl=imax_;    // current imax_
      xFill(date);   // integrate, that changes imax_
      curh=imax_;    // new imax_
      if (curl == curh || date > x0_[imax_]) {
	ss<<"Worldline::getCoord: can't get coordinates for date="<<date;
	throwError(ss.str());
      }
    } else if (date < x0_[imin_]) {
      curh=x_size_-imin_; // trick if line is expanded during xFill()
      xFill(date);   // integrate, that changes imin_
      curh=x_size_-curh;
      curl=imin_;    // new imin_
      if (curl == curh || date < x0_[imin_]) {
	ss<<"Worldline::getCoord: can't get coordinates for date="<<date;
	throwError(ss.str());
      }
    } else if (date >= x0_[curh]) {
      curl=curh;
      curh=imax_;
    } else if (date < x0_[curl]) {
      curh=curl;
      curl=imin_;
    }

    while (curh-curl>1) {
      curm = curl+(curh-curl)/2;
      if (date >= x0_[curm]) curl = curm;
      else curh = curm;
    }

    if (date == x0_[curl]) {
      if (x1)       x1[di] =    x1_[curl];
      if (x2)       x2[di] =    x2_[curl];
      if (x3)       x3[di] =    x3_[curl];
      if (x0dot) x0dot[di] = x0dot_[curl];
      if (x1dot) x1dot[di] = x1dot_[curl];
      if (x2dot) x2dot[di] = x2dot_[curl];
      if (x3dot) x3dot[di] = x3dot_[curl];      
      if (metric_->coordKind() == GYOTO_COORDKIND_SPHERICAL){
	double pos2[8]={0.,0.,x2[di],x3[di],0.,0.,x2dot[di],0.};
	checkPhiTheta(pos2);
	x2[di]=pos2[2];x3[di]=pos2[3];x2dot[di]=pos2[6];
      }
      continue;
    }

    // Attempt to get closer to the specified date using the
    // integrator.
    dtl=date-x0_[curl]; dth=date-x0_[curh];
    Dt=(x0_[curh]-x0_[curl]); Dtm1=1./Dt;
    tauprimel=1./x0dot_[curl]; tauprimeh=1./x0dot_[curh];
    tausecond = (tauprimeh-tauprimel)*Dtm1;
    // tauprime(dt)=tauprime(t0)+tausecond*dt
    dtaul=tauprimel*dtl+0.5*tausecond*dtl*dtl;
    dtauh=tauprimeh*dth+0.5*tausecond*dth*dth;

#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG
	   << "curl=" << curl << ", x0_[curl]=" << x0_[curl]
	   << ", curh=" << curh << ", x0_[curh]=" << x0_[curh]
	   <<endl;
#   endif

    // from below...
    bestl[0] =    x0_[curl];
    bestl[1] =    x1_[curl];
    bestl[2] =    x2_[curl];
    bestl[3] =    x3_[curl];
    bestl[4] = x0dot_[curl];
    bestl[5] = x1dot_[curl];
    bestl[6] = x2dot_[curl];
    bestl[7] = x3dot_[curl];
    state_ -> doStep(bestl, dtaul, resl);

    // from above...
    besth[0] =    x0_[curh];
    besth[1] =    x1_[curh];
    besth[2] =    x2_[curh];
    besth[3] =    x3_[curh];
    besth[4] = x0dot_[curh];
    besth[5] = x1dot_[curh];
    besth[6] = x2dot_[curh];
    besth[7] = x3dot_[curh];
    state_ -> doStep(besth, dtauh, resh);

#   if GYOTO_DEBUG_ENABLED
    GYOTO_IF_DEBUG
      GYOTO_DEBUG_ARRAY(bestl, 8);
      GYOTO_DEBUG_EXPR(dtaul);
      GYOTO_DEBUG_ARRAY(resl, 8);
      GYOTO_DEBUG <<   "tl="         << resl[0] 
		  << ", date="       << date 
		  << ", th="         << resh[0]
		  << ", th-tl="      << resh[0]-resl[0] 
		  << ", Dt="         << Dt 
		  <<   "(th-tl)/Dt=" << (resh[0]-resl[0])*Dtm1
		  << endl;
    GYOTO_ENDIF_DEBUG
#   endif

    // Now sometimes we actually got further away. We have 4 dates
    // well estimated, take the 2 best, 1 above, 1 below.
    if (resl[0]<=date) {
      if (resl[0] > bestl[0]) memcpy(bestl, resl, 8*sizeof(double));
    } else {
      if (resl[0] < besth[0]) memcpy(besth, resl, 8*sizeof(double));
    }

    if (resh[0]<=date) {
      if (resh[0] > bestl[0]) memcpy(bestl, resh, 8*sizeof(double));
    } else {
      if (resh[0] < besth[0]) memcpy(besth, resh, 8*sizeof(double));
    }

#   if GYOTO_DEBUG_ENABLED
    GYOTO_DEBUG
	   << "x0_[curl]=" << x0_[curl]
	   << ", bestl[0]=" << bestl[0]
	   << ", date=" << date
	   << ", besth[0]=" << besth[0]
	   << ", x0_[curh]=" << x0_[curh]
	   << ", besth[0]-bestl[0]=" << besth[0]-bestl[0]
	   << ", Dt=" << Dt
	   << ", (th-tl)/Dt=" << (besth[0]-bestl[0])*Dtm1 << endl;
#   endif

    // Now interpolate as best we can
    if (bestl[0]==date) {
      if (x1)       x1[di] = bestl[1];
      if (x2)       x2[di] = bestl[2];
      if (x3)       x3[di] = bestl[3];
      if (x0dot) x0dot[di] = bestl[4];
      if (x1dot) x1dot[di] = bestl[5];
      if (x2dot) x2dot[di] = bestl[6];
      if (x3dot) x3dot[di] = bestl[7];
    }
    if (besth[0]==date) {
      if (x1)       x1[di] = besth[1];
      if (x2)       x2[di] = besth[2];
      if (x3)       x3[di] = besth[3];
      if (x0dot) x0dot[di] = besth[4];
      if (x1dot) x1dot[di] = besth[5];
      if (x2dot) x2dot[di] = besth[6];
      if (x3dot) x3dot[di] = besth[7];
    }

    dtl=date-bestl[0]; Dt=besth[0]-bestl[0]; Dtm1=1./Dt;
    facth=dtl*Dtm1; factl=1.-facth;
    if (getMass()) {
#     if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "massive particle, interpolating\n";
#     endif
      // Star (massive particle)
      tauprimel=1./bestl[4]; tauprimeh=1./besth[4];

      pos[0] = date;
      
      for (i=1; i<=3; ++i) {
	primel=bestl[i+4]*tauprimel;
	primeh=besth[i+4]*tauprimeh;
	vel[i-1]=primel*factl+primeh*facth;
	second =(primeh-primel)*Dtm1;
	//	pos[i] = resl[i] + primel*dtl + 0.5*second*dtl*dtl;
	pos[i] = bestl[i] + primel*dtl + 0.5*second*dtl*dtl;
      }

      tdot=metric_->SysPrimeToTdot(pos, vel);

      if (x1)       x1[di] = pos[1];
      if (x2)       x2[di] = pos[2];
      if (x3)       x3[di] = pos[3];
      if (x0dot) x0dot[di] = tdot;
      if (x1dot) x1dot[di] = vel[0]*tdot;
      if (x2dot) x2dot[di] = vel[1]*tdot;
      if (x3dot) x3dot[di] = vel[2]*tdot;
    } else {
      // Photon: don't be so elaborate, we certainly don't need it... yet
      if (x1)       x1[di] = bestl[1]*factl + besth[1]*facth;
      if (x2)       x2[di] = bestl[2]*factl + besth[2]*facth;
      if (x3)       x3[di] = bestl[3]*factl + besth[3]*facth;
      if (x0dot) x0dot[di] = bestl[4]*factl + besth[4]*facth;
      if (x1dot) x1dot[di] = bestl[5]*factl + besth[5]*facth;
      if (x2dot) x2dot[di] = bestl[6]*factl + besth[6]*facth;
      if (x3dot) x3dot[di] = bestl[7]*factl + besth[7]*facth;
    }

    /* For spherical-like coordinates,
       transforms theta and phi so that 
       theta is in [0,pi] and phi in [0,2pi] 
       This call is due to interpolation above
       that could lead theta,phi out of their
       correct ranges.
    */
    if (metric_->coordKind() == GYOTO_COORDKIND_SPHERICAL
	&& x2 && x3 && x2dot){
      double pos2[8]={0.,0.,x2[di],x3[di],0.,0.,x2dot[di],0.};
      checkPhiTheta(pos2);
      x2[di]=pos2[2];x3[di]=pos2[3];x2dot[di]=pos2[6];
    }
    
  }

}
void Worldline::getCoord(double *x0dest,
			  double *x1dest, double *x2dest, double *x3dest)
			 const {
  //if (sysco!=sys_)
  //Gyoto::throwError("At this point, coordinate conversion is not implemented");
  size_t ncomp=imax_-imin_+1;
  memcpy(x0dest, x0_+imin_, sizeof(double)*ncomp);
  memcpy(x1dest, x1_+imin_, sizeof(double)*ncomp);
  memcpy(x2dest, x2_+imin_, sizeof(double)*ncomp);
  memcpy(x3dest, x3_+imin_, sizeof(double)*ncomp);
}

void Worldline::checkPhiTheta(double coord[8]) const{
  switch (metric_ -> coordKind()) {
  case GYOTO_COORDKIND_SPHERICAL:
    {
    /* Transforms theta and phi in coord so that 
       theta is in [0,pi] and phi in [0,2pi] */
    double thetatmp=coord[2], phitmp=coord[3];
    while (thetatmp>M_PI) thetatmp-=2.*M_PI;
    while (thetatmp<-M_PI) thetatmp+=2.*M_PI;//then theta in [-pi,pi]
    if (thetatmp<0.) {
      thetatmp*=-1.;//then theta in [0,pi]
      coord[6]*=-1.; //theta -> -theta then idem for derivative
      phitmp+=M_PI;//thus, same point x,y,z
    }
    while (phitmp>2.*M_PI) phitmp-=2.*M_PI;
    while (phitmp<0.) phitmp+=2.*M_PI;//then phi in [0,2pi]
    coord[2]=thetatmp;
    coord[3]=phitmp;
    }
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    throwError("Worldline::checkPhiTheta(): should not be called "
	       "with cartesian-like coordinates");
  default:
    throwError("Worldline::checkPhiTheta(): unknown COORDKIND");
  }
}

void Worldline::get_dot(double *x0dest, double *x1dest, double *x2dest, double *x3dest) const {
  //  if (sysco!=sys_)
  //  Gyoto::throwError("At this point, coordinate conversion is not implemented");
  size_t ncomp=imax_-imin_+1;
  memcpy(x0dest, x0dot_+imin_, sizeof(double)*ncomp);
  memcpy(x1dest, x1dot_+imin_, sizeof(double)*ncomp);
  memcpy(x2dest, x2dot_+imin_, sizeof(double)*ncomp);
  memcpy(x3dest, x3dot_+imin_, sizeof(double)*ncomp);
}

void Worldline::getSkyPos(SmartPointer<Screen> screen, double *dalpha, double *ddelta, double *dD) const {
  double pos[4], skypos[3];
  size_t i, ncomp=imax_-imin_+1;
  for (i=0;i<ncomp;++i) {
    pos[0]=x0_[i+imin_];
    pos[1]=x1_[i+imin_];
    pos[2]=x2_[i+imin_];
    pos[3]=x3_[i+imin_];
    screen -> coordToSky(pos, skypos);
    dalpha[i]=skypos[0];
    ddelta[i]=skypos[1];
    dD[i]    =skypos[2];
  }
}

void Worldline::get_prime(double *x1dest, double *x2dest, double *x3dest) const {
  size_t n;
  double f;
  //if (sysco!=sys_)
  //  Gyoto::throwError("At this point, coordinate conversion is not implemented");
  for (n=imin_;n<=imax_;++n) {
    x1dest[n-imin_]=x1dot_[n]*(f=1./x0dot_[n]);
    x2dest[n-imin_]=x2dot_[n]*f;
    x3dest[n-imin_]=x3dot_[n]*f;
  }
}

void Worldline::save_txyz(char * filename) const {
  size_t n;
  int width=15;
  int prec=12;
  ofstream fichierxyz(filename);
  int coordkind = metric_ -> coordKind();
  switch(coordkind) {
  case GYOTO_COORDKIND_SPHERICAL:
    for (n=imin_;n<=imax_;++n) {
      //      GYOTO_DEBUG<< "dans save imin, coord= " << imin_ << " " << x0_[n] << " " << x1_[n] << " "  << x2_[n] << " " << x3_[n] << endl;
      //fichierxyz << setprecision(prec) << setw(width) << x1_[n] << "  ";
      //fichierxyz << setprecision(prec) << setw(width) << x0[n] << "  ";
      fichierxyz << setprecision(prec) << setw(width) << x0_[n] << "  ";//saving r distance
      fichierxyz << setprecision(prec) << setw(width) << x1_[n]*sin(x2_[n])*cos(x3_[n]) << "  ";
      fichierxyz << setprecision(prec) << setw(width) << x1_[n]*sin(x2_[n])*sin(x3_[n]) << "  ";
      fichierxyz << setprecision(prec) << setw(width) << x1_[n]*cos(x2_[n]) << endl;
    }
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    for (n=imin_;n<=imax_;++n) {
      //fichierxyz << setprecision(prec) << setw(width) << sqrt(x1[n]*x1[n]+x2_[n]*x2_[n]+x3_[n]*x3_[n]) << "  ";
      fichierxyz << setprecision(prec) << setw(width) << x0_[n] << "  ";
      fichierxyz << setprecision(prec) << setw(width) << x1_[n] << "  ";
      fichierxyz << setprecision(prec) << setw(width) << x2_[n] << "  ";
      fichierxyz << setprecision(prec) << setw(width) << x3_[n] << endl;
    }
    break;
  default: Gyoto::throwError("in Worldline::save_xyz: Incompatible coordinate kind");
  }
  
  fichierxyz.close();
}

void Worldline::save_txyz(char * filename, const double t1, const double mass_sun, const double distance_kpc, const string unit, SmartPointer<Screen> sc) {
  xFill(t1);
  size_t nelem = get_nelements(), n=0;
  double * t = new double [nelem];
  double * x = new double [nelem];
  double * y = new double [nelem];
  double * z = new double [nelem];
  double rad2deg = 180./M_PI;
  double f=1.;
  if (getMass()) f /= x0dot_[i0_];
  int width=GYOTO_WIDTH;
  int prec=GYOTO_PREC;
  ofstream fichierxyz(filename);

  get_t(t);
  get_xyz(x, y, z);

  convert(x, nelem, mass_sun, distance_kpc, unit);
  convert(y, nelem, mass_sun, distance_kpc, unit);
  convert(z, nelem, mass_sun, distance_kpc, unit);

  // HEADER

  string metkind = metric_->kind();

  fichierxyz <<   "# Gyoto save file " << endl
	     <<   "# Start Gyoto parameters" << endl
	     <<   "# particle_type = \"" << className_l() << "\"" << endl;
//   if (!metkind.compare("KerrBL")) {
//     SmartPointer<KerrBL> kmet = metric_ ;
//     fichierxyz << "#   metric_type = \"kerr\"" << endl
// 	       << "#          spin = " << kmet -> getSpin() << endl
// 	       << "#          mass = " << mass_sun << endl;
//   } else
    fichierxyz << "#   metric_type = \"" << metkind << "\"" << endl ;
  fichierxyz <<   "#            t0 = " << x0_[i0_] << endl
	     <<   "#            r0 = " << x1_[i0_] << endl
	     <<   "#        theta0 = " << x2_[i0_] << endl
	     <<   "#          phi0 = " << x3_[i0_] << endl
	     <<   "#       rprime0 = " << x1dot_[i0_]*f << endl
	     <<   "#   thetaprime0 = " << x2dot_[i0_]*f << endl
	     <<   "#     phiprime0 = " << x3dot_[i0_]*f << endl
	     <<   "#            t1 = " << t1 << endl;

  if (sc)
    fichierxyz << "#          incl = " << sc -> inclination()*rad2deg << endl
	       << "#          paln = " << sc -> PALN()*rad2deg << endl
	       << "#         phase = " << sc -> argument()*rad2deg << endl
	       << "#      distance = " << distance_kpc << endl;

  fichierxyz <<   "#   length_unit = \"" << unit << "\"" << endl
	     <<   "# End Gyoto parameters" << endl
	     <<   "# Columns are t, x, y, z" << endl;


  // DATA

  for (n=0;n<=nelem;++n) {
    fichierxyz << setprecision(prec) << setw(width) << t[n] << "  "
	       << setprecision(prec) << setw(width) << x[n] << "  "
	       << setprecision(prec) << setw(width) << y[n] << "  "
	       << setprecision(prec) << setw(width) << z[n] << endl ;
  }

  delete [] t; 
  delete [] x; 
  delete [] y; 
  delete [] z; 
}


void Worldline::delta(const double del) { delta_=del; }
void Worldline::delta(double d, const string &unit) {
  delta(Units::ToGeometrical(d, unit, metric_));
}

double Worldline::delta() const { return delta_; }
double Worldline::delta(const string &unit) const {
  return Units::FromGeometrical(delta(), unit, metric_);
}

double Worldline::tMin() const { return tmin_; }
double Worldline::tMin(const string &unit) const {
  return Units::FromGeometricalTime(tMin(), unit, metric_);
}

void Worldline::tMin(double tmin) { tmin_ = tmin; }
void Worldline::tMin(double tmin, const string &unit) {
  tMin(Units::ToGeometricalTime(tmin, unit, metric_));
}

void Worldline::adaptive(bool mode) { adaptive_ = mode; state_->init();}
bool Worldline::adaptive() const { return adaptive_; }

void Worldline::secondary(bool sec) { secondary_ = sec; }
bool Worldline::secondary() const { return secondary_; }

void Worldline::maxiter(size_t miter) { maxiter_ = miter; }
size_t Worldline::maxiter() const { return maxiter_; }

double Worldline::deltaMin() const {return delta_min_;}
double Worldline::deltaMax() const {return delta_max_;}
void  Worldline::deltaMin(double h1) {
  delta_min_=h1;
  state_->init();
}
void  Worldline::deltaMax(double h1) {delta_max_=h1;}
double Worldline::deltaMaxOverR() const { return delta_max_over_r_;}
void Worldline::deltaMaxOverR(double t) {delta_max_over_r_=t;}

double Worldline::absTol() const {return abstol_;}
void Worldline::absTol(double t) {abstol_=t; state_->init();}
double Worldline::relTol() const {return reltol_;}
void Worldline::relTol(double t) {reltol_=t; state_->init();}

double Worldline::deltaMax(double const pos[8], double h1max) const
{
  double h1max_at_r=abs(pos[1]);
  if (metric_ -> coordKind()==GYOTO_COORDKIND_CARTESIAN) {
    double tmp;
    if ((tmp=abs(pos[2]))>h1max_at_r) h1max_at_r=tmp;
	if ((tmp=abs(pos[3]))>h1max_at_r) h1max_at_r=tmp;
  }
  h1max_at_r *= delta_max_over_r_;
  if (h1max > h1max_at_r) h1max = h1max_at_r;
  if (h1max>delta_max_) h1max=delta_max_;
  if (h1max<delta_min_) h1max=delta_min_;
  return h1max;
}

double const * Worldline::getCst() const {
  return cst_;
}

void Worldline::setCst(double const * const cst, const size_t n) {
  if (cst_) delete [] cst_;
  cst_ = new double[n];
  cst_n_ = n;
  for (size_t ii=0;ii<n;ii++) cst_[ii]=cst[ii];
}

void Worldline::getInitialCoord(double coord[8]) const {
  if (imax_<imin_)
    throwError("Worldline::getInitialCoord(): initial coordinate not set yet");
  coord[0] = x0_[i0_];
  coord[1] = x1_[i0_];
  coord[2] = x2_[i0_];
  coord[3] = x3_[i0_];
  coord[4] = x0dot_[i0_];
  coord[5] = x1dot_[i0_];
  coord[6] = x2dot_[i0_];
  coord[7] = x3dot_[i0_];
}
void Worldline::getCoord(size_t index, double coord[8]) const {
  //GYOTO_DEBUG<< "index=" << index << endl;
  //GYOTO_DEBUG<< "x0[index]= " << x1dot_[index] << endl;
  //GYOTO_DEBUG<< "index == " << index << endl;
  if (index<imin_ || index>imax_) {
    cerr << "Indices min curr max= " << imin_ << " " << index << " " << imax_ << endl;
    throwError("Worldline::getCoord: bad index");
  }
  coord[0] = x0_[index];
  coord[1] = x1_[index];
  coord[2] = x2_[index];
  coord[3] = x3_[index];
  coord[4] = x0dot_[index];
  coord[5] = x1dot_[index];
  coord[6] = x2dot_[index];
  coord[7] = x3dot_[index];
}

void Worldline::getCartesianPos(size_t index, double dest[4]) const {
  dest[0] = x0_[index];
  int coordkind = metric_ -> coordKind();
  switch (coordkind) {
  case GYOTO_COORDKIND_SPHERICAL:
    dest[1] = x1_[index]*sin(x2_[index])*cos(x3_[index]);
    dest[2] = x1_[index]*sin(x2_[index])*sin(x3_[index]);
    dest[3] = x1_[index]*cos(x2_[index]);
    break;
  case GYOTO_COORDKIND_CARTESIAN:
    dest[1] = x1_[index];
    dest[2] = x2_[index];
    dest[3] = x3_[index];
    break;
  default: Gyoto::throwError("Worldline::getCartesianPos: Incompatible coordinate kind");
  }
}
