/*
    Copyright 2011-2015, 2017-2020 Frederic Vincent, Thibaut Paumard

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

#ifdef GYOTO_HAVE_BOOST_INTEGRATORS
# define _GYOTO_DEFAULT_INTEGRATOR "runge_kutta_fehlberg78"
#else
# define _GYOTO_DEFAULT_INTEGRATOR "Legacy"
#endif


Worldline::Worldline() : ep0_(NULL), ep1_(NULL), ep2_(NULL), ep3_(NULL),
			 et0_(NULL), et1_(NULL), et2_(NULL), et3_(NULL),
                         stopcond(0), metric_(NULL),
                         imin_(1), i0_(0), imax_(0), adaptive_(1),
			 secondary_(1), 
			 parallel_transport_(false),
			 delta_(GYOTO_DEFAULT_DELTA),
			 tmin_(-DBL_MAX), cst_(NULL), cst_n_(0),
			 wait_pos_(0), init_vel_(NULL),
			 maxiter_(GYOTO_DEFAULT_MAXITER),
			 delta_min_(GYOTO_DEFAULT_DELTA_MIN),
			 delta_max_(GYOTO_DEFAULT_DELTA_MAX),
			 delta_max_over_r_(GYOTO_DEFAULT_DELTA_MAX_OVER_R),
			 abstol_(GYOTO_DEFAULT_ABSTOL),
			 reltol_(GYOTO_DEFAULT_RELTOL),
			 maxCrossEqplane_(DBL_MAX),
			 state_(NULL)
{ 
  xAllocate();
  integrator(_GYOTO_DEFAULT_INTEGRATOR);
}

Worldline::Worldline(const Worldline& orig) :
  stopcond(orig.stopcond),
  ep0_(NULL), ep1_(NULL), ep2_(NULL), ep3_(NULL),
  et0_(NULL), et1_(NULL), et2_(NULL), et3_(NULL),
  metric_(NULL),
  x_size_(orig.x_size_), imin_(orig.imin_), i0_(orig.i0_), imax_(orig.imax_),
  adaptive_(orig.adaptive_), secondary_(orig.secondary_),
  parallel_transport_(orig.parallel_transport_),
  delta_(orig.delta_), tmin_(orig.tmin_), cst_(NULL), cst_n_(orig.cst_n_),
  wait_pos_(orig.wait_pos_), init_vel_(NULL),
  maxiter_(orig.maxiter_),
  delta_min_(orig.delta_min_),
  delta_max_(orig.delta_max_),
  delta_max_over_r_(orig.delta_max_over_r_),
  abstol_(orig.abstol_),
  reltol_(orig.reltol_),
  maxCrossEqplane_(orig.maxCrossEqplane_),
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
  memcpy(tau_+imin_, orig.tau_+imin_, sz);
  memcpy(x0_+imin_, orig.x0_+imin_, sz);
  memcpy(x1_+imin_, orig.x1_+imin_, sz);
  memcpy(x2_+imin_, orig.x2_+imin_, sz);
  memcpy(x3_+imin_, orig.x3_+imin_, sz);
  memcpy(x0dot_+imin_, orig.x0dot_+imin_, sz);
  memcpy(x1dot_+imin_, orig.x1dot_+imin_, sz);
  memcpy(x2dot_+imin_, orig.x2dot_+imin_, sz);
  memcpy(x3dot_+imin_, orig.x3dot_+imin_, sz);
  if (parallel_transport_) {
    memcpy(ep0_+imin_, orig.ep0_+imin_, sz);
    memcpy(ep1_+imin_, orig.ep1_+imin_, sz);
    memcpy(ep2_+imin_, orig.ep2_+imin_, sz);
    memcpy(ep3_+imin_, orig.ep3_+imin_, sz);
    memcpy(et0_+imin_, orig.et0_+imin_, sz);
    memcpy(et1_+imin_, orig.et1_+imin_, sz);
    memcpy(et2_+imin_, orig.et2_+imin_, sz);
    memcpy(et3_+imin_, orig.et3_+imin_, sz);
  }
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
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done\n";
# endif
}

Worldline::Worldline(Worldline *orig, size_t i0, int dir, double step_max) :
  ep0_(NULL), ep1_(NULL), ep2_(NULL), ep3_(NULL),
  et0_(NULL), et1_(NULL), et2_(NULL), et3_(NULL),
  metric_(orig->metric_),
//  x_size_(orig.x_size_), imin_(orig.imin_), i0_(orig.i0_), imax_(orig.imax_),
  adaptive_(orig->adaptive_), secondary_(orig->secondary_),
  parallel_transport_(orig->parallel_transport_),
  delta_(orig->delta_), tmin_(orig->tmin_), cst_(NULL), cst_n_(orig->cst_n_),
  wait_pos_(orig->wait_pos_), init_vel_(NULL),
  maxiter_(orig->maxiter_),
  delta_min_(orig->delta_min_),
  delta_max_(orig->delta_max_),
  delta_max_over_r_(orig->delta_max_over_r_),
  abstol_(orig->abstol_),
  reltol_(orig->reltol_),
  maxCrossEqplane_(orig->maxCrossEqplane_),
  state_(NULL)
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

  orig->getCoord(x0_, x_size_, x1_, x2_, x3_, x0dot_, x1dot_, x2dot_, x3dot_, ep0_, ep1_, ep2_, ep3_, et0_, et1_, et2_, et3_, tau_);

# if GYOTO_DEBUG_ENABLED
  GYOTO_IF_DEBUG
    {
      GYOTO_DEBUG << "(Worldline*, "<<i0<<", "<<dir<<", "<<step_max<<")"<<endl;
      GYOTO_DEBUG << "d1="<<d1<<", d2="<<d2<<endl;
      GYOTO_DEBUG_ARRAY(tau_, x_size_);
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
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "done\n";
# endif
}

Worldline::~Worldline(){
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << endl;
# endif
  if (metric_) metric_ -> unhook(this);
  delete[] tau_;
  delete[] x0_;
  delete[] x1_;
  delete[] x2_;
  delete[] x3_;
  delete[] x0dot_;
  delete[] x1dot_;
  delete[] x2dot_;
  delete[] x3dot_;
  eDeallocate();
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
  tau_ = new double[x_size_];
  x0_ = new double[x_size_];
  x1_ = new double[x_size_];
  x2_ = new double[x_size_];
  x3_ = new double[x_size_];
  x0dot_ = new double[x_size_];
  x1dot_ = new double[x_size_];
  x2dot_ = new double[x_size_];
  x3dot_ = new double[x_size_];
  eAllocate();
}

void Worldline::xExpand(double* &x, int dir) {
  double * old;
  size_t offset=(dir==1)?0:x_size_;
  size_t i;
  size_t nsize=2*x_size_;

  old=x;
  x=new double[nsize];
  for (i=imin_;i<=imax_;++i) x[i+offset]=old[i];
  GYOTO_IF_DEBUG
  GYOTO_DEBUG_EXPR(imin_);
  GYOTO_DEBUG_EXPR(imin_+offset);
  GYOTO_DEBUG_EXPR(old[imin_]);
  GYOTO_DEBUG_EXPR(x[imin_+offset]);
  GYOTO_ENDIF_DEBUG
  delete [] old;
}

size_t Worldline::xExpand(int dir) {
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(dir);
  GYOTO_DEBUG << endl << "massive=" << getMass()
	      << ", size=" << x_size_
	      << ", imin_=" << imin_
	      << ", i0_=" << i0_
	      << ", imax_=" << imax_
	      << endl;
# endif

  xExpand(tau_, dir);
  xExpand(x0_, dir);
  xExpand(x1_, dir);
  xExpand(x2_, dir);
  xExpand(x3_, dir);
  xExpand(x0dot_, dir);
  xExpand(x1dot_, dir);
  xExpand(x2dot_, dir);
  xExpand(x3dot_, dir);
  eExpand(dir);

  size_t retval=(dir==1)?(x_size_-1):x_size_;
  size_t offset=(dir==1)?0:x_size_;

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << "retval=" << retval
	      << ", offset=" << offset
		   << ", dir=" << dir;
# endif

  x_size_*=2;

  imin_+=offset;
  i0_+=offset;
  imax_+=offset;

# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG << ", x_size_=" << x_size_
	      << ", imin_=" << imin_
	      << ", i0_=" << i0_
	      << ", imax_=" << imax_
	      << endl;
# endif

  return retval;
}

void Worldline::eAllocate()
{
# if GYOTO_DEBUG_ENABLED
  GYOTO_DEBUG_EXPR(x_size_);
# endif
  if (!x_size_ || !parallel_transport_) return;
  ep0_ = new double[x_size_];
  ep1_ = new double[x_size_];
  ep2_ = new double[x_size_];
  ep3_ = new double[x_size_];
  et0_ = new double[x_size_];
  et1_ = new double[x_size_];
  et2_ = new double[x_size_];
  et3_ = new double[x_size_];
}

void Worldline::eDeallocate() {
  if (ep0_) {delete[] ep0_; ep0_=NULL;}
  if (ep1_) {delete[] ep1_; ep1_=NULL;}
  if (ep2_) {delete[] ep2_; ep2_=NULL;}
  if (ep3_) {delete[] ep3_; ep3_=NULL;}
  if (et0_) {delete[] et0_; et0_=NULL;}
  if (et1_) {delete[] et1_; et1_=NULL;}
  if (et2_) {delete[] et2_; et2_=NULL;}
  if (et3_) {delete[] et3_; et3_=NULL;}
}

void Worldline::eExpand(int dir) {
  if (ep0_) xExpand(ep0_, dir);
  if (ep1_) xExpand(ep1_, dir);
  if (ep2_) xExpand(ep2_, dir);
  if (ep3_) xExpand(ep3_, dir);
  if (et0_) xExpand(et0_, dir);
  if (et1_) xExpand(et1_, dir);
  if (et2_) xExpand(et2_, dir);
  if (et3_) xExpand(et3_, dir);
}

void Worldline::metric(SmartPointer<Metric::Generic> gg) {
  // Unhook from previous metric
  if (metric_) metric_ -> unhook(this);

  // Set the Metric
  metric_=gg;

  // Hook to new metric
  if (metric_) metric_ -> hook(this);

  // Reset integration
  reInit();

}

void Worldline::tell(Gyoto::Hook::Teller* msg) {
  if (msg != metric_) {
    GYOTO_ERROR("Worldline::tell(): wrong Teller");
  }
  reInit();
}

void Worldline::integrator(std::string const &type) {
  if (type=="Legacy") state_ = new IntegState::Legacy(this);
#ifdef GYOTO_HAVE_BOOST_INTEGRATORS
  else state_ = new IntegState::Boost(this, type);
#else
  else GYOTO_ERROR("unrecognized integrator (recompile with boost?)");
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
    GYOTO_ERROR("Worldline::initCoord() requires an 8-element vector");
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
  coord[4] = (i0_<=imax_ && imin_<=i0_)?x0dot_[i0_]:0.;
  coord[5] = x1dot_[i0_];
  coord[6] = x2dot_[i0_];
  coord[7] = x3dot_[i0_];
  return coord;
}

void Worldline::setInitCoord(const double coord[8], int dir,
			     double const Ephi[8], double const Etheta[8]) {
  GYOTO_DEBUG_ARRAY(coord, 8);

  // If dir is not forced and Worldline has never been initialize,
  // make a default direction depending on particle mass
  if (dir==0 && i0_==0 && imin_==1 && imax_==0) dir = getMass() ? 1 : -1;

  switch (dir) {
  case 0: break; // don't move i0_
  case -1: i0_=x_size_-1; break;// integrate backwards
  case 1: i0_=0; // integrate forwards
  }

  imin_=imax_=i0_;
  tau_[i0_]=0.;
  x0_[i0_]=coord[0];
  x1_[i0_]=coord[1];
  x2_[i0_]=coord[2];
  x3_[i0_]=coord[3];
  x0dot_[i0_]=coord[4];
  x1dot_[i0_]=coord[5];
  x2dot_[i0_]=coord[6];
  x3dot_[i0_]=coord[7];
  if (parallel_transport_) {
    ep0_[i0_] = Ephi[0];
    ep1_[i0_] = Ephi[1];
    ep2_[i0_] = Ephi[2];
    ep3_[i0_] = Ephi[3];
    et0_[i0_] = Etheta[0];
    et1_[i0_] = Etheta[1];
    et2_[i0_] = Etheta[2];
    et3_[i0_] = Etheta[3];
  }
  reInit();
}

void Worldline::setInitCoord(const double coord[8], int dir) {
  double const zeroes[4] = {0., 0., 0., 0.};
  setInitCoord(coord, dir, zeroes, zeroes);
}

void Worldline::setInitCoord(double const pos[4], double const v[3], int dir) {
  if (!getMass())
    GYOTO_ERROR("Worldline::setInitCoord(pos, vel) "
	       "only makes sense for massive particles");
  if (!metric_)
    GYOTO_ERROR("Please set metric before calling "
	       "Worldline::setInitCoord(double pos[4], double vel[3])");
  double coord[8]={pos[0], pos[1], pos[2], pos[3],
		   1., v[0], v[1], v[2]};
  metric_ -> normalizeFourVel(coord);
  GYOTO_DEBUG_ARRAY(coord, 8);
  setInitCoord(coord, dir);
}

void Worldline::setInitialCondition(SmartPointer<Metric::Generic> met,
				    const double coord[8],
				    const int dir,
				    double const Ephi[4], double const Etheta[4])
{
  metric(met);
  setInitCoord(coord, dir, Ephi, Etheta);
}

void Worldline::setInitialCondition(SmartPointer<Metric::Generic> met,
				    const double coord[8],
				    const int dir)
{
  metric(met);
  setInitCoord(coord, dir);
}

void Worldline::setPosition(double const pos[4]) {
  double vel[] = {0., 0., 0.};
  setInitCoord(pos, vel);
}

void Worldline::setVelocity(double const vel[3]) {
  state_t coord;
  getInitialCoord(coord);
  setInitCoord(&coord[0], vel);
}



void Worldline::reset() { if (imin_<=imax_) imin_=imax_=i0_; }
void Worldline::reInit() {
  if (imin_ <= imax_) {
    reset();
    state_t coord;
    getInitialCoord(coord);
    GYOTO_DEBUG_ARRAY(coord, 8);
    if (metric_) {
      if ( metric_() -> coordKind() == GYOTO_COORDKIND_SPHERICAL 
	   && x2_[i0_]==0. ) {
	if (verbose() >= GYOTO_SEVERE_VERBOSITY)
	  cerr << "SEVERE: Worldline::reInit(): Kicking particle off z axis\n";
	x2_[i0_]=coord[2]=1e-10;
      }
      metric_ -> setParticleProperties(this,&coord[0]);
    }
  }
  state_ -> init();
}


void Worldline::xStore(size_t ind, state_t const &coord, double tau)
{
  tau_[ind]= tau;
  x0_[ind] = coord[0];
  x1_[ind] = coord[1];
  x2_[ind] = coord[2];
  x3_[ind] = coord[3];
  x0dot_[ind] = coord[4];
  x1dot_[ind] = coord[5];
  x2dot_[ind] = coord[6];
  x3dot_[ind] = coord[7];
  if (parallel_transport_) {
    ep0_[ind] = coord[ 8];
    ep1_[ind] = coord[ 9];
    ep2_[ind] = coord[10];
    ep3_[ind] = coord[11];
    et0_[ind] = coord[12];
    et1_[ind] = coord[13];
    et2_[ind] = coord[14];
    et3_[ind] = coord[15];
  }
}

void Worldline::xFill(double tlim, bool proper) {

  int dir;
  stopcond=0;
  size_t ind;
  
  // Check whether anything needs to be done,
  // Determine direction,
  // Allocate memory.
  double * time_ = proper?tau_:x0_;
  GYOTO_IF_DEBUG
  GYOTO_DEBUG << "x_size_=" << x_size_
	      <<", imin_=" << imin_
	      << ", i0_=" << i0_
	      << ", imax_=" << imax_
	      << ", tlim=" << tlim
	      << ", time_[imin_]=" << time_[imin_]
	      << ", time_[i0_]=" << time_[i0_]
	      << ", time_[imax_]=" << time_[imax_]
	      << ", proper=" << proper
	      << ", (time_==tau_)=" << (time_==tau_)
	      << endl;
  GYOTO_DEBUG << "Possibly extending worldline" << endl;
  GYOTO_ENDIF_DEBUG
  size_t old_x_size=x_size_;
  if (tlim > time_[imax_]) {
    // integrate forward
    dir = 1;
    ind = (imax_==x_size_-1)?xExpand(1):imax_;
  } else if (tlim < time_[imin_]) {
    // integrate backward
    dir = -1;
    ind = (imin_==0)?xExpand(-1):imin_;
  } else return ; // nothing to do
  if (x_size_ != old_x_size) {
    time_ = proper?tau_:x0_;
    old_x_size=x_size_;
  }
  GYOTO_DEBUG << "x_size_=" << x_size_
	      <<", imin_=" << imin_
	      << ", i0_=" << i0_
	      << ", imax_=" << imax_
	      << ", dir=" << dir << ", ind=" << ind << ", tlim=" << tlim
	      << ", time_[imin_]=" << time_[imin_]
	      << ", time_[i0_]=" << time_[i0_]
	      << ", time_[imax_]=" << time_[imax_]
	      << ", proper=" << proper
	      << ", (time_==tau_)=" << (time_==tau_)
	      << endl;
  
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
    GYOTO_ERROR("In Worldline.C Unrecognized mass.");
    //GYOTO_DEBUG<< "of unrecognized mass (!!) particule ....." << endl;
    //equations of geodesics written for a mass=1 star
  }

  state_t coord(parallel_transport_?16:8);
  getCoord(ind, coord);
  double tau=tau_[ind];
  
  GYOTO_DEBUG << "IntegState initialization" << endl;
  
  state_->init(this, coord, dir*delta_);
  //delta_ = initial integration step (defaults to 0.01)
  
  GYOTO_DEBUG << "IntegState initialized" << endl;
  
  size_t mycount=0;// to prevent infinite integration
  
  while (!stopcond) {
    mycount++;
    
    stopcond= state_ -> nextStep(coord, tau);
    
    if (coord[0] == x0_[ind]) { // here, ind denotes previous step
      stopcond=1;
#     if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "time did not evolve, break." << endl;
#     endif
      break;
    }
    if(metric_->isStopCondition(&coord[0])) {
#     if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG << "stopcond set by metric"<<endl;
#     endif
      //coord[0]=1.01*tlim;
      //xStore(ind, coord);
      break;
    }
    
    //if (stopcond && debug()) cout << "stopcond from integrator" << endl;
    if (mycount==maxiter_) {
      stopcond=1;
      Error ( "***WARNING STOP: in Worldline.C unexplained stop !!!" );
    }
    // Check stop condition and whether we need to expand the arrays
    if (dir==1) {
      if ((proper?tau:coord[0]) > tlim) stopcond=1;
      if (ind==x_size_-1) {
	imax_=x_size_-1;
	ind=xExpand(1);
      }
    } else {
      if ((proper?tau:coord[0]) < tlim) {
	stopcond=1;
      }
      if (ind==0) {
	imin_=0;
	ind=xExpand(-1);
      }
    }
    if (x_size_ != old_x_size) {
      time_ = proper?tau_:x0_;
      old_x_size=x_size_;
    }

    // store particle's trajectory for later use
    GYOTO_DEBUG_EXPR(ind);
    ind +=dir;
    xStore(ind, coord, tau);
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

void Worldline::get_tau(double *dest) const
{ memcpy(dest, tau_+imin_, sizeof(double)*(imax_-imin_+1)); }

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
  default: GYOTO_ERROR("in Worldline::get_xyz: Incompatible coordinate kind");
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
    GYOTO_ERROR("in Worldline::get_xyz: unknown coordinate kind");
    x1=x2=x3=NULL; // fix warning
  }

  getCoord(dates, n_dates, x1, x2, x3, x0dot, x1dot, x2dot, x3dot);

  // for (di=0; di<n_dates; ++di) {
  //   double xcur=x1[di], ycur=x2[di], zcur=x3[di], xpcur=x1dot[di]/x0dot[di],
  //     ypcur=x2dot[di]/x0dot[di],
  //     zpcur=x3dot[di]/x0dot[di],
  //     rcur=sqrt(xcur*xcur+ycur*ycur+zcur*zcur),
  //     costhcur=zcur/rcur, sinthcur=sqrt(1-costhcur*costhcur),
  //     cosphcur=xcur/(rcur*sinthcur), sinphcur=ycur/(rcur*sinthcur),
  //     rpcur=xpcur*sinthcur*cosphcur + ypcur*sinphcur*sinthcur + zpcur*costhcur,
  //     phpcur=(ypcur*cosphcur - xpcur*sinphcur)/(rcur*sinthcur),
  //     dtKS_over_dtau=x0dot[di],
  //     gpp=rcur*rcur*sinthcur*sinthcur,
  //     Lcst=gpp*phpcur*dtKS_over_dtau;
  //     cout << "LL= " << Lcst << endl;
  // }
  
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
	    + r * phiprime * sintheta * cosphi;
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
  default: GYOTO_ERROR("in Worldline::get_xyz: unknown coordinate kind");
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
			 double * const x2dot, double * const x3dot,
			 double * ep0, double * ep1, double * ep2, double * ep3,
			 double * et0, double * et1, double * et2, double * et3,
			 double * otime, bool proper)
{
  size_t curl=imin_, curm, curh=imax_;
       // Current indices for the binary search. 
       // To avoid overflows, compute curm as curl+(curh-curl)/2
  size_t di=0; // current date index
  double date; // current date

  double * time_ = proper?tau_:x0_;  // tau_ or x0_
  double * otime_ = proper?x0_:tau_; // x0_ or tau_

  // For the interpolation
  int sz = parallel_transport_?16:8;
  state_t bestl(sz), besth(sz), resl(sz), resh(sz); // i/o for myrk4
  double factl, facth, bestaul, bestauh, restaul, restauh;
  double tausecond, dtaul, dtauh, dtl, dth, Dt, Dtm1, tauprimel, tauprimeh;
  double second, primel, primeh, pos[4], vel[3], tdot;
  int i;
  stringstream ss;
  GYOTO_IF_DEBUG
  GYOTO_DEBUG_EXPR(dates[0]);
  GYOTO_DEBUG_EXPR(time_[imin_]);
  GYOTO_DEBUG_EXPR(time_[imax_]);
  GYOTO_ENDIF_DEBUG

  for (di=0; di<n_dates; ++di) {
    date = dates[di];
    if (date == time_[imax_]) {
      double pos2[8]={0.,0.,x2_[imax_],x3_[imax_],0.,0.,x2dot_[imax_],0.};
      if (metric_->coordKind() == GYOTO_COORDKIND_SPHERICAL)
	checkPhiTheta(pos2);
      if (otime)     otime[di] = otime_[imax_];
      if (x1)       x1[di] =    x1_[imax_];
      if (x2)       x2[di] =   pos2[2];
      if (x3)       x3[di] =   pos2[3];
      if (x0dot) x0dot[di] = x0dot_[imax_];
      if (x1dot) x1dot[di] = x1dot_[imax_];
      if (x2dot) x2dot[di] =   pos2[6];
      if (x3dot) x3dot[di] = x3dot_[imax_];
      if (parallel_transport_) {
	if (ep0)     ep0[di] =   ep0_[imax_];
	if (ep1)     ep1[di] =   ep1_[imax_];
	if (ep2)     ep2[di] =   ep2_[imax_];
	if (ep3)     ep3[di] =   ep3_[imax_];
	if (et0)     et0[di] =   et0_[imax_];
	if (et1)     et1[di] =   et1_[imax_];
	if (et2)     et2[di] =   et2_[imax_];
	if (et3)     et3[di] =   et3_[imax_];
      }
      continue;
    } else if (date > time_[imax_]) {
      GYOTO_DEBUG << "Extending worldline towards future" << endl;
      curl=imax_;    // current imax_
      xFill(date, proper);   // integrate, that changes imax_, tau_, x0_...
      time_ = proper?tau_:x0_;
      otime_ = proper?x0_:tau_;
      curh=imax_;    // new imax_
      if (curl == curh || date > time_[imax_]) {
	ss<<"Worldline::getCoord: can't get coordinates for date="<<date;
	GYOTO_ERROR(ss.str());
      }
    } else if (date < time_[imin_]) {
      GYOTO_DEBUG << "Extending worldline towards past" << endl;
      GYOTO_DEBUG << "imin_=" << imin_ << ", i0_=" << i0_ << ", imax_=" << imax_ << endl;
      curh=x_size_-imin_; // trick if line is expanded during xFill()
      xFill(date, proper);   // integrate, that changes imin_
      time_ = proper?tau_:x0_;
      otime_ = proper?x0_:tau_;
      GYOTO_IF_DEBUG
      GYOTO_DEBUG_THIS << "imin_=" << imin_ << ", i0_=" << i0_ << ", imax_=" << imax_ << endl;
      GYOTO_DEBUG_THIS_EXPR(time_[imin_]);
      GYOTO_DEBUG_THIS_EXPR(time_[imin_+1]);
      GYOTO_DEBUG_THIS_EXPR(time_[imin_+2]);
      GYOTO_DEBUG_THIS_EXPR(time_[imin_+3]);
      GYOTO_DEBUG_THIS_EXPR(time_[imin_+4]);
      GYOTO_DEBUG_THIS_EXPR(time_[imin_+5]);
      GYOTO_DEBUG_THIS_EXPR(time_[imin_+6]);
      GYOTO_DEBUG_THIS_EXPR(time_[imin_+7]);
      GYOTO_DEBUG_THIS_EXPR(time_[imin_+8]);
      GYOTO_DEBUG_THIS_EXPR(time_[imax_-5]);
      GYOTO_DEBUG_THIS_EXPR(time_[imax_-4]);
      GYOTO_DEBUG_THIS_EXPR(time_[imax_-3]);
      GYOTO_DEBUG_THIS_EXPR(time_[imax_-2]);
      GYOTO_DEBUG_THIS_EXPR(time_[imax_-1]);
      GYOTO_DEBUG_THIS_EXPR(time_[imax_]);
      GYOTO_ENDIF_DEBUG
      curh=x_size_-curh;
      curl=imin_;    // new imin_
      if (curl == curh || date < time_[imin_]) {
	ss<<"Worldline::getCoord: can't get coordinates for date="<<date;
	GYOTO_ERROR(ss.str());
      }
    } else if (date >= time_[curh]) {
      curl=curh;
      curh=imax_;
    } else if (date < time_[curl]) {
      curh=curl;
      curl=imin_;
    }

    while (curh-curl>1) {
      curm = curl+(curh-curl)/2;
      if (date >= time_[curm]) curl = curm;
      else curh = curm;
    }

    if (date == time_[curl]) {
      double pos2[8]={0.,0.,x2_[curl],x3_[curl],0.,0.,x2dot_[curl],0.};
      if (metric_->coordKind() == GYOTO_COORDKIND_SPHERICAL)
	checkPhiTheta(pos2);
      if (otime)     otime[di] = otime_[curl];
      if (x1)       x1[di] =    x1_[curl];
      if (x2)       x2[di] =   pos2[2];
      if (x3)       x3[di] =   pos2[3];
      if (x0dot) x0dot[di] = x0dot_[curl];
      if (x1dot) x1dot[di] = x1dot_[curl];
      if (x2dot) x2dot[di] =   pos2[6];
      if (x3dot) x3dot[di] = x3dot_[curl];
      if (parallel_transport_) {
	if (ep0)     ep0[di] =   ep0_[curl];
	if (ep1)     ep1[di] =   ep1_[curl];
	if (ep2)     ep2[di] =   ep2_[curl];
	if (ep3)     ep3[di] =   ep3_[curl];
	if (et0)     et0[di] =   et0_[curl];
	if (et1)     et1[di] =   et1_[curl];
	if (et2)     et2[di] =   et2_[curl];
	if (et3)     et3[di] =   et3_[curl];
      }
      continue;
    }

    // Attempt to get closer to the specified date using the
    // integrator.
    if (proper) {
      dtaul=date-tau_[curl];
      dtauh=date-tau_[curh];
    } else {
      dtl=date-x0_[curl]; dth=date-x0_[curh];
      Dt=(x0_[curh]-x0_[curl]); Dtm1=1./Dt;
      tauprimel=1./x0dot_[curl]; tauprimeh=1./x0dot_[curh];
      tausecond = (tauprimeh-tauprimel)*Dtm1;
      // tauprime(dt)=tauprime(t0)+tausecond*dt
      dtaul=tauprimel*dtl+0.5*tausecond*dtl*dtl;
      dtauh=tauprimeh*dth+0.5*tausecond*dth*dth;
    }

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
    if (parallel_transport_) {
      bestl[8]  =   ep0_[curl];
      bestl[9]  =   ep1_[curl];
      bestl[10] =   ep2_[curl];
      bestl[11] =   ep3_[curl];
      bestl[12] =   et0_[curl];
      bestl[13] =   et1_[curl];
      bestl[14] =   et2_[curl];
      bestl[15] =   et3_[curl];
    }
    bestaul=tau_[curl];
    restaul=bestaul+dtaul;
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
    if (parallel_transport_) {
      besth[8]  =   ep0_[curh];
      besth[9]  =   ep1_[curh];
      besth[10] =   ep2_[curh];
      besth[11] =   ep3_[curh];
      besth[12] =   et0_[curh];
      besth[13] =   et1_[curh];
      besth[14] =   et2_[curh];
      besth[15] =   et3_[curh];
    }

    bestauh=tau_[curh];
    restauh=bestauh+dtauh;
    state_ -> doStep(besth, dtauh, resh);

    if (proper) {
      facth=dtaul/(bestauh-bestaul);
      factl=1.-facth;
      if (otime)     otime[di] = factl*resl[0]+facth*resh[0];
      if (x1)       x1[di] = factl*resl[1]+facth*resh[1];
      if (x2)       x2[di] = factl*resl[2]+facth*resh[2];
      if (x3)       x3[di] = factl*resl[3]+facth*resh[3];
      if (x0dot) x0dot[di] = factl*resl[4]+facth*resh[4];
      if (x1dot) x1dot[di] = factl*resl[5]+facth*resh[5];
      if (x2dot) x2dot[di] = factl*resl[6]+facth*resh[6];
      if (x3dot) x3dot[di] = factl*resl[7]+facth*resh[7];
      if (parallel_transport_) {
	if (ep0)     ep0[di] =   factl*resl[ 8]+facth*resh[ 8];
	if (ep1)     ep1[di] =   factl*resl[ 9]+facth*resh[ 9];
	if (ep2)     ep2[di] =   factl*resl[10]+facth*resh[10];
	if (ep3)     ep3[di] =   factl*resl[11]+facth*resh[11];
	if (et0)     et0[di] =   factl*resl[12]+facth*resh[12];
	if (et1)     et1[di] =   factl*resl[13]+facth*resh[13];
	if (et2)     et2[di] =   factl*resl[14]+facth*resh[14];
	if (et3)     et3[di] =   factl*resl[15]+facth*resh[15];
      }
      continue;
    }

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
      if (resl[0] > bestl[0]) {
	memcpy(&bestl[0], &resl[0], sz*sizeof(double));
	bestaul=restaul;
      }
    } else {
      if (resl[0] < besth[0]) {
	memcpy(&besth[0], &resl[0], sz*sizeof(double));
	bestauh=restaul;
      }
    }

    if (resh[0]<=date) {
      if (resh[0] > bestl[0]) {
	memcpy(&bestl[0], &resh[0], sz*sizeof(double));
	bestaul=restauh;
      }
    } else {
      if (resh[0] < besth[0]) {
	memcpy(&besth[0], &resh[0], sz*sizeof(double));
	bestauh=restauh;
      }
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
      if (otime)     otime[di] = bestaul;
      if (x1)       x1[di] = bestl[1];
      if (x2)       x2[di] = bestl[2];
      if (x3)       x3[di] = bestl[3];
      if (x0dot) x0dot[di] = bestl[4];
      if (x1dot) x1dot[di] = bestl[5];
      if (x2dot) x2dot[di] = bestl[6];
      if (x3dot) x3dot[di] = bestl[7];
      if (parallel_transport_) {
	if (ep0)     ep0[di] =   bestl[ 8];
	if (ep1)     ep1[di] =   bestl[ 9];
	if (ep2)     ep2[di] =   bestl[10];
	if (ep3)     ep3[di] =   bestl[11];
	if (et0)     et0[di] =   bestl[12];
	if (et1)     et1[di] =   bestl[13];
	if (et2)     et2[di] =   bestl[14];
	if (et3)     et3[di] =   bestl[15];
      }
    }
    if (besth[0]==date) {
      if (otime)     otime[di] = bestauh;
      if (x1)       x1[di] = besth[1];
      if (x2)       x2[di] = besth[2];
      if (x3)       x3[di] = besth[3];
      if (x0dot) x0dot[di] = besth[4];
      if (x1dot) x1dot[di] = besth[5];
      if (x2dot) x2dot[di] = besth[6];
      if (x3dot) x3dot[di] = besth[7];
      if (parallel_transport_) {
	if (ep0)     ep0[di] =   besth[ 8];
	if (ep1)     ep1[di] =   besth[ 9];
	if (ep2)     ep2[di] =   besth[10];
	if (ep3)     ep3[di] =   besth[11];
	if (et0)     et0[di] =   besth[12];
	if (et1)     et1[di] =   besth[13];
	if (et2)     et2[di] =   besth[14];
	if (et3)     et3[di] =   besth[15];
      }
    }

    dtl=date-bestl[0]; Dt=besth[0]-bestl[0]; Dtm1=1./Dt;
    facth=dtl*Dtm1; factl=1.-facth;
    if (otime)       otime[di] = bestaul*factl + bestauh*facth;
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
      if (parallel_transport_) GYOTO_ERROR("Parallel transport not implemented for massive particles.");
    } else {
      // Photon: don't be so elaborate, we certainly don't need it... yet
      if (x1)       x1[di] = bestl[ 1]*factl + besth[ 1]*facth;
      if (x2)       x2[di] = bestl[ 2]*factl + besth[ 2]*facth;
      if (x3)       x3[di] = bestl[ 3]*factl + besth[ 3]*facth;
      if (x0dot) x0dot[di] = bestl[ 4]*factl + besth[ 4]*facth;
      if (x1dot) x1dot[di] = bestl[ 5]*factl + besth[ 5]*facth;
      if (x2dot) x2dot[di] = bestl[ 6]*factl + besth[ 6]*facth;
      if (x3dot) x3dot[di] = bestl[ 7]*factl + besth[ 7]*facth;
      if (parallel_transport_) {
	if (ep0)   ep0[di] = bestl[ 8]*factl + besth[ 8]*facth;
	if (ep1)   ep1[di] = bestl[ 9]*factl + besth[ 9]*facth;
	if (ep2)   ep2[di] = bestl[10]*factl + besth[10]*facth;
	if (ep3)   ep3[di] = bestl[11]*factl + besth[11]*facth;
	if (et0)   et0[di] = bestl[12]*factl + besth[12]*facth;
	if (et1)   et1[di] = bestl[13]*factl + besth[13]*facth;
	if (et2)   et2[di] = bestl[14]*factl + besth[14]*facth;
	if (et3)   et3[di] = bestl[15]*factl + besth[15]*facth;
      }
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
  //GYOTO_ERROR("At this point, coordinate conversion is not implemented");
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
    GYOTO_ERROR("Worldline::checkPhiTheta(): should not be called "
	       "with cartesian-like coordinates");
  default:
    GYOTO_ERROR("Worldline::checkPhiTheta(): unknown COORDKIND");
  }
}

void Worldline::get_dot(double *x0dest, double *x1dest, double *x2dest, double *x3dest) const {
  //  if (sysco!=sys_)
  //  GYOTO_ERROR("At this point, coordinate conversion is not implemented");
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
  //  GYOTO_ERROR("At this point, coordinate conversion is not implemented");
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
  default: GYOTO_ERROR("in Worldline::save_xyz: Incompatible coordinate kind");
  }
  
  fichierxyz.close();
}

void Worldline::save_txyz(char * filename, const double t1, const double mass_sun, const double distance_kpc, const string unit, SmartPointer<Screen> sc) {
  xFill(t1, false);
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

void Worldline::integ31(bool integ) {
  state_->integ31(integ);
  reInit();
}
bool Worldline::integ31() const{ return state_->integ31(); }

void Worldline::parallelTransport(bool pt) {
  bool reinit = pt && !parallel_transport_;
  bool deinit = !pt && parallel_transport_;
  parallel_transport_ = pt;
  if (reinit) {
    eAllocate();
    reInit();
  } else if (deinit) {
    eDeallocate();
  }
}
bool Worldline::parallelTransport() const { return parallel_transport_; }

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

double Worldline::maxCrossEqplane() const {return maxCrossEqplane_;}
void Worldline::maxCrossEqplane(double max) {maxCrossEqplane_=max;}

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

void Worldline::constantsOfMotion(std::vector<double> const cstv) {
  setCst(cstv.data(), cstv.size());
}

std::vector<double> Worldline::constantsOfMotion() const {
  std::vector<double> out ;
  out.assign(cst_, cst_+cst_n_);
  return out;
}

void Worldline::getInitialCoord(state_t &coord) const {
  if (imax_<imin_)
    GYOTO_ERROR("Worldline::getInitialCoord(): initial coordinate not set yet");
  getCoord(i0_, coord);
}

void Worldline::getCoord(size_t index, state_t &coord) {
  const_cast<Worldline const *>(this)->getCoord(index, coord);
}

void Worldline::getCoord(size_t index, state_t &coord) const {
  //GYOTO_DEBUG<< "index=" << index << endl;
  //GYOTO_DEBUG<< "x0[index]= " << x1dot_[index] << endl;
  //GYOTO_DEBUG<< "index == " << index << endl;
  if (index<imin_ || index>imax_) {
    cerr << "Indices min curr max= " << imin_ << " " << index << " " << imax_ << endl;
    GYOTO_ERROR("Worldline::getCoord: bad index");
  }
  size_t sz=coord.size();
  if (sz==0) {
    sz = parallel_transport_?16:8;
    coord.resize(sz);
  }
  coord[0] = x0_[index];
  coord[1] = x1_[index];
  coord[2] = x2_[index];
  coord[3] = x3_[index];
  if (sz>=4) {
    coord[4] = x0dot_[index];
    coord[5] = x1dot_[index];
    coord[6] = x2dot_[index];
    coord[7] = x3dot_[index];
    if (parallel_transport_ && sz==16) {
      coord[8] = ep0_[index];
      coord[9] = ep1_[index];
      coord[10] = ep2_[index];
      coord[11] = ep3_[index];
      coord[12] = et0_[index];
      coord[13] = et1_[index];
      coord[14] = et2_[index];
      coord[15] = et3_[index];
    }
  }
}

void Worldline::getCoord(double t, state_t &coord, bool proper) {
  size_t sz=coord.size();
  if (sz==0) {
    sz = parallel_transport_?16:8;
    coord.resize(sz);
  }
  double otime;
  double * x1=&coord[1];
  double * x2=&coord[2];
  double * x3=&coord[3];
  double * x0dot=(sz>4)?&coord[4]:NULL;
  double * x1dot=(sz>4)?&coord[5]:NULL;
  double * x2dot=(sz>4)?&coord[6]:NULL;
  double * x3dot=(sz>4)?&coord[7]:NULL;
  double * ephi0=(sz>8)?&coord[8]:NULL;
  double * ephi1=(sz>8)?&coord[9]:NULL;
  double * ephi2=(sz>8)?&coord[10]:NULL;
  double * ephi3=(sz>8)?&coord[11]:NULL;
  double * etheta0=(sz>8)?&coord[12]:NULL;
  double * etheta1=(sz>8)?&coord[13]:NULL;
  double * etheta2=(sz>8)?&coord[14]:NULL;
  double * etheta3=(sz>8)?&coord[15]:NULL;
  getCoord(&t, 1, x1, x2, x3,                  // x
	   x0dot, x1dot, x2dot, x3dot,         // k
	   ephi0, ephi1, ephi2, ephi3,         // Ephi
	   etheta0, etheta1, etheta2, etheta3, // Etheta
	   &otime, proper);
  coord[0] = proper?otime:t;
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
  default: GYOTO_ERROR("Worldline::getCartesianPos: Incompatible coordinate kind");
  }
}
