/*
    Copyright 2011-2012, 2014-2016, 2018-2020 Thibaut Paumard & Frédéric Vincent

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
#include "GyotoError.h"
#include "GyotoPhoton.h"
#include <cmath>
#include <cstdlib>
#include <clocale>

#include "GyotoScenery.h"
#include "GyotoSpectrum.h"
#include "GyotoMetric.h"
#include "GyotoAstrobj.h"
#include "GyotoSpectrometer.h"
#include "GyotoScreen.h"

using namespace Gyoto;
using namespace std;

static int gyoto_debug=GYOTO_DEFAULT_DEBUG_MODE;
#if GYOTO_DEFAULT_DEBUG_MODE
static int gyoto_verbosity=GYOTO_DEBUG_VERBOSITY;
static int gyoto_prev_verbosity=GYOTO_DEBUG_VERBOSITY;
#else
static int gyoto_verbosity=GYOTO_DEFAULT_VERBOSITY;
static int gyoto_prev_verbosity=GYOTO_DEBUG_VERBOSITY;
#endif

#if defined GYOTO_USE_ARBLIB
# include <acb_hypgeom.h>
#elif defined GYOTO_USE_AEAE
# include <complex>
# include <iostream>
# define SIGN(a) (((a) < 0) ? (-1) : (1))
# include "complex_functions.H"
# include "hyp_2F1.cpp"
#endif

void Gyoto::debug(int mode) {
  if (mode != gyoto_debug) {
    if (mode) {
      gyoto_prev_verbosity=verbose();
      verbose(GYOTO_DEBUG_VERBOSITY);
    } else {
      verbose(gyoto_prev_verbosity);
    }
    gyoto_debug=mode;
  } 
}
int Gyoto::debug() { return gyoto_debug; }

void Gyoto::verbose(int mode) { gyoto_verbosity=mode; }
int Gyoto::verbose() { return gyoto_verbosity; }

void Gyoto::convert(double * const x, const size_t nelem, const double mass_sun, const double distance_kpc, const string unit) {
  /// Convert lengths
  
  double distance = distance_kpc*GYOTO_KPC;  // m
  double fact   = mass_sun * GYOTO_SUN_MASS * GYOTO_G_OVER_C_SQUARE; // m
  size_t i =0;


  if (!unit.compare("geometrical"))     return ;
  else if (!unit.compare("m"))          ;
  else if (!unit.compare("km"))         fact *=  1e-3 ;
  else if (!unit.compare("sun radius")) fact *=  1.      / GYOTO_SUN_RADIUS;
  else if (!unit.compare("rad"))        fact *=  1.      / (distance);
  else if (!unit.compare("degree"))     fact *=  180.    / (distance*M_PI);
  else if (!unit.compare("arcmin"))     fact *=  1.08e4  / (distance*M_PI);
  else if (!unit.compare("arcsec"))     fact *=  6.48e5  / (distance*M_PI);
  else if (!unit.compare("mas"))        fact *=  6.48e8  / (distance*M_PI);
  else if (!unit.compare("uas"))        fact *=  6.48e11 / (distance*M_PI);
  else GYOTO_ERROR("Unknown unit.");

  for (i=0; i<nelem; ++i) x[i] *= fact ;

}

double Gyoto::atof(const char * str)
{
  GYOTO_DEBUG << "Gyoto::atof(\"" << str << "\")";
  ptrdiff_t offset=0;
  while (isspace(str[offset])) ++offset;
  bool positive=true;
  double retval=0.;
  if (str[offset] == '-') {
    positive=false;
    ++offset;
  }
  if (str[offset++]=='D' && str[offset++]=='B' && str[offset++]=='L' &&
      str[offset++]=='_' && str[offset++]=='M') {
    if (str[offset]=='A' && str[offset+1]=='X') {
      if (positive) retval = DBL_MAX;
      else retval = -DBL_MAX;
    } else if (str[offset]=='I' && str[offset+1]=='N') {
      if (positive) retval = DBL_MIN;
      else retval = -DBL_MIN;
    } else GYOTO_ERROR("unrecognize double representation");
  } else {
    std::string loc(setlocale(LC_NUMERIC, NULL));
    setlocale(LC_NUMERIC, "C");
    retval = std::atof(str);
    setlocale(LC_NUMERIC, loc.c_str());
  }

  GYOTO_DEBUG << "==" << retval << endl;

  return retval;
}

void Gyoto::help(std::string class_name) {
  std::vector<std::string> plugins;
  if (class_name.substr(0, 7)=="Gyoto::")
    class_name=class_name.substr(7);
  if (class_name=="Scenery") {Scenery().help(); return;}
  if (class_name=="Screen") {Screen().help(); return;}
  if (class_name=="Photon") {Photon().help(); return;}
  size_t pos=class_name.find("::");
  if (pos==0 || pos+2==class_name.size())
    GYOTO_ERROR("Not a valid class name: "+class_name);
  if (pos > 0 && pos != string::npos) {
    string nspace = class_name.substr(0, pos);
    class_name = class_name.substr(pos+2);
    if (nspace=="Astrobj") {
      (*Astrobj::getSubcontractor(class_name, plugins))
	(NULL, plugins)->help();
      return;
    }
    if (nspace=="Metric") {
      (*Metric::getSubcontractor(class_name, plugins))
	(NULL, plugins)->help();
      return;
    }
    if (nspace=="Spectrum") {
      (*Spectrum::getSubcontractor(class_name, plugins))
	(NULL, plugins)->help();
      return;
    }
    if (nspace=="Spectrometer") {
      (*Spectrometer::getSubcontractor(class_name, plugins))
	(NULL, plugins)->help();
      return;
    }
    GYOTO_ERROR("Unrecognized namespace: "+nspace);
  }
  GYOTO_ERROR("Help string not implemented (yet) for "+class_name);
}

std::vector<std::string> Gyoto::split(std::string const &src, std::string const &delim) {
  std::vector<std::string> res;
  size_t pos=0, fpos=0, sz=src.length();
  std::string tmp("");
  while (fpos != string::npos && pos < sz) {
    fpos = src.find_first_of(delim, pos);
    if (fpos==pos) {++pos; continue;}
    res.push_back(src.substr(pos, fpos-pos));
    pos = fpos+1;
  }
  return res;
}

// Bessel functions
double Gyoto::bessi0(double xx) {
  double ax,ans,y;
  if((ax=fabs(xx))< 3.75){
    y=xx/3.75;
    y*=y;
    ans=1.0+y*(3.5156229+y*(3.0899424+y*(1.2067492
					 +y*(0.2659732
					     +y*(0.360768e-1
						 +y*0.45813e-2)))));
  }else{
    y=3.75/ax;
    ans=(exp(ax)/sqrt(ax))
      *(0.39894228
	+y*(0.1328592e-1
	    +y*(0.225319e-2
		+y*(-0.157565e-2
		    +y*(0.916281e-2
			+y*(-0.2057706e-1
			    +y*(0.2635537e-1
				+y*(-0.1647633e-1
				    +y*0.392377e-2))))))));
  }
  return ans;
}
double Gyoto::bessk0(double xx) {
  double ans,y;
  if(xx<=2.0){
    y=xx*xx/4.0;
    ans=(-log(xx/2.0)*bessi0(xx))
      +(-0.57721566
	+y*(0.42278420
	    +y*(0.23069756+y*(0.3488590e-1
			      +y*(0.262698e-2
				  +y*(0.10750e-3+y*0.74e-5))))));
  }else{
    y=2.0/xx;
    ans=(exp(-xx)/sqrt(xx))*(1.25331414
			     +y*(-0.7832358e-1
				 +y*(0.2189568e-1
				     +y*(-0.1062446e-1
					 +y*(0.587872e-2
					     +y*(-0.251540e-2
						 +y*0.53208e-3))))));
  }
  return ans;
}
double Gyoto::bessi1(double xx) {
  double ax,ans,y;
  if((ax=fabs(xx))< 3.75){
    y=xx/3.75;
    y*=y;
    ans=ax*(0.5+y*(0.87890594
		   +y*(0.51498869
		       +y*(0.15084934
			   +y*(0.2658733e-1
			       +y*(0.301532e-2+y*0.32411e-3))))));
  }else{
    y=3.75/ax;
    ans=0.2282967e-1+y*(-0.2895312e-1+y*(0.1787654e-1
					 -y*0.420059e-2));
    ans=0.39894228+y*(-0.3988024e-1+y*(-0.362018e-2
				       +y*(0.163801e-2
					   +y*(-0.1031555e-1+y*ans))));
    ans*=(exp(ax)/sqrt(ax));
  }
  return xx<0.0 ? -ans : ans;
}
double Gyoto::bessk1(double xx) {
  double yy,ans;
  if(xx<=2.0){
    yy=xx*xx/4.0;
    ans=(log(xx/2.0)*bessi1(xx))
      +(1.0/xx)*(1.0+yy*(0.15443144
			 +yy*(-0.67278579
			      +yy*(-0.18156897
				   +yy*(-0.1919402e-1
					+yy*(-0.110404e-2
					     +yy*(-0.4686e-4)))))));
  }else{
    yy=2.0/xx;
    ans=(exp(-xx)/sqrt(xx))*(1.25331414
			     +yy*(0.23498619
				  +yy*(-0.3655620e-1
				       +yy*(0.1504268e-1
					    +yy*(-0.780353e-2
						 +yy*(0.325614e-2
						      +yy*(-0.68245e-3)))))));
  }
  return ans;
}
double Gyoto::bessk(int nn,double xx) {
  double bk,bkm,bkp,tox;
  if(nn< 2) GYOTO_ERROR("In Utils::besselk n>2!");
  tox=2.0/xx;
  bkm=bessk0(xx);
  bk=bessk1(xx);
  for(int j=1;j<nn;j++){
    bkp=bkm+j*tox*bk;
    bkm=bk;
    bk=bkp;
  }
  return bk;
}
// End Bessel functions

double Gyoto::hypergeom (double kappaIndex, double thetae) {
#if defined GYOTO_USE_ARBLIB
  // See documentation: http://arblib.org/acb_hypgeom.html#c.acb_hypgeom_2f1
  acb_t FF, aa, bb, cc, zed;
  acb_init(FF);
  acb_init(aa);
  acb_init(bb);
  acb_init(cc);
  acb_init(zed);
  acb_set_d_d(aa,   kappaIndex-1./3.,  0.);
  acb_set_d_d(bb,   kappaIndex+1.,     0.);
  acb_set_d_d(cc,   kappaIndex+2./3.,  0.);
  acb_set_d_d(zed, -kappaIndex*thetae, 0.);
  slong prec=53; // 53 for double precision
  acb_hypgeom_2f1(FF, aa, bb, cc, zed, ACB_HYPGEOM_2F1_AC, prec);
  double hypergeom = arf_get_d(&acb_realref(FF)->mid, ARF_RND_NEAR);
  // uncertainty
  // double rad = mag_get_d(&acb_realref(FF)->rad);
  acb_clear(FF);
  acb_clear(aa);
  acb_clear(bb);
  acb_clear(cc);
  acb_clear(zed);
  return hypergeom;
#elif defined GYOTO_USE_AEAE
  complex<double> aa=kappaIndex-1./3., bb=kappaIndex+1.,
    cc=kappaIndex+2./3., zed=-kappaIndex*thetae;
  return hyp_2F1(aa,bb,cc,zed).real();
#else
  GYOTO_ERROR("Utils::_hypergeom() is not functional, please recompile Gyoto with either ARBLIB or AEAE");
  return 0.;
#endif
}

// Coordinate transforms
void Gyoto::cartesianToSpherical(double const cpos[3], double spos[3]) {
  spos[0]=sqrt(cpos[0]*cpos[0]+cpos[1]*cpos[1]+cpos[2]*cpos[2]);
  spos[1]=acos(cpos[2]/spos[0]);
  spos[2]=atan2(cpos[1],cpos[0]);
}

void Gyoto::sphericalToCartesian(double const spos[3], double cpos[3]) {
  double c1, s1; sincos(spos[1], &s1, &c1);
  double c2, s2; sincos(spos[2], &s2, &c2);
  cpos[0] = spos[0]*s1*c2;
  cpos[1] = spos[0]*s1*s2;
  cpos[2] = spos[0]*c1;
}

// Matrix inversion
void Gyoto::matrix4Invert(double Am1[4][4], double const A[4][4]) {
  // Invert 4×4 matrix using Gauss pivot
  double tmp[4][4];
  int i, j, jj, jl;
  double fact;

  // Initialize Am1 as identity matrix
  // copy A into tmp
  for (j=0; j<4; ++j) {
    for (i=0; i<4; ++i) {
      Am1[i][j] = (i==j);
      tmp[i][j] = A[i][j];
    }
  }

  // Turn tmp in upper diagonal matrix with ones on the diagonal
  for (j=0; j<4; ++j) {
    // exchange row j with later row with largest coeff in column j
    // (to deal with case of zero or small value on the diagonal
    fact=tmp[j][j];
    jl=j;
    for (jj=j+1; jj<4; ++jj) {
      if (fabs(tmp[j][jj])>fabs(fact)) {
	jl=jj;
	fact=tmp[j][jj];
      }
    }
    fact=1./fact;
    if (jl!=j) {
      double c;
      for (i=j; i<4; ++i) {
	c=tmp[i][jl];
	tmp[i][jl]=tmp[i][j];
	tmp[i][j]=c*fact;
      }
      for (i=0; i<4; ++i) {
	c=Am1[i][jl];
	Am1[i][jl]=Am1[i][j];
	Am1[i][j]=c*fact;
      }
    } else {
      for (i=j; i<4; ++i) tmp[i][j] *= fact;
      for (i=0; i<4; ++i) Am1[i][j] *= fact;
    }
    // subtract this row times of factor to each later row
    for (jj = j+1; jj < 4; ++jj) {
      fact = tmp[j][jj];
      for (i=j+1; i<4; ++i) tmp[i][jj] -= fact*tmp[i][j];
      for (i=0; i<4; ++i) Am1[i][jj] -= fact*Am1[i][j];
    }
  }

  // Cancel upper triangle to get finally the identity matrix
  for (j=3; j>=0; --j) {
    for (jj = j-1; jj >= 0; --jj) {
      fact = tmp[j][jj];
      for (i=0; i<4; ++i) tmp[i][jj] -= fact*tmp[i][j];
      for (i=0; i<4; ++i) Am1[i][jj] -= fact*Am1[i][j];
    }
  }
}


void Gyoto::matrix4CircularInvert(double Am1[4][4], double const A[4][4]) {
  // Works for a metric where gtr, gttheta and grtheta are 0 (and symmetrical...)
  double a=A[0][0], b=A[1][1], c=A[2][2], d=A[3][3], t=A[0][3];
  double t2=t*t;
  double X=d-t2/a;
  double aX=a*X;

  Am1[0][0]=(aX+t2)/(a*aX);
  Am1[1][1]=1./b;
  Am1[2][2]=1./c;
  Am1[3][3]=1/X;
  Am1[0][3]=Am1[3][0]=-t/aX;
  Am1[0][1]=Am1[1][0]=0.;
  Am1[0][2]=Am1[2][0]=0.;
  Am1[1][2]=Am1[2][1]=0.;
  Am1[1][3]=Am1[3][1]=0.;
  Am1[2][3]=Am1[3][2]=0.;
}
