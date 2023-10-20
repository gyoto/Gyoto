/**
 *  \file GyotoPhoton.h
 *  \brief A single light ray
 *
 *   Gyoto::Photon is a class for ray-tracing.
 *
 */

/*
    Copyright 2011-2016 Frederic Vincent, Thibaut Paumard

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

#ifndef __GyotoPhoton_H_ 
#define __GyotoPhoton_H_ 

#include "GyotoFunctors.h"

namespace Gyoto{
  class Photon;
  namespace Astrobj { class Generic; }
}

#include <GyotoDefs.h>
#include <GyotoMetric.h>
#include <GyotoScreen.h>
#include <GyotoWorldline.h>
#include <GyotoObject.h>

#include <float.h>
#include <eigen3/Eigen/Dense>

/**
 * \class Gyoto::Photon
 * \brief A null geodesic transporting light
 *
 * This is the central object for ray-tracing. 
 */
class Gyoto::Photon
: public Gyoto::Worldline,
  public Gyoto::SmartPointee,
  public Gyoto::Object
{
  friend class Gyoto::SmartPointer<Gyoto::Photon>;
  // Data : 
  // -----

 protected:
  /// The astronomical target
  /**
   * The (only) possible origin for this Photon.
   */
  SmartPointer<Gyoto::Astrobj::Generic> object_;

  /// Photon's frequency in observer's frame
  /**
   * In Hz, for quantity Emission
   */
  double freq_obs_;

  /// Integrated optical transmission
  /**
   * At Photon::freq_obs_, between current position and observer.
   */
  double transmission_freqobs_;

  /// Integrated optical transmission matrix (polarization)
  /**
   * At Photon::freq_obs_, between current position and observer.
   */
  Eigen::Matrix4d transmissionMatrix_freqobs_;

  /// Observer's spectrometer
  /**
   * Conveying observation frequencies for quantities Spectrum and
   * BinSpectrum.
   */
  SmartPointer<Spectrometer::Generic> spectro_;

  /// Integrated optical transmissions
  /**
   * For each frequency in Photon::spectro_->getMidpoints(), between
   * current position and observer.
   */
  double * transmission_;

  /// Nb of crossings of equatorial plane z=0, theta=pi/2
  int nb_cross_eqplane_;

  /// Integrated optical transmission matrix (polarization)
  /**
   * For each frequency in Photon::spectro_->getMidpoints(), between
   * current position and observer.
   */
  Eigen::Matrix4d * transmissionMatrix_;

  // Constructors - Destructor
  // -------------------------

 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;
  GYOTO_WORLDLINE;
  virtual std::string className() const ; ///< "Photon"
  virtual std::string className_l() const ; ///< "photon"

  /**
   * Allocates the default size (GYOTO_DEFAULT_X_SIZE) for x1, x2 etc.
   */
  Photon() ; ///< Default constructor
  Photon(const Photon& ) ;                ///< Copy constructor
  Photon* clone() const ; ///< Cloner
 protected:  
  Photon(Photon* orig, size_t i0, int dir, double step_max);
  ///< Used by Photon::Refined::Refined()
 public:
  /// Same as Photon() followed by setInitialCondition(SmartPointer<Metric::Generic> gg, SmartPointer<Astrobj::Generic> obj, const double coord[8])
  Photon(SmartPointer<Metric::Generic> gg, SmartPointer<Astrobj::Generic> obj, 
	 double* coord) ;

  /// Same as Photon() followed by setInitialCondition(SmartPointer<Metric::Generic> gg, SmartPointer<Astrobj::Generic> obj, SmartPointer<Screen> screen, double d_alpha, double d_delta)
  Photon(SmartPointer<Metric::Generic> gg, SmartPointer<Astrobj::Generic> obj, 
	 SmartPointer<Screen> screen, double d_alpha, double d_delta);

  virtual ~Photon() ;                        ///< Destructor
  
  virtual double getMass() const ; ///< Return 0.

  /// Set Photon::object_
  void astrobj(SmartPointer<Astrobj::Generic>);
  /// Get Photon::object_
  SmartPointer<Astrobj::Generic> astrobj() const ;

  using Worldline::metric;
  // Set Worldline::metric_ as well as the Astrobj metric
  void metric(SmartPointer<Metric::Generic>);


  /// Set Photon::spectro_
  void spectrometer(SmartPointer<Spectrometer::Generic> spr);
  /// Get Photon::spectro_
  SmartPointer<Spectrometer::Generic> spectrometer() const ;

  /// Set Photon::freq_obs__
  void freqObs(double);
  /// Get Photon::freq_obs__
  double freqObs() const;

  /// Set Photon::nb_cross_eqplane_
  void nb_cross_eqplane(int);
  /// Get Photon::nb_cross_eqplane_
  int nb_cross_eqplane() const;


  // Mutators / assignment
  // ---------------------
 public:
  using Worldline::setInitialCondition;

  /// Set or re-set the initial condition prior to integration.
  /**
   * Set initial condition for this Photon. Several versions of this
   * function, with and without base vectors to parallel-transport,
   * are provided. If polarization is to be treated, use the version
   * with the base vectors.
   *
   * \param gg    Gyoto::SmartPointer to the Gyoto::Metric in this universe;
   * \param obj   Gyoto::SmartPointer to the target Gyoto::Astrobj;
   * \param coord 8-element array containing the initial condition,
   *        i.e. the 4-position and the 4-velocity of the Photon at
   *        the receiving end;
   *
   */
  void setInitialCondition(SmartPointer<Metric::Generic> gg,
			   SmartPointer<Astrobj::Generic> obj,
			   const double coord[8]) ;

  /// Set or re-set the initial condition prior to integration.
  /**
   * Set initial condition for this Photon. Several versions of this
   * function, with and without base vectors to parallel-transport,
   * are provided. If polarization is to be treated, use the version
   * with the base vectors. Initialize those vectors to the two base
   * vectors for the Stokes parameters in the observer screen.
   *
   * \param gg    Gyoto::SmartPointer to the Gyoto::Metric in this universe;
   * \param obj   Gyoto::SmartPointer to the target Gyoto::Astrobj;
   * \param coord 8-element array containing the initial condition,
   *        i.e. the 4-position and the 4-velocity of the Photon at
   *        the receiving end;
   * \param Ephi initial value of base vector to parallel-transport.
   *        Ignored if #parallel_transport_ is false.
   * \param Etheta initial value of base vector to parallel-transport.
   *        Ignored if #parallel_transport_ is false.
   *
   */
  void setInitialCondition(SmartPointer<Metric::Generic> gg,
			   SmartPointer<Astrobj::Generic> obj,
			   const double coord[8],
			   const double Ephi[4],
			   const double Etheta[4]) ;

  /// Set or re-set the initial condition prior to integration.
  /**
   * Set initial condition for this Photon. This version of this
   * function automatically takes care of initializing the base
   * vectors for polarization when needed.
   *
   * \param gg       Gyoto::SmartPointer to the Gyoto::Metric in this universe;
   * \param obj      Gyoto::SmartPointer to the target Gyoto::Astrobj;
   * \param screen   Observer's screen;
   * \param d_alpha  Direction of arrival (RA offset) in radians;
   * \param d_delta  Direction of arrival (Dec offset) in radians.
   */
  void setInitialCondition(SmartPointer<Metric::Generic> gg,
			   SmartPointer<Astrobj::Generic> obj,
			   SmartPointer<Screen> screen,
			   double d_alpha,
			   double d_delta);

  /// Integrate the geodesic
  /**
   * \param[in,out] data Optional Astrobj::Properties to fill with
   * observational quantities.
   * \return 1 if object was hit, else 0.
   */
  int hit(Astrobj::Properties *data=NULL);

  /**
   * \brief Find minimum of photon--object distance
   *
   * Return the minimum of (*object)(this->getCoord())
   * between t1 and t2. The date of this minimum is returned in tmin.
   *
   * \param[in] object
   *             the distance to minimize is given by
   *             object->operator()(). This method is in particular
   *             implemented by the subclasses of Astrobj::Standard.
   * \param[in]  t1   date
   * \param[in]  t2   date
   * \param[out] tmin on output, date correspondig to the minimum
   * \param[in]  threshold stop searching for a minimum if a value <
   *             threshold is found (very often, we just want to find
   *             a date below the threshold, not the accurate
   *             minimum).
   */
  double findMin(Functor::Double_constDoubleArray* object,
		 double t1, double t2, double &tmin,
		 double threshold = DBL_MIN) ;

  /// Find date for which the photon is at a given distance from the object
  /**
   * \param[in] object Object, must implement operator()
   *        (e.g. Astrobj::Standard, ThinDisk::Standard)
   * \param[in] value The value to find
   * \param[in] tinside A date for which
   *        object->Astrobj::operator()(Photon::getCoord()) is < value
   * \param[in,out] toutside On input: a date for which
   *        object->Astrobj::operator()(Photon::getCoord()) is >
   *        value.  on output, (*object)(getCoord(toutside)) is <
   *        value, very close to value. toutside is closer to tinside
   *        on output than on input.
   */
  void findValue(Functor::Double_constDoubleArray* object,
		 double value,
		 double tinside, double &toutside) ;

#ifdef GYOTO_USE_XERCES
  virtual void setParameters(FactoryMessenger *fmp) ;
  static SmartPointer<Photon> Subcontractor(Gyoto::FactoryMessenger*);
#endif

    /* transmission stuff */
 public:
  /// Set transmission to 1 for each channel as well as scalar transmission
  void resetTransmission() ;
  /// Set transmission matrix to identity matrix for each frequency
  void resetTransmissionMatrix();

  /// Get transmission
  /**
   * Get either Photon::transmission_freqobs_ (with i=-1) or
   * Photon::transmission_[i].
   *
   * \param i channel number of the requested frequency, -1 for
   * Photon::freq_obs_.
   */ 
  double getTransmission(size_t i) const ;

  /// Get transmission matrix
  /**
   * Get either Photon::transmission_freqobs_ (with i=-1) or
   * Photon::transmission_[i].
   *
   * \param i channel number of the requested frequency, -1 for
   * Photon::freq_obs_.
   */ 
  Eigen::Matrix4d getTransmissionMatrix(size_t i) const ;

  /// Get maximum transmission;
  /**
   * Get current maximum of all the transmissions,
   * Photon::transmission_freqobs_ or one elements of the
   * Photon::transmission_ array.
   *
   */ 
  double getTransmissionMax() const ;

  /// Get Photon::transmission_
  /**
   * getTansmission()[i] == getTransmission(size_t i)
   */
  double const * getTransmission() const ;
  /// Get Photon::transmissionMatrix_
  /**
   * getTansmissionMatrix()[i] == getTransmissionMatrix(size_t i)
   */
  Eigen::Matrix4d const * getTransmissionMatrix() const ;

  /// Update transmission in a given channel
  /**
   * getTransmission(size_t i) *= t.
   *
   * \param i channel number. -1 for scalar Photon::transmission_freqobs_.
   * \param t transmission of this fluid element.
   */
  virtual void transmit(size_t i, double t);

  /// Update transmission matrix in a given channel
  /**
   * getTransmissionMatrix(size_t i) *= mat.
   *
   * \param i channel number. -1 for Matrix Photon::transmissionMatrix_freqobs_.
   * \param mat transmission Matrix of this fluid element.
   */
  virtual void transmit(size_t i, Eigen::Matrix4d mat);
  virtual void transfer(double * Inu, double * Qnu, double * Unu, double * Vnu, Eigen::Matrix4d * Onu);
  ///< Perform one step of polarized radiative transfert and update transmission matrix

 private:
  /// Allocate Photon::transmission_
  void _allocateTransmission();
  void _allocateTransmissionMatrix();

 public:
  class Refined;

};

/**
 * \class Gyoto::Photon::Refined
 * \brief Refine last step of integration in a Photon
 *
 * The integration step of a Photon's geodesic is adaptive. This is
 * computationally efficient, but sometimes it is necessary to get the
 * position of a Photon with a finer
 * step. Gyoto::ComplexAstrobj::Impact() is a typical use case.
 *
 * A Refined photon is linked to its parent. In particular, care is
 * taken so that the parent's to update the parent's transmissions
 * whenever the Refined transmissions are touched.
 *
 * Don't use this class blindly: what's guaranteed to work is what is
 * used in Gyoto::ComplexAstrobj::Impact().
 *
 * XML description corresponding to this class is &lt;Photon/&gt;. It
 * supports all the parameters supported by the Gyoto::Worldline class
 * plus an optional &lt;Astrobj/&gt; section to attach a instance of a
 * Gyoto::Astrobj::Generic sub-class.
 */
class Gyoto::Photon::Refined : public Gyoto::Photon {
 protected:
  Photon * parent_; ///< Parent Photon.
 public:
  Refined(Photon *parent, size_t i, int dir, double step_max);
  ///< Constructor
  virtual void transmit(size_t i, double t);
  ///< Update transmission both in *this and in *parent_
  virtual void transmit(size_t i, Eigen::Matrix4d mat);
  ///< Update transmission Matrix both in *this and in *parent_
  virtual void transfer(double * Inu, double * Qnu, double * Unu, double * Vnu, Eigen::Matrix4d * Onu);
  ///< Perform one step of polarized radiative transfert and update transmission matrix
};


#endif
