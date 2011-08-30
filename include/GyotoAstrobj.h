/**
 * \file GyotoAstrobj.h
 * \brief Astronomical object
 *
 *  The target of ray-traced Gyoto::Photon
 */

/*
    Copyright 2011 Thibaut Paumard, Frederic Vincent

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


#ifndef __GyotoAstrobj_H_ 
#define __GyotoAstrobj_H_ 

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

namespace Gyoto{
  class Photon;
  class Astrobj;
  class AstrobjProperties;
  class Metric;
}

#include <GyotoDefs.h>
#include <GyotoSmartPointer.h>
#include <GyotoRegister.h>

/**
 * \class Gyoto::Astrobj
 * \brief Astronomical object
 *
 *  Objects which are supposed to be the target of the ray-tracing
 *  code should inherit from this class.
 *
 *  When implementing a new object, you must:
 *    - make sure the object can be loaded from XML by providing a
 *      subcontractor;
 *    - make sure this subcontractor is registerred in the initialization
 *      routine of your plug-in;
 *    - make sure Impact() works (see below).
 *
 *  In addition, you should make sure that your object plays nicely in
 *  the Yorick plug-in, which means:
 *    - implement the copy constructor and the clone() method;
 *    - implement the fillElement method, used for printing and saving to
 *      XML.
 *
 *  There are basically two ways of making Impact() work: either by
 *  providing your own Impact() function, or by implementing a bunch
 *  of lower level, simpler functions which are used by the generic
 *  Astrobj::Impact(). Those lower functions are not pure virtual, but
 *  the default throws a runtime error:
 *    - operator()() yields a distance or potential defining the interior
 *      of the object;
 *    - getVelocity() yields the velocity field of the fluid ;
 *    - processHitQuantities() fills the Spectrum etc. quantities in the
 *      data parameter of Impact().
 *
 *  processHitQuantities() itself is like Impact() in that you have
 *  two choices: either reimplement it or implement a second lot of
 *  small, low-level functions:
 *    - emission();
 *    - integrateEmission();
 *    - transmission().
 *
 */
class Gyoto::Astrobj : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Astrobj>;
  
  // Data : 
  // -----
 protected:

  /**
   * The Metric in this end of the Universe
   */
  SmartPointer<Gyoto::Metric> gg_;


  /**
   * Maximum distance from the center of the coordinate system at
   * which a photon may hit the object.  Child classes may choose to
   * update rmax at all time or to use it to cache the value, for
   * instance when getRmax() is called. External classes (Photons in
   * particular) must use getRmax() to access this information.
   *
   * rmax_set_==1 means that rmax_ was set using setRmax() or the
   * constructor. In this case, getRmax() must always return this
   * value, not recompute it.
   *
   */                                         
  double rmax_; ///< Maximum distance to the center of the coordinate system.

  /**
   * rmax_set_==1 means that rmax_ was set using setRmax() or the
   * constructor. In this case, getRmax() must always return this
   * value, not recompute it.
   *
   * Use unsetRmax() to reset rmax_set_ to 0.
   *
   */                                         
  int rmax_set_; ///< Never recompute rmax: it was externally set

  /**
   * The kind should match the name of the class, e.g. "Star" for a
   * Gyoto::Star.
   */
  std::string kind_; ///< Kind of object (e.g. "Star"...)

  int flag_radtransf_; ///< 1 if radiative transfer inside Astrobj, else 0

  double critical_value_; ///< see operator()(double const coord[4]) const
  double safety_value_; ///< see operator()(double const coord[4]) const

  // Constructors - Destructor
  // -------------------------
 public:
  /**
   *  kind_ =  "Default", rmax_ = 0., rmax_set_ = 0.
   */
  Astrobj(); ///< Default constructor.

  /**
   *  kind_ =  "Default", rmax_ = radmax, rmax_set_ = 1.
   */
  Astrobj(double radmax); ///< Set rmax in constructor.

  /**
   *  kind_ =  kind, rmax_ = 0., rmax_set_ = 0.
   */
  Astrobj(std::string kind); ///< Set kind in constructor.

  Astrobj(const Astrobj& ) ; ///< Copy constructor.
  virtual Astrobj* clone() const; ///< "Virtual" copy constructor
  
  virtual ~Astrobj() ; ///< Destructor: does nothing.

  // Accessors
  // ---------
 public:
  /**
   * Get the Metric
   */
  virtual SmartPointer<Metric> getMetric() const;

  /**
   * Set the Metric
   */
  virtual void setMetric(SmartPointer<Metric>) ;

  /**
   *  Get maximal distance from center of coordinate system at which a
   *  Photon may hit the object.
   *  
   *  Child classes may use the rmax_ member to cache this value.
   *
   *  It can also be set using setRmax(). If setRmax has been used
   *  to set rmax_, getRmax() must not recompute it.
   */
  virtual double getRmax(); ///< Get maximal distance from center of coordinate system

  const std::string getKind() const; ///< Get the kind of the Astrobj (e.g. "Star")

  /**
   *  Set maximal distance from center of coordinate system at which a
   *  Photon may hit the object.
   *  
   *  Side effect: set rmax_set_ to 1.
   */
  virtual void setRmax(double val); ///< Set maximal distance from center of coordinate system

  /**
   * getRmax() will then be free to recompute rmax_. Astrobjs
   * which cache rmax_ may need to update it when unsetRmax() is
   * called.
   */
  virtual void unsetRmax() ; ///< Set rmax_set_ to 0.

  void setFlag_radtransf(int flag);
  int getFlag_radtransf() const ;

  /**
   * Return a Gyoto::Quantity_t suitable as input to
   * Gyoto::Scenery::setRequestedQuantities() to set de default
   * quantities to compute for this object. The default of these
   * defaults GYOTO_QUANTITY_INTENSITY.
   */
  virtual Quantity_t getDefaultQuantities();

  //XML I/O
 public:
#ifdef GYOTO_USE_XERCES
  /**
   * Astrobj implementations should impement fillElement to save their
   * parameters to XML and call the generic implementation to save
   * generic parts such as Flag_radtrans: Astrobj::fillElement(fmp).
   */

  virtual void fillElement(factoryMessenger *fmp) const ;
                                             /// < called from Factory
  void setGenericParameter(std::string name, std::string content) ;
  ///< To be called by fillElement()

#endif
  
  // Outputs
  // -------
 public:
  /**
   * Impact() checks whether a Photon impacts the object between two
   * integration steps of the photon's trajectory (those two steps are
   * photon->getCoord(index, coord1) and photon->getCoord(index+1,
   * coord2)). Impact returns 1 if the photon impacts the object
   * between these two steps, else 0.
   *
   * Impact will compute observable properties on demand: if the data
   * pointer is non-NULL, the object will look in it for pointers to
   * properties which apply to its kind. If a pointer to a property
   * known to this object is present, then the property is computed
   * and store at the pointed-to adress. For instance, all objects
   * know the "intensity" property. If data->intensity != NULL, the
   * instensity is computed and stored in *data->intensity.
   *
   * \param ph   Gyoto::Photon aimed at the object;
   * \param index    Index of the last photon step;
   * \param data     Pointer to a structure to hold the observables at impact.
   *
   * \return 1 if impact, 0 if not.
   */
  virtual int Impact(Gyoto::Photon* ph, size_t index,
		     AstrobjProperties *data=NULL)  ;
  ///< does a photon at these coordinates impact the object?


  /**
   * A potential, distance, or whatever function such that
   * operator()(double coord[4]) < critical_value_ if and only if
   * coord is inside the object. This function is used by the default
   * implmenetation of Impact(). If Impact() is overloaded, it is not
   * necessary to overload operator()(double coord[4]). The default
   * implementation throws an error.
   */
  virtual double operator()(double const coord[4]) ;
  
 protected:
  /*
    THOSE ARE NOT PART OF THE API.

    THEY _MAY_ BE OVERLOADED AND CALLED BY IMPACT OR CALL ONE-ANOTHER
   */

  /**
   * Not part of the API. This function is not pure virtual since it
   * doesn't need to be implemented if Impact() is overloaded. The
   * generic implementation throws a runtime error.
   *
   * Used by the generic Impact().
   *
   * Fill vel with the 4-vector velocity of the fluid at 4-position pos.
   *
   * \param pos input, 4-position at which to compute velocity;
   * \param vel output, 4-velocity at pos.
   */
  virtual void getVelocity(double const pos[4], double vel[4]) ;

  /**
   * Not part of the API.
   *
   * processHitQuantities fills the requested data in Impact. To use
   * it, you need to call it in the Impact() method for you object in
   * case of hit. It will fill Redshift, Intensity, Spectrum,
   * BinSpectrum.
   *
   * You can overload it for your Astrobj. The generic implementation
   * calls emission() below.
   */
  virtual void processHitQuantities(Photon* ph, double* coord_ph_hit,
				    double* coord_obj_hit, double dt,
				    AstrobjProperties* data) const;

  /**
   * Not part of the API, called by the default implementation for
   * processHitQuantities().
   *
   * emission() computes the intensity I_nu emitted by the small
   * volume of length dsem. It should take self-absorption along dsem
   * into account.
   *
   * Reminder :
   *  - intensity = I_nu [erg cm^-2 s^-1 ster^-1 Hz^-1];
   *  - invariant intensity = I_nu/nu^3, which has the same value in any frame;
   *  - emission coefficient = j_nu [erg cm^-3 s^-1 ster^-1 Hz^-1] ,
   *            defined by dI_nu = j_nu*ds, where ds is the distance
   *            travelled by the photon inside the object;
   *  - invariant emission coef = j_nu/nu^2, which has the same value
   *            in any frame.
   *
   * The equation used for radiative transfer (without absorption) is:
   *                      d(I_nu/nu^3)/dlambda = (j_nu/nu^2)  [*]
   *  where lambda is the integration parameter along the null geodesic.
   *
   *      NB: Let us consider a particular observer, 
   *          with nu being the frequency measured by this observer,
   *          and ds being the proper distance (as measured by the observer) 
   *          that the photon travels as it moves
   *          from lambda to lambda+dlambda along its geodesic.
   *          Then it can be shown that :
   *                          dlambda = ds/nu
   *          This shows that Eq. [*] is homogeneous.
   *
   * \param nu_em Frequency at emission
   * \param dsem length over which to integrate inside the object
   * \param coord_ph Photon coordinate
   * \param coord_obj Emitter coordinate at current photon position
   */
  virtual double emission(double nu_em, double dsem, double coord_ph[8],
			  double coord_obj[8]=NULL)
    const =0; ///< INVARIANT emission j_{\nu}/\nu^{2}

  virtual double integrateEmission(double nu1, double nu2, double dsem,
				   double c_ph[8], double c_obj[8]=NULL) const;
    ///< \sum_nu1^nu2 I_nu dnu (or j_nu)

  /**
   * Not part of the API, called by the default implementation for
   * processHitQuantities().
   *
   * transmission() computes the transmission of this fluid element or
   * 0 if optically thick.
   *
   * \param nuem frequency in the fluid's frame
   * \param coord Photon coordinate
   * \param dsem geometrical length in geometrical units
   */
  virtual double transmission(double nuem, double dsem, double coord[8]) const ;
     ///< Transmission: exp( \alpha_{\nu} * dsem )
  
  // Display
  //friend std::ostream& operator<<(std::ostream& , const Astrobj& ) ;

  ///// REGISTER STUFF ///////
  public:
  typedef Gyoto::SmartPointer<Gyoto::Astrobj>
    Subcontractor_t(Gyoto::factoryMessenger*);

  static Gyoto::Register::Entry * Register_;
  static void initRegister(); 
  static void Register(std::string name, Gyoto::Astrobj::Subcontractor_t* scp);
  static Gyoto::Astrobj::Subcontractor_t* getSubcontractor(std::string name);
};

/**
 * \class Gyoto::AstrobjProperties
 * \brief Observable properties of an Astronomical object
 *
 *  The sort of properties one wants to measure on a ray-traced Gyoto::Photon which hits a Gyoto::Astrobj. Not all Astrobj are able to fill all of these properties.
 *
 */
class Gyoto::AstrobjProperties : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::AstrobjProperties>;
 public:
  double *intensity; ///< Apparent intensity (takes beaming into account); 
  double *time; ///< Date of impact (= date of emission of the photon);
  double *distance; ///< Behaves like the square of the closest distance between Photon and Astrobj (but not exactly that). Initialize it to DBL_MAX from float.h.;
  double * first_dmin; ///< first local minimum in distance from object
  int first_dmin_found; ///< first_dmin will be set to the first local minimum and first_dmin_found will be set to 1 if a local minimum in distance is found. Initialize it to 0.
  double *redshift; ///< redshift factor nuobs/nuem (necessary for emission lines computation)
  double *rimpact; ///< radial coordinate at impact (necessary for emission lines computation)
  double *spectrum; ///< I_nu (nu) (observed specific intensity)
  double *binspectrum; ///< I_nu1^nu2, the integral of I_nu over each spectral channel (i.e. what a spectrometer would measure)
  int offset; ///< spectra elements are separated by offset doubles in memory. In other words, the ith spectral element is a spectrum[i*offset].
  double *x, *y, *z; ///< Cartesian coordinates of the Photon at impact;
  double *user1, *user2, *user3, *user4, *user5; ///< Quantities specific to Astrobj
 public:
  AstrobjProperties(); ///< Default constructor (everything is set to NULL);
  AstrobjProperties (double*, double*); ///<< set intensity and time pointers.
  void init(size_t nbnuobs=0);
  AstrobjProperties operator++();
};

#endif
