/**
 * \file GyotoScreen.h
 * \brief Description of the observer screen
 * 
 */

/*
    Copyright 2011-2013 Thibaut Paumard, Frederic Vincent

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

#ifndef __GyotoScreen_H_
#define __GyotoScreen_H_ 

#include <iostream>
#include <fstream>
#include <string>

namespace Gyoto {
  class Screen;
}

#include <GyotoDefs.h>
#include <GyotoUtils.h>
#include <GyotoSmartPointer.h>
#include <GyotoMetric.h>
#include <GyotoSpectrometer.h>

/**
 * \class Gyoto::Screen
 * \brief The camera with which the Astrobj is observed
 *
 * In the observer-centric point-of-view, the center of the Metric's
 * coordinate system is positioned relatively to the observing Screen
 * using three Euler angles and the distance (in meters). The three
 * Euler angles are:
 *   - position angle of the line of nodes (North of East);
 *   - inclination (0 = face-on);
 *   - argument of the X axis of the Metric's coordinate system.
 *  We use the z-x-z convention.
 *  See http://en.wikipedia.org/wiki/Euler_angles
 *
 * In addition, the Screen conveys:
 *   - the observing date (in geometrical units, but expect it to
 *     change to seconds in a future version);
 *   - the field-of-view of the image;
 *   - the resolution of the camera: number of pixels on each side
 *     (the camera is square);
 *   - the observing frequency.
 *
 * The scalar FreqObs defines the observing frequency for Scenery
 * quantity Intensity.
 *
 * Likewise, a Gyoto::Spectrometer defines for which frequencies
 * spectra are computed (when the Quantity Spectrum is requested in
 * the Scenery).
 * 
 * For the sake of theoreticians, there is an alternate way of
 * specifying the relative position of the Screen and Metric, by
 * specifying the 4-coordinates of the Screen in the Metric's
 * coordinate system (in that case, eerything is specified in
 * geometrical units).
 *
 * So an XML stanza for a Screen may look like that:
 * \code
 *  <Screen>
 *    <Time>       1000.      </Time>
 *    <FieldOfView>   0.3141592653589793 </FieldOfView>
 *    <Resolution>  128       </Resolution>
 *    <Distance>      1e30    </Distance>
 *    <PALN>          3.14159 </PALN>
 *    <Inclination>   2.0944  </Inclination>
 *    <Argument>     -2.0944  </Argument>
 *    <Spectrometer kind="freqlog" nsamples="10"> 17. 23. </Spectrometer> 
 *    <FreqObs>       1e20    </FreqObs>
 *  </Screen>
 * \endcode
 *
 * or like that:
 *
 * \code
 *  <Screen>
 *    <Position> 1000. 1000. 0.15. 0.</Position>
 *    <FieldOfView>   0.3141592653589793 </FieldOfView>
 *    <Resolution>  128 </Resolution>
 *    <Spectrometer kind="freqlog" nsamples="10"> 17. 23. </Spectrometer> 
 *    <FreqObs>       1e20    </FreqObs>
 *  </Screen>
 * \endcode
 *
 *
 * Units can be specified using the unit attribute in the XML file,
 * for instance:
 * 
 * \code
 *   <Distance unit="kpc"> 8 </Distance>
 * \endcode
 *
 * Possible units are (with [] noting the default):
 *  - distance: [m], geometrical, cm, km, AU, ly, pc, kpc, Mpc;
 *  - PALN, inclination, argument: [rad], deg.
 *  - frequency: [Hz], µm, GeV...
 *
 * When the distance is really large and most of the ray-tracing would
 * happen de facto in flat space, the camera is transported to a
 * location at a reasonable distance from the metric and the images
 * are scaled accordingly. The default value for this distance should
 * be fine, but it can be customized using the "dmax" attribute of the
 * "Distance" element. "dmax" is always expressed in geometrical
 * units:
 *
 * \code
 *    <Distance unit="kpc" dmax="1e7"> 8 </Distance>
 * \endcode
 *
 * Symptoms when dmax is too large include pixelization of the image
 * (neighbouring photons are numerically identical) and other
 * numerical overflows. dmax is too small when it is apparent that
 * changing it yields projection effects. dmax must be large compared
 * to rmax in the Astrobj and ideally, changing it by an order of
 * magnitude should not yield significant changes in the ray-traced
 * image.
 *
 */
class Gyoto::Screen : protected Gyoto::SmartPointee {
  friend class Gyoto::SmartPointer<Gyoto::Screen>;

 private:
  double tobs_; ///< Observing date in s
  double fov_;  ///< Field-of-view in rad
  //  double tmin_;
  size_t npix_; ///< Resolution in pixels

  double distance_; ///< Distance to the observer in m
  double dmax_; ///< Maximum distance from which the photons are launched (geometrical units) 

  int anglekind_; ///< Screen angles kind (0: equatorial, 1: spherical)
  
  /**
   * The angles are position angle of the line of nodes (North of
   * East), inclination (0 = face-on), argument of X axis. We use the
   * z-x-z convention. See http://en.wikipedia.org/wiki/Euler_angles
   */
  double euler_[3]; ///< Euler angles
  double ex_[3]; ///< Sky coordinate of base X vector
  double ey_[3]; ///< Sky coordinate of base Y vector
  double ez_[3]; ///< Sky coordinate of base Z vector

  double fourvel_[4]; ///< Observer's 4-velocity
  double screen1_[4]; ///< Screen e1 vector
  double screen2_[4]; ///< Screen e2 vector
  double screen3_[4]; ///< Screen e3 vector (normal)

  double alpha0_; ///< Screen orientation (0,0) is right towards the BH
  double delta0_; ///< Screen orientation (0,0) is right towards the BH
  SmartPointer<Metric::Generic> gg_; ///< The Metric in this end of the Universe

  /**
   * \brief Gyoto::Spectrometer::Generic subclass instance used for quantities Spectrum and BinSpectrum
   */
  SmartPointer<Spectrometer::Generic> spectro_;

  /**
   * \brief Frequency at which the observer observes
   *
   * For the quantity Intensity
   */
  double freq_obs_;

 public:
   
  // Constructors - Destructor
  // -------------------------
  Screen() ; ///< Default constructor
  Screen(const Screen& ) ;                ///< Copy constructor
  Screen * clone() const; ///< Cloner

  virtual ~Screen() ;                        ///< Destructor
  
  // Mutators / assignment
  // ---------------------

  /// Set inclination etc.
  void setProjection(const double paln,
		     const double inclination,
		     const double argument);
  /// Set distance, inclination etc.
  void setProjection(const double distance,
		     const double paln,
		     const double inclination,
		     const double argument);

  /// Set distance from observer
  /**
   * \param dist Distance in meters.
   */
  void setDistance(double dist);

  /// Set ray-tracing maximum distance
  /**
   * \param dist Distance in geometrical units.
   */
  void setDmax(double dist);

  /// Set distance from observer
  /**
   * \param dist the distance expressed in the specified unit;
   * \param unit convertible to meters
   */
  void setDistance(double dist, const std::string unit);

  /// Set inclination relative to line-of-sight
  /**
   * Inclination of z-axis relative to line-of-sight, or inclination
   * of equatorial plane relative to plane of the sky, in radians
   */
  void setInclination(double);

  /// Set inclination relative to line-of-sight
  /**
   * Inclination of z-axis relative to line-of-sight, or inclination
   * of equatorial plane relative to plane of the sky, in specified unit.
   */
  void setInclination(double, const std::string &unit);

  void setPALN(double);
           ///< Set position angle of the line of nodes
  void setPALN(double, const std::string &unit);
           ///< Set position angle of the line of nodes
  void setArgument(double);
           ///< Set angle beetwen line of nodes and X axis of object
  void setArgument(double, const std::string &unit);
           ///< Set angle beetwen line of nodes and X axis of object
  void setSpectrometer(SmartPointer<Spectrometer::Generic> spectro);
           ///< Set Screen::spectro_
  SmartPointer<Spectrometer::Generic> getSpectrometer() const ;
           ///< Get Screen::spectro_

  /**
   * \brief Set freq_obs_
   * \param fo double: observing frequency in Hz
   */
  void setFreqObs(double fo);


  /**
   * \brief Set freq_obs_
   * \param fo double: observing frequency (or wavelength) in "unit"
   * \param unit string: unit in which fo is expressed, convertable to
   * Herz or meters or energy.
   */
  void setFreqObs(double fo, const std::string &unit);

  /**
   * \brief Get freq_obs_.
   */
  double getFreqObs() const ;

  /**
   * \brief Get freq_obs_.
   * \param unit string: unit in which freq_obs_ should be returned is
   * expressed, convertable to Herz or meters or energy.
   */
  double getFreqObs(const std::string &unit) const;

  /// Alternative way to set projection
  /**
   * Beware : paln can not be set this way, setting later other
   * parameters change the observer's coordinates. For observationnal
   * ray-tracing purposes, prefer setProjection().
   *
   * \param[in] pos position of observer in Screen's coordinate
   * system. Content is copied.
   */
  void setObserverPos(const double pos[4]);

  void setFourVel(const double coord[4]);
  ///< Sets the observer's 4-velocity
  void setScreen1(const double coord[4]);
  ///< Sets the screen vector e1
  void setScreen2(const double coord[4]);
  ///< Sets the screen vector e2
  void setScreen3(const double coord[4]);
  ///< Sets the screen vector e3 (normal)

  // Accessors
  // ---------

  /// Get coordinate kind
  /**
   * From Screen::gg_.
   */
  int getCoordKind() const;

  /// Get distance from observer
  /**
   * In meters.
   */
  double getDistance() const;

  /// Get distance from observer
  /**
   * In specified unit.
   */
  double getDistance(const std::string&) const;	 ///< Get distance from observer

  /// Get maximum ray-tracing distance
  /**
   * In geometrical units.
   */
  double getDmax() const;

  /// Get inclination relative to line-of-sight
  /**
   * Inclination of z-axis relative to line-of-sight, or inclination
   * of equatorial plane relative to plane of the sky, in radians.
   */
  double getInclination() const;

  /// Get inclination relative to line-of-sight
  /**
   * Inclination of z-axis relative to line-of-sight, or inclination
   * of equatorial plane relative to plane of the sky, in specified unit.
   */
  double getInclination(const std::string&) const;

  double getPALN() const;	 ///< Get position angle of the line of nodes
  double getPALN(const std::string&) const;	 ///< Get position angle of the line of nodes
  double getArgument() const;	 ///< Get angle between line of nodes and X axis of object
  double getArgument(const std::string&) const;	 ///< Get angle between line of nodes and X axis of object

  SmartPointer<Metric::Generic> getMetric() const; ///< Get Screen::gg_
  void setMetric(SmartPointer<Metric::Generic> gg); ///< Set Screen::gg_

  /// Get observing date in seconds
  double getTime();

  /// Get observing date in seconds
  double getTime(const std::string &);

  /// Set observing date in specified unit
  void setTime(double, const std::string &);

  /// Set observing date in seconds
  void setTime(double);

  /// Get Screen::fov_ in radians
  double getFieldOfView();

  /// Get Screen::fov_ in specified unit
  double getFieldOfView(std::string unit);

  /// Set Screen::fov_ in radians
  void setFieldOfView(double);

  /// Set Screen::fov_ in specified unit
  void setFieldOfView(double, const std::string &unit);

  /// Set direction of the line-of-view
  void setAlpha0(double);
  /// Set direction of the line-of-view
  void setDelta0(double);

  /// Set Screen::anglekind_
  void setAnglekind(int);

  /// Get Screen::npix_
  size_t getResolution();
  /// Set Screen::npix_
  void setResolution(size_t);

  /// 4-Position of the observer relative to the metric
  /**
   * A Screen is positioned relative to the observer with four elements:
   * Screen::distance, Screen::inclination, Screen::paln and
   * Screen::argument.
   *
   * This function returns the position of the observer relative to
   * the metric system in Screen::gg_, using these parameters. The
   * output parameter is coord.
   *
   * \param[out] coord position of the observer. Must be preallocated.
   */
  void getObserverPos(double coord[]) const;

  /// Get copy of Screen::fourvel_
  /**
   * \param[out] fourvel preallocated 4-element array
   */
  void getFourVel(double fourvel[]) const;

  /// Get copy of Screen::screen1_
  /**
   * \param[out] output preallocated 4-element array
   */
  void getScreen1(double output[]) const;

  /// Get copy of Screen::screen2_
  /**
   * \param[out] output preallocated 4-element array
   */
  void getScreen2(double output[]) const;

  /// Get copy of Screen::screen3_
  /**
   * \param[out] output preallocated 4-element array
   */
  void getScreen3(double output[]) const;

  /// Get 8-coordinate of Photon hitting screen from a given direction
  /**
   * Similar to Screen::getObserverPos() but will return in addition
   * the 4-velocity of a photon corresponding to the sky direction
   * given by x and y.
   * \param[in] x    RA (d_alpha*cos(delta)) offset in radians;
   * \param[in] y    Dec offset (d_delta) in radians; 
   * \param[out] coord position-velocity of the observer Photon. Preallocated.
   * 
   */
  void getRayCoord(double x, double y, double coord[]) const;

  /// Get 8-coordinate of Photon hitting screen pixel
  /**
   * Similar to Screen::getObserverPos() but will return in addition
   * the 4-velocity of a photon corresponding to the sky direction
   * given by x and y.
   * \param[in] i, j pixel coordinates   
   * \param[out] coord position-velocity of the Photon. Preallocated.
   * 
   */
  void getRayCoord(const size_t i, const size_t j, double coord[]) const;
  
  void coordToSky(const double pos[4], double skypos[3]) const;
  ///< Convert 4-position to 3-sky position

  void coordToXYZ(const double pos[4], double xyz[3]) const;
  ///< Convert 4-position to 3-cartesian coordinates

  void computeBaseVectors() ;
  ///< Compute base vectors according to projection parameters

  /// Display
  //  friend std::ostream& operator<<(std::ostream& , const Screen& ) ;
  std::ostream& print(std::ostream&) const ; ///< Debug helper
  std::ostream& printBaseVectors(std::ostream&) const ; ///< Debug helper

  // UDUNITS
# ifdef HAVE_UDUNITS
  /// Map "pix" and "pixel" to angular pixel width in unit system
  /**
   * "pix" or "pixel" can then be used in units.
   *
   * There is only one unit system in Gyoto: "pix" can therefore be
   * registered only for one Screen at a time. See Gyoto::Units.
   * 
   * The unit must later be unmapped with unmapPixUnit().
   */
  void mapPixUnit();

  /// Unmap "pix" and "pixel" from unit system
  /**
   * See also mapPixUnit().
   */
  void unmapPixUnit();
# endif


#ifdef GYOTO_USE_XERCES
 public:
    void fillElement(FactoryMessenger *fmp); ///< called from Factory
    /// Instanciate a Screen from XML entity 
    static   SmartPointer<Screen> Subcontractor(FactoryMessenger* fmp);
#endif


};

#endif
