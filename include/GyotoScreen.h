/**
 * \file GyotoScreen.h
 * \brief Description of the observer screen
 * 
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
 *   - a description of the Spectrometer.
 *
 * The Spectrometer defines for which frequencies spectra are computed
 * (when the Quantity Spectrum is requested in the Scenery). The
 * frequency axis can be linear in frequency or in wavelength, or
 * logarithmic in frequency or in wavelength. The spectrometr is thus
 * defines by:
 *   - a kind: one of "none", "freqlog", "freq", "wavelog" or "wave";
 *   - a number of spectral channels;
 *   - the boundaries of the spectral band.
 *
 * Each spectral channel will be of equal width in the space specified
 * by the kind of the Spectrometer. The boundaries are specified in Hz
 * for a "freq" Spectrometer, in Log(Hz) for freqlog, in meters for wave
 * and Log(meter) for wavelog. The boundaries correspond to the extremum
 * boundaries of the outermost spectral channels. When computing the
 * Quantity "Spectrum", the frequency considered is the midpoint of
 * each spectral channel.
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
 *    <Time>       1000. </Time>
 *    <FieldOfView>   0.3141592653589793 </FieldOfView>
 *    <Resolution>  128 </Resolution>
 *    <Distance>      1e30 </Distance>
 *    <PALN>          3.14159 </PALN>
 *    <Inclination>   2.0944 </Inclination>
 *    <Argument>     -2.0944 </Argument>
 *    <Spectrometer kind="freqlog" nsamples="10"> 17. 23. </Spectrometer> 
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
 *  </Screen>
 * \endcode
 *
 *
 * Units can be specified using the unit attribute in the XML file,
 * for instance:
 * 
 * \code
 *   <Distance unit="kpc"> 8 </Distance>
 * \endcode.
 *
 * The possible units are (with [] noting the default):
 *  - distance: [m], geometrical, cm, km, AU, ly, pc, kpc, Mpc;
 *  - PALN, inclination, argument: [rad], deg.
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
 * \encode
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
  size_t npix_; ///< resolution in pixels

  double distance_; ///< Distance to the observer in m
  double dmax_; ///< Maximum distance from which the photons are launched (geometrical units) 
  
  /**
   * The angles are position angle of the line of nodes (North of
   * East), inclination (0 = face-on), argument of X axis. We use the
   * z-x-z convention. See http://en.wikipedia.org/wiki/Euler_angles
   */
  double euler_[3]; ///< Euler angles
  double ex_[3]; ///< Sky coordinate of base X vector
  double ey_[3]; ///< Sky coordinate of base Y vector
  double ez_[3]; ///< Sky coordinate of base Z vector

  SmartPointer<Metric> gg_; ///< Metric in which the screen is placed ; necessary for unitLength

  SmartPointer<Spectrometer> spectro_;

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

  /**
   * \param dist the distance in meters.
   */
  void setDistance(double dist);    ///< Set distance from observer
  void setDmax(double dist);    ///< Set ray-tracing maximum distance
  /**
   * \param dist the distance expressed in the specified unit;
   * \param unit one of: [m], geometrical, cm, km, sunradius, AU, ly,
   *             pc, kpc, Mpc.
   */
  void setDistance(double dist, const std::string unit);
           ///< Set distance from observer
  void setInclination(double);
           ///< Set inclination relative to line-of-sight
  void setInclination(double, const std::string &unit);
  void setPALN(double);
           ///< Set position angle of the line of nodes
  void setPALN(double, const std::string &unit);
  void setArgument(double);
           ///< Set angle beetwen line of nodes and X axis of object
  void setArgument(double, const std::string &unit);
  void setSpectrometer(SmartPointer<Spectrometer> spectro);
  SmartPointer<Spectrometer> getSpectrometer() const ;

  /// Alternative way to set projection
  /**
   * Beware : paln can not be set this way, setting later other
   * parameters change the observer's coordinates. For observationnal
   * ray-tracing purposes, prefer setProjection().
   *
   * \param pos[4] position of observer in Screen's coordinate system
   * \param sys    coordinate system used for pos
   */
  void setObserverPos(const double pos[4]);
  ///< Sets the orientation of Screen relative to observer

  // Accessors
  // ---------

  int getCoordKind() const;      ///< Get coordinate kind
  double getDistance() const;	 ///< Get distance from observer
  double getDmax() const;	 ///< Get maximum ray-tracing distance
  double getInclination() const; ///< Get inclination relative to line-of-sight
  double getPALN() const;	 ///< Get position angle of the line of nodes
  double getArgument() const;	 ///< Get angle beetwen line of nodes and X axis of object


  SmartPointer<Metric> getMetric() const; ///< Get gg_;
  void setMetric(SmartPointer<Metric> gg);

  double getTime();
  void setTime(double, const std::string &);
  void setTime(double);
  //  double getMinimumTime();
  //  void setMinimumTime(double);
  double getFieldOfView();
  void setFieldOfView(double);
  void setFieldOfView(double, const std::string &unit);
  size_t getResolution();
  void setResolution(size_t);


  /**
   *A Screen is positioned relative to the observer with four elements:
   * Screen::distance, Screen::inclination, Screen::paln and
   * Screen::argument.
   *
   * This function returns the position of the observer relative to
   * the Screen, using these parameters. The output parameter is
   * coord.
   *
   * \param tobs observing time;
   * \param coord[4] output: position of the observer;
   *
   */
  void getObserverPos(double coord[]) const;
  ///< 4-Position of the observer relative to the Screen

  
  /**
   * Similar to Screen::getObserverPos() but will return in addition
   * the 4-velocity of a photon corresponding to the sky direction
   * given by x and y.
   * \param tobs observing time;
   * \param x    RA (d_alpha*cos(delta)) offset in radians;
   * \param y    Dec offset (d_delta) in radians; 
   * \param coord[8] output: position-velocity of the observer in system sys;
   * 
   */
  void getRayCoord(double x, double y, double coord[]) const;
  void getRayCoord(const size_t i, const size_t j, double coord[]) const;
  
  void coordToSky(const double pos[4], double skypos[3]) const;
  ///< Convert 4-position to 3-sky position

  void coordToXYZ(const double pos[4], double xyz[3]) const;
  ///< Convert 4-position to 3-cartesian coordinates

  void computeBaseVectors() ;
  ///< Compute base vectors according to projection parameters




  /// Display
  friend std::ostream& operator<<(std::ostream& , const Screen& ) ;
  std::ostream& print(std::ostream&) const ;
  std::ostream& printBaseVectors(std::ostream&) const ;

#ifdef GYOTO_USE_XERCES
 public:
    void fillElement(FactoryMessenger *fmp); /// < called from Factory
    static   SmartPointer<Screen> Subcontractor(FactoryMessenger* fmp);
#endif


};

#endif
