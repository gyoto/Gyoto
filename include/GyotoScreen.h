/**
 * \file GyotoScreen.h
 * \brief Description of the observer screen
 * 
 */

/*
    Copyright 2011-2019 Thibaut Paumard, Frederic Vincent

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
#if defined HAVE_BOOST_ARRAY_HPP
# include <boost/array.hpp>
# define GYOTO_ARRAY boost::array
# if defined HAVE_MPI
#  include <boost/version.hpp>
#  if BOOST_VERSION >= 106400 
#   include <boost/serialization/boost_array.hpp>
#   include <boost/serialization/array_wrapper.hpp>
#  endif
# endif
#else
template <typename T, size_t sz> class GYOTO_ARRAY {
 private:
  T buf[sz];
 public:
  T& operator[](size_t c) { return buf[c] ; }
};
#endif

namespace Gyoto {
  class Screen;
}

#include <GyotoDefs.h>
#include <GyotoUtils.h>
#include <GyotoSmartPointer.h>
#include <GyotoObject.h>
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
 * A mask may be used to limit ray-tracing to only some portions of
 * the field. The Scenery checks whether a mask is to be used using
 * Screen::operator()(size_t i, size_t j). The mask can be loaded from
 * a FITS file as a square image of doubles:
 * \code
 *    <Mask>maskfile.fits</Mask>
 * \endcode
 * The mask needs to be have the same size as the Screen itself, so
 * loading a mask also sets the resolution, and changing the
 * resolution after setting a mask also removes the mask. The content
 * of the Mask entity is parsed by Factory::fullPath(), so it can be
 * an absolute path, a path relative to where the XML file is stored,
 * or relative to the current working directory if prefixed with
 * "`pwd`/".
 *
 */
class Gyoto::Screen
: public Gyoto::SmartPointee,
  public Gyoto::Object
{
  friend class Gyoto::SmartPointer<Gyoto::Screen>;

 private:
  double tobs_; ///< Observing date in s
  double fov_;  ///< Field-of-view in rad
  double azimuthal_fov_; ///< Azimuthal field-of-view for Spherical Angles images. Maximal extent of image in the azimuthal b-angle direction.
  //  double tmin_;
  size_t npix_; ///< Resolution in pixels

  /**
   * \brief Mask with 0 where the ray-tracing should not be performed
   */
  double * mask_;

  /**
   * \brief Last read or written FITS file
   *
   * Used when saving to XML: if the mask was saved or loaded from
   * FITS file, output this file name in the XML.
   */
  std::string mask_filename_;

  double distance_; ///< Distance to the observer in m
  double dmax_; ///< Maximum distance from which the photons are launched (geometrical units) 

  enum anglekind_e { equatorial_angles=0, rectilinear=1, spherical_angles=2};
  typedef int anglekind_t;

  anglekind_t anglekind_; ///< Screen angles kind (0: equatorial, 1: spherical)
  
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

  double dangle1_; ///< Increment to first position angle of Screen; can be typically alpha if in Equatorial Angles, or a if in Spherical Angles
  double dangle2_; ///< Increment to second position angle of Screen; can be typically delta if in Equatorial Angles, or b if in Spherical Angles
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

  /**
   * \brief What kind of observer are we considering? (At infinity, ZAMO...)
   *
   */
  obskind_t observerkind_;

 public:
  GYOTO_OBJECT;
  GYOTO_OBJECT_THREAD_SAFETY;

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
  void distance(double dist);

  /// Set ray-tracing maximum distance
  /**
   * \param dist Distance in geometrical units.
   */
  void dMax(double dist);

  /// Set distance from observer
  /**
   * \param dist the distance expressed in the specified unit;
   * \param unit convertible to meters
   */
  void distance(double dist, const std::string &unit);

  /// Set inclination relative to line-of-sight
  /**
   * Inclination of z-axis relative to line-of-sight, or inclination
   * of equatorial plane relative to plane of the sky, in radians
   */
  void inclination(double);

  /// Set inclination relative to line-of-sight
  /**
   * Inclination of z-axis relative to line-of-sight, or inclination
   * of equatorial plane relative to plane of the sky, in specified unit.
   */
  void inclination(double, const std::string &unit);

  void PALN(double);
           ///< Set position angle of the line of nodes
  void PALN(double, const std::string &unit);
           ///< Set position angle of the line of nodes
  void argument(double);
           ///< Set angle beetwen line of nodes and X axis of object
  void argument(double, const std::string &unit);
           ///< Set angle beetwen line of nodes and X axis of object
  void spectrometer(SmartPointer<Spectrometer::Generic> spectro);
           ///< Set Screen::spectro_
  SmartPointer<Spectrometer::Generic> spectrometer() const ;
           ///< Get Screen::spectro_

  /**
   * \brief Set freq_obs_
   * \param fo double: observing frequency in Hz
   */
  void freqObs(double fo);


  /**
   * \brief Set freq_obs_
   * \param fo double: observing frequency (or wavelength) in "unit"
   * \param unit string: unit in which fo is expressed, convertible to
   * Herz or meters or energy.
   */
  void freqObs(double fo, const std::string &unit);

  /**
   * \brief Get freq_obs_.
   */
  double freqObs() const ;

  /**
   * \brief Get freq_obs_.
   * \param unit string: unit in which freq_obs_ should be returned is
   * expressed, convertible to Herz or meters or energy.
   */
  double freqObs(const std::string &unit) const;

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
  void observerKind(const std::string &kind);
  std::string observerKind() const;
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
  int coordKind() const;

  /// Get distance from observer
  /**
   * In meters.
   */
  double distance() const;

  /// Get distance from observer
  /**
   * In specified unit.
   */
  double distance(const std::string&) const;	 ///< Get distance from observer

  /// Get maximum ray-tracing distance
  /**
   * In geometrical units.
   */
  double dMax() const;

  /// Get inclination relative to line-of-sight
  /**
   * Inclination of z-axis relative to line-of-sight, or inclination
   * of equatorial plane relative to plane of the sky, in radians.
   */
  double inclination() const;

  /// Get inclination relative to line-of-sight
  /**
   * Inclination of z-axis relative to line-of-sight, or inclination
   * of equatorial plane relative to plane of the sky, in specified unit.
   */
  double inclination(const std::string&) const;

  double PALN() const;	 ///< Get position angle of the line of nodes
  double PALN(const std::string&) const;	 ///< Get position angle of the line of nodes
  double argument() const;	 ///< Get angle between line of nodes and X axis of object
  double argument(const std::string&) const;	 ///< Get angle between line of nodes and X axis of object

  SmartPointer<Metric::Generic> metric() const; ///< Get Screen::gg_
  void metric(SmartPointer<Metric::Generic> gg); ///< Set Screen::gg_

  /// Get observing date in seconds
  double time() const;

  /// Get observing date in seconds
  double time(const std::string &) const;

  /// Set observing date in specified unit
  void time(double, const std::string &);

  /// Set observing date in seconds
  void time(double);

  /// Get Screen::fov_ in radians
  double fieldOfView() const;

  /// Get Screen::fov_ in specified unit
  double fieldOfView(std::string const &unit) const;

  /// Set Screen::fov_ in radians
  void fieldOfView(double);

  /// Set Screen::fov_ in specified unit
  void fieldOfView(double, const std::string &unit);

  /// Get Screen::azimuthal_fov_
  double azimuthalFieldOfView() const;

  /// Set Screen::azimuthal_fov_
  void azimuthalFieldOfView(double ff);

  /// Set increment to first position angle
  void dangle1(double);
  /// Set increment to first position angle in specified unit
  void dangle1(double, const std::string &unit);
  /// Get increment to first position angle
  double dangle1() const;
  /// Get increment to first position angle in specified unit
  double dangle1(std::string const &unit)const;
  /// Set increment to second position angle
  void dangle2(double);
  /// Set increment to second position angle in specified unit
  void dangle2(double, const std::string &unit);
  /// Get increment to second position angle
  double dangle2() const;
  /// Get increment to second position angle in specified unit
  double dangle2(std::string const &unit)const;

  /// Set Screen::anglekind_
  void anglekind(int);
  void anglekind(std::string const&);
  std::string anglekind() const;

  /// Get Screen::npix_
  size_t resolution() const;
  /// Set Screen::npix_
  void resolution(size_t);

  /// Set mask_ from array
  /**
   * mm will be copied. mm must be a square resolution x resolution
   * array. If mm==NULL, just deallocate mask_.
   */
  void mask(double const * const mm, size_t resolution=0);

  /// Retrieve const pointer to mask_
  double const * mask() const ;
  void maskFile(std::string const &fname);
  std::string maskFile() const;
# ifdef GYOTO_USE_CFITSIO

  /// Read mask_ from FITS file
  void fitsReadMask(std::string const &fname);

  /// Save mask_ from FITS file
  void fitsWriteMask(std::string const &fname);
# endif

  /// Whether this pixel should be ray-traced
  /**
   * If mask_ is not set, always true. Else, true for non-zero cells
   * in mask_.
   */
  bool operator()(size_t, size_t);


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
  void getObserverPos(double dest[4]) const;

  /// Get copy of Screen::fourvel_
  /**
   * \param[out] fourvel preallocated 4-element array
   */
  void getFourVel(double dest[4]) const;

  void fourVel(std::vector<double> const &);
  std::vector<double> fourVel() const;
  void screenVector1(std::vector<double> const &);
  std::vector<double> screenVector1() const;
  void screenVector2(std::vector<double> const &);
  std::vector<double> screenVector2() const;
  void screenVector3(std::vector<double> const &);
  std::vector<double> screenVector3() const;

  /// Get copy of Screen::screen1_
  /**
   * \param[out] dest preallocated 4-element array
   */
  void getScreen1(double dest[4]) const;

  /// Get copy of Screen::screen2_
  /**
   * \param[out] dest preallocated 4-element array
   */
  void getScreen2(double dest[4]) const;

  /// Get copy of Screen::screen3_
  /**
   * \param[out] dest preallocated 4-element array
   */
  void getScreen3(double dest[4]) const;

  /// Get 8-coordinate of Photon hitting screen from a given direction and polarization basis if needed
  /**
   * Similar to Screen::getObserverPos() but will return in addition
   * the 4-velocity of a photon corresponding to the sky direction
   * given by x and y.
   * \param[in] x    RA (d_alpha*cos(delta)) offset in radians;
   * \param[in] y    Dec offset (d_delta) in radians; 
   * \param[out] dest position-velocity of the observer Photon. Preallocated.
   * \param[in] compute_polar_basis True if polarization basis Ephi,Etheta is needed
   * \param[out] Ephi first polarisation direction. Preallocated. Default: NULL.
   * \param[out] Etheta second polarisation direction. Preallocated. Default: NULL.
   * 
   */
  void getRayTriad(double x, double y,
		   double dest[8],
		   bool compute_polar_basis=false,
		   double Ephi[4]=NULL, double Etheta[4]=NULL) const;

  /// Get 8-coordinate of Photon hitting screen pixel and polarization basis if needed
  /**
   * Similar to Screen::getObserverPos() but will return in addition
   * the 4-velocity of a photon corresponding to the sky direction
   * given by x and y.
   * \param[in] i, j pixel coordinates  
   * \param[out] dest position-velocity of the Photon. Preallocated.
   * \param[in] compute_polar_basis True if polarization basis Ephi,Etheta is needed
   * \param[out] Ephi first polarisation direction. Preallocated. Default: NULL.
   * \param[out] Etheta second polarisation direction. Preallocated. Default: NULL.
   * 
   */
  void getRayTriad(const size_t i, const size_t j,
		   double dest[8],
		   bool compute_polar_basis=false,
		   double Ephi[4]=NULL, double Etheta[4]=NULL) const;
  
  /** \brief Convert metric 4-position to sky 3-position
   *
   * \param[in] pos 4-position in metric coordinates.
   * \param[in] dest 3-position in plane of the sky: Cartesian East, North, front.
   * \param[in] geometrical: if true, #dest will be in geometrical units instead of meters.
   */
  void coordToSky(const double pos[4], double dest[3], bool geometrical=false) const;

  /** \brief Convert sky 3-position to metric 4-position
   *
   * \param[in] sky 3-position in plane of the sky.
   * \param[in] dest 4-position in metric coordinates (dest[0] is not modified).
   * \param[in] geometrical: set to true if #sky is in geometrical units instead of meters.
   */
  void skyToCoord(const double sky[3], double dest[4], bool geometrical=false) const;

  void coordToXYZ(const double pos[4], double dest[3]) const;
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
  void fillProperty(Gyoto::FactoryMessenger *fmp, Property const &p) const;

    /// Instanciate a Screen from XML entity 
    static   SmartPointer<Screen> Subcontractor(FactoryMessenger* fmp);
#endif

    /// Enum to specify whether a coordinate set (Coord1dSet or Coord2dSet) holds pixel values or angles
    enum CoordType_e {angle, pixel};

    /// Set of 1-d coordinates: indices or angles
    /**
     * Acts like a container (array-like) of either size_t (pixel
     * coordinate) or double (angle) values. This container can be
     * iterated-through using the operator++(), derefenced using the
     * operator*() (if containing pixel coordinates) or angle() (in
     * containing angles).
     */
    class Coord1dSet {
    public:
      /// Whether this specifier represents angles or pixels
      const CoordType_e kind;
    public:
      /// Set kind during initialization
      Coord1dSet(CoordType_e k);
      /// Virtual destructor
      virtual ~Coord1dSet();
      /// Reset specifier to point to the first value
      virtual void begin() =0;
      /// True if pointing to something, false if end has been reached.
      virtual bool valid() =0;
      /// Number of values in this container
      virtual size_t size()=0;
      /// Get size_t value currently pointed to
      virtual size_t operator*() const ;
      /// Get double value currently pointed to
      virtual double angle() const ;
      /// Increment iterator (point to next value)
      virtual Coord1dSet& operator++()=0;
      /// Get index of value currently pointed to
      /**
       * Starts at 0 and is implemented each time operator++ is
       * called. Depending on the implementation, this may be a real
       * index or computed on demand.
       */
      virtual size_t index() const=0;
    };

    /// Class to specify a set of points on the Screen
    /**
     * Container (array-like) holding several 2D points. Can be a 2D
     * grid of pixel coordinates or a vector of floating-point (alpha,
     * delta) pairs, for instance.
     */
    class Coord2dSet {
    public:
      /// Whether this set holds pixels or angle specifications
      const CoordType_e kind;
      /// Set kind at initialisation
      Coord2dSet(CoordType_e k);
      /// Virtual destructor
      virtual ~Coord2dSet();
      /// Increment pointer
      virtual Coord2dSet& operator++()    =0;
      /// Get pixel coordinates
      virtual GYOTO_ARRAY<size_t, 2> operator*  () const;
      /// Get angle coordinates
      virtual GYOTO_ARRAY<double, 2> angles() const ;
      /// Reset pointer
      virtual void begin() =0;
      /// Whether the end has not been passed
      virtual bool valid() =0;
      /// Number of positions contained
      virtual size_t size()=0;
    };

    /// Class containing 2D-points organized in a grid
    class Grid: public Coord2dSet {
    protected:
    protected:
      /// If non-NULL, cout j each tims it is incremented.
      char * prefix_;
      Coord1dSet &iset_;
      Coord1dSet &jset_;
    public:
      Grid(Coord1dSet &iset, Coord1dSet &jset, const char * const p=NULL);
      virtual ~Grid();
      virtual Coord2dSet& operator++();
      virtual GYOTO_ARRAY<size_t, 2> operator*  () const;
      virtual void begin();
      virtual bool valid();
      virtual size_t size();
    };

    /// Class containing arbitrary 2D-points 
    /**
     * ispec_ and jspec_ must be the same size.
     */
    class Bucket : public Coord2dSet {
    protected:
      Coord1dSet &alpha_;
      Coord1dSet &delta_;
    public:
      Bucket(Coord1dSet &iset, Coord1dSet &jset);
      virtual Coord2dSet& operator++();
      virtual GYOTO_ARRAY<double, 2> angles() const;
      virtual GYOTO_ARRAY<size_t, 2> operator*() const;
      virtual void begin();
      virtual bool valid();
      virtual size_t size();
    };

    /// A dummy, empty 2D set.
    class Empty: public Coord2dSet {
    public:
      Empty();
      virtual Coord2dSet& operator++();
      virtual void begin();
      virtual bool valid();
      virtual size_t size();
    };

    /// 1D coordinated specifier for a range
    class Range : public Coord1dSet {
    protected:
      const size_t mi_, ma_, d_, sz_;
      size_t cur_;
    public:
      /// Specify min, max and step of this range.
      Range(size_t mi, size_t ma, size_t d);
      void begin();
      bool valid();
      size_t size();
      Coord1dSet& operator++();
      size_t operator*() const ;
      virtual size_t index() const ;
    };

    /// 1D specifier for an arbitrary pixel coordinate set.
    class Indices : public Coord1dSet {
    protected:
      size_t * indices_;
      size_t const sz_;
      size_t i_;
    public:
      Indices (size_t const*const buf, size_t sz);
      ~Indices();
      void begin();
      bool valid();
      size_t size();
      Coord1dSet& operator++();
      size_t operator*() const ;
      virtual size_t index() const ;
    };

    /// 1D specifier for an arbitrary angle coordinate set.
    class Angles : public Coord1dSet {
    protected:
      double * buf_;
      size_t const sz_;
      size_t i_;
    public:
      Angles (double const*const buf, size_t sz);
      ~Angles();
      void begin();
      bool valid();
      size_t size();
      Coord1dSet& operator++();
      double angle() const ;
      virtual size_t index() const ;
    };

    /// 1D specifier for an angle that is repeated.
    class RepeatAngle : public Coord1dSet {
    protected:
      double const val_;
      size_t const sz_;
      size_t i_;
    public:
      RepeatAngle (double val, size_t sz);
      void begin();
      bool valid();
      size_t size();
      Coord1dSet& operator++();
      double angle() const ;
      virtual size_t index() const ;
    };
};

#endif
