/**
 *  \file GyotoKerr2PN.h
 *  \brief Kerr spacetime at second post-Newtonian order
 *         in Cartesian harmonic coordinates (t,x,y,z).
 */

/*
    Copyright 2026 Karim Abd El Dayem, Frederic Vincent & Thibaut Paumard

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

#ifndef __GyotoKerr2PN_h
#define __GyotoKerr2PN_h

#include <GyotoMetric.h>

namespace Gyoto {
  namespace Metric {
    class Kerr2PN;
  };
};

/**
  \class Gyoto::Metric::Kerr2PN
  \brief Kerr metric in harmonic Cartesian coordinates, truncated at 1PN / 2PN.

  Both gmunu() and christoffel() are approximated at the second
  post-Newtonian order.

  Conventions (see the document "Christoffel Symbols in the Kerr Metric,
  Harmonic Coordinates, Second-Order Post-Newtonian Approximation",
  K. Abd El Dayem):

   - Geometrized units G = c = M = 1.
   - Cartesian harmonic coordinates x^mu = (t, x, y, z), pos[0..3].
   - Spin axis along z: e = (0,0,1), J = (0,0,spin_).
   - n = x/r, n.e = z/r.

  PERFORMANCE NOTES
  -----------------
  christoffel() is the hot path of the integrator (called several times
  per step), gmunu(g,pos) is called by every norm / scalar-product check.
  Both are therefore written as straight-line code:

   * one sqrt and one division per call (everything else is multiplication
     by powers of ir = 1/r),
   * no Levi-Civita symbol, no function calls, no branches inside loops:
     with J = spin*z_hat the frame-dragging terms have a closed form
     (derived below) and only the non-zero components are written,
   * the PN order and spin tests are hoisted out of all loops.

  Closed form of the gravitomagnetic part, h_i := g_{0i} = 2 (x X J)_i / r^3
  with J = (0,0,chi):
      h = ( 2 chi y / r^3 , -2 chi x / r^3 , 0 )
      d_x h_x = -6 chi x y / r^5          d_y h_x =  2chi/r^3 - 6 chi y^2/r^5
      d_z h_x = -6 chi y z / r^5          d_x h_y = -2chi/r^3 + 6 chi x^2/r^5
      d_y h_y =  6 chi x y / r^5          d_z h_y =  6 chi x z / r^5
      d_m h_z = 0
  and
      Gamma^0_ij = -1/2 (d_i h_j + d_j h_i)
      Gamma^i_0j =  1/2 (d_j h_i - d_i h_j)      (antisymmetric in i,j)
  ---------------------------------------------------------------------
*/
class Gyoto::Metric::Kerr2PN
: public Gyoto::Metric::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Metric::Kerr2PN>;

protected:
  double spin_; ///< Geometrized units spin parameter
  double Q_; ///< Geometrized units mass quadrupole parameter
  bool quadrupole_is_kerr_; ///< Boolean: True iff the quadrupole is tied to its Kerr value -spin^2 in geometrized units
  unsigned int pnorder_; ///< Post-Newtonian truncation order (1 or 2)

public:
  GYOTO_OBJECT;
  Kerr2PN();
  Kerr2PN(const Kerr2PN & orig);
  virtual ~Kerr2PN();
  Kerr2PN * clone() const override;

  void spin(const double charge); ///< Sets spin
  double spin() const ; ///< Returns spin
  void quadrupole(const double charge); ///< Sets quadrupole
  double quadrupole() const ; ///< Returns quadrupole
  void kerrQuadrupole(bool t) ; ///< Sets boolean quadrupole_is_kerr_
  bool kerrQuadrupole() const ; ///< Returns boolean quadrupole_is_kerr_
  void pnOrder(double n); ///< Sets PN truncation order (1 or 2)
  double pnOrder() const; ///< Returns PN truncation order

  using Generic::gmunu;
  double gmunu(double const x[4], int mu, int nu) const override;
  void gmunu(double ARGOUT_ARRAY2[4][4], double const IN_ARRAY1[4]) const override;
  using Generic::christoffel;
  /**
   * \brief Chistoffel symbol
   *
   * Value of Christoffel symbol
   * &Gamma;<SUP>&alpha;</SUP><SUB>&mu;&nu;</SUB> at point
   * (x<SUB>1</SUB>, x<SUB>2</SUB>, x<SUB>3</SUB>).
   *
   * The expression is truncated at the 2PN order
   *
   * \return 1 on error, 0 otherwise
   */
  int christoffel(double dst[4][4][4], double const pos[4]) const override;
  int isStopCondition(double const coord[8]) const override;
  //void circularVelocity(double const * coor, double* vel, double dir) const override;
  //double getPotential(double const pos[4], double l_cst) const override;
  //double getSpecificAngularMomentum(double rr) const override;
};
#endif
