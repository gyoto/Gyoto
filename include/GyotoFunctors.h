/**
 * \file GyotoFunctors.h
 * \brief Classes with an operator() method
 */

/*
    Copyright 2011 Thibaut Paumard

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


#ifndef __GyotoFunctors_H_
#define __GyotoFunctors_H_

namespace Gyoto {
  /**
   * \namespace Gyoto::Functor
   * \brief Classes with an operator() method
   */
  namespace Functor {
    class Double_constDoubleArray;
    class Double_Double_const;
  }
}

/**
 * \brief A functor like double (func) (double const data[])
 */
class Gyoto::Functor::Double_constDoubleArray
{
 public:
  virtual ~Double_constDoubleArray();
  /**
   * \brief The actual function
   */
  virtual double operator()(double const data[]) = 0;
};


/**
 * \brief A functor like double (func) (double) const
 */
class Gyoto::Functor::Double_Double_const
{
 public:
  virtual ~Double_Double_const();
  /**
   * \brief Exit status code of "various" methods (at least secant() !)
   */
  int status;

  /**
   * \brief The actual function
   */
  virtual double operator()(double) const = 0;

  /**
   * \brief Ridder's root-finding method applied on operator()()
   * \param from, to boundaries for root-searching
   * \return the root
   */
  double ridders(double from, double to) const;

  /**
   * \brief Secant root-finding method applied on operator()()
   *
   * Sets status to
   *  -0 in case of convergence
   *  -1 if two distinct inputs evaluated to the same output
   *  -2 if maximum number of iterations (20) reached
   *
   * \param from, to boundaries for root-finding
   * \return the root
   */
  double secant(double from, double to);
};

#endif
