/**
 * \file GyotoPython.h
 * \brief Extending Gyoto using Python
 *
 */

/*
    Copyright 20015 Thibaut Paumard

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

#ifndef __GyotoPython_H_ 
#define __GyotoPython_H_ 
#include <GyotoSpectrum.h>
#include <Python.h>

namespace Gyoto {
  namespace Spectrum {
    class Python;
  }
}


/**
 * \class Gyoto::Spectrum::Python
 * \brief Spectrum using Python
 *
 *  Light emitted by e.g. a Star.
 *
 *  XML stanza:
 *  \code
 *    <Spectrum kind="Python">
 *      <Parameters> 0. 1. 2. ... </Parameters>
 *      <Function> 
 *      </Function>
 *    </Spectrum>
 *  \endcode
 */
class Gyoto::Spectrum::Python : public Gyoto::Spectrum::Generic {
  friend class Gyoto::SmartPointer<Gyoto::Spectrum::Python>;
 protected:
  std::string module_;
  std::string class_;
  PyObject * pModule_;
  PyObject * pClass_;
  PyObject * pInstance_;
  std::vector<double> parameters_;

 public:
  GYOTO_OBJECT;

  Python();

  Python(const Python&);

  virtual Python * clone() const;

  ~Python();
  
  std::string module() const ;
  void module(const std::string&);

  std::string klass() const ;
  void klass(const std::string&);

  std::vector<double> parameters() const;
  void parameters(const std::vector<double>&);

  using Gyoto::Spectrum::Generic::operator();
  virtual double operator()(double nu) const;

};

#endif
