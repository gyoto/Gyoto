/**
 * \file GyotoError.h
 * \brief Error handling
 *
 *  Every Gyoto method should check for possible error conditions and
 *  throw adequate Gyoto::Error exceptions. For instance:
 *  \code
 *  if (error_condition) throw Gyoto::Error("Useful error message");
 *  \endcode
 *
 *  The main code can then catch these exceptions and act appropriately,
 *  for instance:
 *  \code
 *  try { gyoto_code ; }
 *  catch (Gyoto::Error err)
 *  {
 *     err.Report();
 *     abort();
 *  }
 * \endcode
 *
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

#ifndef __GyotoError_H_
#define __GyotoError_H_

/**
 * \namespace Gyoto
 * \brief Namespace for the Gyoto library
 */

#include <string>

namespace Gyoto {
  class Error;
}

/**
 *  \class Gyoto::Error
 *  \brief Class for thowing exceptions.
 *
 *  Every Gyoto method should check for possible error conditions and
 *  throw adequate Gyoto::Error exceptions. For instance:
 *  \code
 *  if (error_condition) throw Gyoto::Error("Useful error message");
 *  \endcode
 */
class Gyoto::Error
{
 private:
  //const char* message; /*!< error message */
  const std::string message; /*!< error message */
  const int errcode;
 public:

  /**
   * \brief Constructor with an error message
   * \param m : pointer (char*) to the error message
   */
  //Error( const char* m );
  Error( const std::string m );

  /**
   * \brief Constructor with an error code
   * \param int errcode : error code
   */
  Error( const int errcode );


  /**
   * \brief Constructor with both an error message and an error code
   * \param m : pointer (char*) to the error message
   * \param int errcode : error code
   */
  Error( const char* m , const int errcode );

  /**
   * \brief Print-out error message on standard error
   */
  void Report() const ;  /*!<  */

  /**
   * \brief Retrieve error code
   * \return Error code
   */
  int getErrcode() const ;

  /**
   * \brief Retrieve error message for custom handling of the exception
   * \return char* message : pointer to the error message
   */
  //char const * const get_message() const ;
  std::string get_message() const ;
};

typedef void GyotoErrorHandler_t (const char*);

namespace Gyoto {
  void setErrorHandler( GyotoErrorHandler_t* );
  void throwError( std::string );
}

#endif
