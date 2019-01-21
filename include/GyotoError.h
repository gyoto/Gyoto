/**
 * \file GyotoError.h
 * \brief Error handling
 *
 *  Gyoto dlopens its plug-ins. The throw/catch C++ mechanism cannot
 *  pass the dlopen boundary. The Gyoto::Error mechanism alleviates
 *  this C++ language limitation.
 *
 *  Every Gyoto method (either in the main Gyoto library or in a Gyoto
 *  plug-in) should check for possible error conditions and throw
 *  adequate Gyoto::Error exceptions through the
 *  Gyoto::Error::throw() function. For instance:
 *  \code
 *  if (error_condition) Gyoto::Error::throw("Useful error message");
 *  \endcode
 *
 *  If the main code has set Gyoto::Error::handler_t error handler
 *  using Gyoto::Error::setHandler(), these errors will then be passed
 *  to it. Else, the Error is C++-thrown at the main Gyoto library
 *  level, above the dlopen boundary.
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
    Copyright 2011, 2013 Thibaut Paumard

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
#include "GyotoDefs.h"

namespace Gyoto {
  class Error;
}

/**
 *  \class Gyoto::Error
 *  \brief Class for thowing exceptions.
 *
 *  Gyoto dlopens its plug-ins. The throw/catch C++ mechanism cannot
 *  pass the dlopen boundary. The Gyoto::Error mechanism alleviates
 *  this C++ language limitation.
 *
 *  Every Gyoto method (either in the main Gyoto library or in a Gyoto
 *  plug-in) should check for possible error conditions and throw
 *  adequate Gyoto::Error exceptions through the GYOTO_ERROR macro
 *  (which calls the Gyoto::throwError() function). For instance:
 *  \code
 *  if (error_condition) GYOTO_ERROR("Useful error message");
 *  \endcode
 *
 *  If the main code has set Gyoto::Error::handler_t error handler
 *  using Gyoto::Error::setHandler(), these errors will then be passed
 *  to it. Else, the Error is C++-thrown at the main Gyoto library
 *  level, above the dlopen boundary.
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
 */
class Gyoto::Error
{
 private:
  /// Error message.
  const std::string message;

  /// Error code.
  /**
   * Default value is EXIT_FAILURE from cstdlib. Currently not used in
   * practice.
   */
  const int errcode;

 public:

  /// Constructor with an error message.
  Error( const std::string m );

  // Copy constructor
  Error( const Gyoto::Error &o);

  /// Print-out error message on standard error.
  void Report() const ;

  /// Retrieve error code.
  /**
   * See also operator const char * () const and get_message().
   * \return Error code
   */
  int getErrcode() const ;

  /// Cast to const char *.
  /**
   * Retrieve error message as a C string. See also get_message() and
   * gerErrcode().
   */
  operator const char * () const;

  /// Retrieve error message for custom handling of the exception.
  /**
   * See also operator const char * () const and getErrCode().
   * \return char* message : pointer to the error message
   */
  std::string get_message() const ;

  /// Error handler type.
  /**
   * Instead of catching Gyoto errors directly (for instance if gyoto
   * itself is dlopened), you can set a Handler_t error handler using
   * setHandler().
   *
   * A very simple handler could be:
   * \code
   * void applicationErrorHandler(const Gyoto::Error e) {
   *   e.Report();
   *   exit ( e.getErrCode() );
   * }
   * \endcode
   */
  typedef void Handler_t (const Error);

  /// Set application error handler.
  /**
   * Instead of catching Gyoto errors directly (for instance if gyoto
   * itself is dlopened), you can set an Error::Handler_t error
   * handler using setHandler().
   *
   * \code
   * void applicationErrorHandler(const Gyoto::Error e) {
   *   e.Report();
   *   exit ( e.getErrCode() );
   * }
   * int main() {
   *   Gyoto::Error::setHandler(&applicationErrorHandler);
   * }
   * \endcode
   * \param phandler Function pointer to the handler.
   */
  static void setHandler( Gyoto::Error::Handler_t* phandler);

};

namespace Gyoto {
  /// Throw a Gyoto::Error
  /**
   * Most code should use the GYOTO_ERROR macro instead
   */
  void throwError( std::string );
}


/// Throw a Gyoto::Error nicely
/**
 * Throw an Error, prepending current function name. Calls
 * Gyoto::throwError(std::string).
 */
#define GYOTO_ERROR(msg) Gyoto::throwError(std::string(__FILE__ ":" GYOTO_STRINGIFY(__LINE__) " in ")+ __PRETTY_FUNCTION__ + ": " + msg)

#endif
