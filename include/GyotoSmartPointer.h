/**
 * \file GyotoSmartPointer.h
 * \brief Reference-counting pointers

  Template class Gyoto::SmartPointer<T> performs reference counting on
  its pointee and makes sure the pointed to object is actually deleted
  when the number of references to it drops to zero. The pointed-to
  object must inherit from Gyoto::SmartPointee for this to work and be
  referenced to only by SmartPointers.

  @code
  class Gyoto::Metric : public Gyoto::SmartPointee {...};
  SmartPointer<Gyoto::Metric::Generic> ObjPtr (new Gyoto::Metric(...));
  @endcode


 */

/*
    Copyright 2011-2014, 2016, 2020 Thibaut Paumard

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


#ifndef __GyotoSmartPointer_H_
#define __GyotoSmartPointer_H_

#include "GyotoUtils.h"
#ifdef HAVE_PTHREAD
#include <pthread.h>
#endif

namespace Gyoto {
  class SmartPointee;
  class FactoryMessenger;
  template <class T> class SmartPointer;
}

#include <GyotoError.h>
#include <stddef.h>
#include <iostream>
#include <typeinfo>
#include <string>
#include <vector>

/**
 * \brief Can be pointed to by a SmartPointer
 *
 * A class can be pointed to by a SmartPointer when it inherits from
 * class SmartPointee.
 *
 * The SmartPointee methods need to be public to be accessed by all
 * the SmartPointer < T > classes. However, it is a bad idea to
 * manipulate the counter directly. To protect these methods inside
 * your derive object, you can do as in the following example:
 *
 * @code
 * class Gyoto::Metric : protected Gyoto::SmartPointee
 * {
 *    friend class Gyoto::SmartPointer<Gyoto::Metric::Generic>;
 *    ...
 * };
 * @endcode
 *
 */
class Gyoto::SmartPointee
{
 private:
  int refCount; ///< Reference counter.

# ifdef HAVE_PTHREAD
  /**
   * When compiled with libpthread
   */
  pthread_mutex_t mutex_; ///< A mutex
#endif

 public:
  SmartPointee () ;
  virtual ~SmartPointee() ;
  SmartPointee (const   SmartPointee&) ; ///< Copy constructor
  void incRefCount () ; ///< Increment the reference counter. Warning: Don't mess with the counter.
  int decRefCount () ;  ///< Decrement the reference counter and return current value. Warning: Don't mess with the counter.
  int getRefCount () ;  ///< Get the current number of references
  /**
   * Various classes need to provide a subcontractor to be able to
   * instantiate themselves upon order from the Factory. A
   * subcontractor is a function (often a static member function)
   * which accepts a pointer to a FactoryMessenger as unique
   * parameter, communicates with the Factory using this messenger to
   * read an XML description of the object to build, and returns this
   * objet. SmartPointee::Subcontractor_t* is just generic enough a
   * typedef to cast to and from other subcontractor types:
   * Astrobj::Subcontractor_t, Metric::Subcontractor_t,
   * Spectrum::Subcontractor_t. A subcontractor needs to be registered
   * using the relevant Register() function: Astrobj::Register(),
   * Metric::Register(), Spectrum::Register().
   */
  typedef Gyoto::SmartPointer<Gyoto::SmartPointee>
    Subcontractor_t(Gyoto::FactoryMessenger*, std::vector<std::string> const &);
  ///< A subcontractor builds an object upon order from the Factory

};


/**
 * \brief Pointers performing reference counting
 *
 * Pointee must inherit from class SmartPointee.
 *
 * To create an object and a SmartPointer pointing to it:
 *
 * \code
 * SmartPointer<Gyoto::Metric::Generic> ObjPtr (new Gyoto::Metric(...));
 * \endcode
 *
 * \tparam T Sub-class of Gyoto::SmartPointee.
 */
template< class T >
class Gyoto::SmartPointer
{
 private:

  /**
   * \brief Real pointer, don't mess with it.
   */
  T *obj;

 private:
  /**
   * \brief Decrement the reference counter. Warning: don't mess with it.
   */
  void decRef ()
  {
    if (obj && obj->decRefCount() == 0) {
#     if GYOTO_DEBUG_ENABLED
      GYOTO_DEBUG_EXPR(obj);
#     endif
      delete obj;
      obj = NULL;
    }
  }
  
 public:
  /**
   * \brief Constructor from a standard pointer-to-class
   *
   * \param orig : a pointer to an instance of class T, created using new T().
   *
   * Example:
   * \code 
   * SmartPointer<Gyoto::Metric::Generic> ObjPtr (new Gyoto::Metric(...)); // create SmartPointer ObjPtr
   * \endcode
   */
  SmartPointer (T *orig = NULL) : obj(orig)
  {
    if (obj)
      obj->incRefCount();
  } 

  /**
   * \brief Copy constructor from same type
   *
   * \param orig : a SmartPointer to an instance of class T
   *
   * Example:
   * \code 
   * SmartPointer<Gyoto::Metric::Generic> ObjPtr (new Gyoto::Metric(...)); // create SmartPointer ObjPtr
   * SmartPointer<Gyoto::Metric::Generic> ObjPtr2 = ObjPtr; // create SmartPointer ObjPtr2
   * \endcode
   *
   * ObjPtr and ObjPtr2 point to the same instance of class T. Copying
   * increments the reference counter.
   */
  SmartPointer (const SmartPointer< T > &orig)
    {
      obj = orig.obj;
      if (obj)
	obj->incRefCount ();
    }

  /**
   * \brief Copy constructor from compatible type (used for casting)
   *
   * \param orig : a SmartPointer to an instance of another class U
   *
   * Example: MetricPtr is a SmartPoiter<Metric>, but really points to
   * an instance of the child class Gyoto::Kerr:
   *
   * \code 
   * SmartPointer<Gyoto::Kerr> KerrPtr (MetricPtr);
   * \endcode
   *
   * MetricPtr and KerrPtr point to the same instance of class
   * Kerr. The methods specific to class Kerr are available only to
   * KerrPtr.
   */
  template<class U>
    SmartPointer(const SmartPointer<U>& orig)
    {
      obj = dynamic_cast<T*>(const_cast<U*>(orig()));
      if (obj)
	obj->incRefCount ();
    }

  /**
   * \brief Dereference operator "*"
   *
   * \return address of the pointed-to-object
   */
  T& operator* ()
    {
      if (!obj)
	Gyoto::throwError("Null Gyoto::SmartPointer dereference in operator*");
      return *obj;
    }

  /**
   * \brief Dereference operator "*"
   *
   * \return address of the pointed-to-object
   */
  const T& operator* () const
    {
      if (!obj)
	Gyoto::throwError("Null Gyoto::SmartPointer dereference in operator*");
      return *obj;
    }

  /**
   * \brief Dereference operator "->"
   *
   * Access to the pointed-to-object's members.
   */
  T* operator-> ()
    {
      if (!obj)
	Gyoto::throwError("Null Gyoto::SmartPointer dereference in operator->");
      return obj;
    }

  /**
   * \brief Dereference operator "->" (const)
   *
   * Access to the pointed-to-object's members.
   */
  T* operator-> () const
    {
      if (!obj)
	Gyoto::throwError("Null Gyoto::SmartPointer dereference in operator->");
      return obj;
    }

  /**
   * \brief Comparison operator between two SmartPointer of same kind
   */
  bool operator== (const SmartPointer< T >&right) { return obj == right.obj; }
  
  /**
   * \brief Comparison operator between two SmartPointer of same kind
   */
  bool operator!= (const SmartPointer< T >&right) { return obj != right.obj; }
  
  /**
   * \brief Copy a SmartPointer to another (already defined)
   * SmartPointer of same kind
   */
  SmartPointer< T > &operator= (SmartPointer< T > &right)
    {
      if (this == &right)
	return *this;
      
      if (right.obj)
	right.obj->incRefCount ();
      
      decRef ();
      
      obj = right.obj;
      
      return *this;
    }

  /**
   * \brief Copy a normal pointer to an (already defined)
   * SmartPointer of same kind
   */
  SmartPointer< T > &operator= (T* right)
    {
      if (obj == right)
	return *this;
      
      decRef ();
      
      obj = right;

      if (obj) obj->incRefCount();
      
      return *this;
    }

  /**
   * \brief Cast SmartPointer to normal pointer
   */
  //  template<class U> operator U() { return obj; }
  operator T*() const { return obj; }

#if 0
  operator const T*() { return obj; }

  /**
   * \brief Check whether SmartPointer is valid
   * \return 0 if SmartPointer is NULL, 1 if valid.
   */
  operator bool () const { return obj != NULL; }

  /**
   * \brief Check whether SmartPointer is valid
   * \return 1 if SmartPointer is NULL, 0 if valid.
   */
  bool operator! () const { return obj == NULL; }
#endif

  ~SmartPointer< T > () { decRef(); }

 public:
  /**
   * \brief Get standard, non-smart pointer to object. Use with care.
   *
   * This public method is needed to cast from a SmartPointer flavor
   * to another. It should almost certainly never be used in user
   * code, except perhaps for debugging purposes.
   *
   * \return usual pointer to object, T*.
   *
   */
  //  const T* address() const { return obj; }
  const T* operator()() const { return obj; }

};

#endif
