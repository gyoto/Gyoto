/*
    Copyright 2011-2012, 2014, 2020 Thibaut Paumard

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

#include <cstddef>
#include <GyotoSmartPointer.h>
#include <GyotoError.h>
#include <string>
#include <cstring>

Gyoto::SmartPointee::SmartPointee() :
  refCount (0)
{
#ifdef HAVE_PTHREAD
  pthread_mutex_init(&mutex_, NULL);
#endif
}

Gyoto::SmartPointee::~SmartPointee() {
  GYOTO_DEBUG << typeid(*this).name() << ": refCount=" << refCount << std::endl;
}

Gyoto::SmartPointee::SmartPointee(const SmartPointee&) :
  refCount (0)
{
#ifdef HAVE_PTHREAD
  pthread_mutex_init(&mutex_, NULL);
#endif
}

Gyoto::SmartPointee& Gyoto::SmartPointee::operator=(const SmartPointee& right){
  // Check for self-assignment
  if (this == &right) {
    return *this;
  }

  // initialize the refcound to 0: the two objects do not use the same
  // reference counter
  refCount = 0;

  // Reinitialize the mutex if needed
#ifdef HAVE_PTHREAD
  pthread_mutex_destroy(&mutex_); // Destroy the old mutex
  pthread_mutex_init(&mutex_, NULL); // Initialize a new mutex
#endif

  // Return a reference to the current object
  return *this;
}

void Gyoto::SmartPointee::incRefCount () {
#ifdef HAVE_PTHREAD
 pthread_mutex_lock(&mutex_);
#endif
 refCount++;
 GYOTO_DEBUG << typeid(*this).name() << ": refCount=" << refCount << std::endl;
#ifdef HAVE_PTHREAD
 pthread_mutex_unlock(&mutex_);
#endif
}
int Gyoto::SmartPointee::decRefCount () {
#ifdef HAVE_PTHREAD
 pthread_mutex_lock(&mutex_);
 int n = --refCount;
 GYOTO_DEBUG << typeid(*this).name() << ": refCount=" << refCount << std::endl;
 pthread_mutex_unlock(&mutex_);
 return n;
#else
 GYOTO_DEBUG << typeid(*this).name() << ": refCount=" << refCount-1 << std::endl;
 return --refCount;
#endif
}
int Gyoto::SmartPointee::getRefCount() {
#ifdef HAVE_PTHREAD
 pthread_mutex_lock(&mutex_);
 int n = refCount;
 pthread_mutex_unlock(&mutex_);
 return n;
#else
 return refCount;
#endif
}
