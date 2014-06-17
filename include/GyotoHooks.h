/**
 * \file GyotoHooks.h 
 * \brief Tellers tell Listeners when they mutate
   *
   * A Listener can hook() to a Teller. The Teller will tell it when
   * it mutates using Listener::tell(), usually through the highter
   * lever Teller::tellListeners(). The Listener can later
   * unhook(). The Listener must therefore implement Listener::tell()
   * to react when it is told.
 */

/*
    Copyright 2012-2013 Thibaut Paumard

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

#ifndef __GyotoHooks_h
#define __GyotoHooks_h

/**
 * \namespace Gyoto::Hook
 * \brief Listeners attach to Tellers
 *
 * A Listener can hook() to a Teller. The Teller will tell it when
 * it mutates using Listener::tell(), usually through the highter
 * lever Teller::tellListeners(). The Listener can later
 * unhook(). The Listener must therefore implement Listener::tell()
 * to react when it is told.
 */
namespace Gyoto {
  namespace Hook {
    /**
     * \class Teller
     * \brief Listen to me and I'll warn you when I change
     *
     * Listen to me by calling my hook() method.
     */
    class Teller;

    /**
     * \class Listener
     * \brief I might listen to a Teller
     *
     * Whisper to my ear by using my tell() method.
     */
    class Listener;
  }
}

class Gyoto::Hook::Listener {
  friend class Gyoto::Hook::Teller;
 public:
  Listener(); ///< Constructor
  virtual ~Listener(); ///< Destructor

 protected:
  /**
   * \brief This is how a Teller tells
   *
   * A teller will basically call listener->tell(this).
   *
   * \param msg Teller* the Teller who is telling... Useful if the
   * Listener listens to several Tellers.
   */
  virtual void tell(Gyoto::Hook::Teller *msg);
};

class Gyoto::Hook::Teller {
  friend class Gyoto::Hook::Listener;
 private:
  /**
   * \class ListenerItem
   * \brief Private (undocumented) class to hold listeners_
   */
  class ListenerItem;

  /**
   * \brief Linked list of Listener items
   */
  ListenerItem *listeners_;

 public:
  Teller(); ///< Default constructor
  Teller(const Teller &); ///< Copy constructor
  virtual ~Teller(); ///< Destructor

  /**
   * \brief Start listening
   *
   * Use from a Hook::Listener object method:
   * \code
   * teller->hook(this)
   * \endcode
   * where "this" is a Listener and "teller" is a Teller.
   *
   * Use unhook() later to stop listening to a given Teller. 
   *
   * \param listener pointer to the new listener
   */
  virtual void hook (Listener * listener);

  /**
   * \brief Stop listening
   *
   * Use from a Hook::Listener object method:
   * \code
   * teller->unhook(this)
   * \endcode

   * where "this" is a Listener, "teller" is a Teller, and "this" has
   * called teller->hook(this) previously.
   *
   * \param listener pointer to the listener
   */
  virtual void unhook (Listener * listener);

 protected:
  /**
   * \brief Call tell() on each hooked Listener
   *
   * Whenever a Teller mutates, it should warn any Listener hooked to
   * it using tellListeners().
   */
  virtual void tellListeners();
};

#endif
