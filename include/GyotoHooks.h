/**
 * \file GyotoHooks.h 
 * \brief Tellers tell Listeneres when they mutate
 */

/*
    Copyright 2012 Thibaut Paumard

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

namespace Gyoto {
  /**
   * \namespace Hook
   * \brief Listeners attach to Tellers
   */
  namespace Hook {
    /**
     * \class Teller
     * \brief Listen to me and I'll warn you when I change
     */
    class Teller;

    /**
     * \class Listener
     * \brief I might listen to a Teller
     */
    class Listener;
  };
};

class Gyoto::Hook::Listener {
  friend class Gyoto::Hook::Teller;
 public:
  Listener();
  ~Listener();

 protected:
  /**
   * \brief This is how a Teller tells
   *
   * A teller will basically call listener->tell(this)
   *
   * \param msg Teller* the Teller who is telling... Useful if the
   * Listener listens to several Tellers.
   */
  virtual void tell(Gyoto::Hook::Teller *msg); ///< Called by Teller
};

class Gyoto::Hook::Teller {
  friend class Gyoto::Hook::Listener;
 private:
  class ListenerItem;

  /**
   * \brief Linked list of Listener items
   */
  ListenerItem *listeners_;

 public:
  Teller();
  Teller(const Teller &);
  ~Teller();

  /**
   * \brief Start listening
   *
   * Use as teller->hook(this)
   *
   * \param listener pointer to the new listener
   */
  virtual void hook (Listener *);

  /**
   * \brief Stop listening
   *
   * Use as teller->hook(this)
   *
   * \param listener pointer to the new listener
   */
  virtual void unhook (Listener *);

 protected:
  /**
   * \brief Call tell() on each hooked Listener 
   */
  virtual void tellListeners();
};

#endif
