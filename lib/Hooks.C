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

#include "GyotoHooks.h"
#include "GyotoDefs.h"
#include "GyotoUtils.h"
#include <iostream>
#include <cstddef>

using namespace Gyoto;
using namespace Gyoto::Hook;

//////// LOCAL TO THIS FILE ////////
/* LISTENER LIST */
class Gyoto::Hook::Teller::ListenerItem {
  friend  class Gyoto::Hook::Teller;
  friend class Gyoto::Hook::Listener;
 protected:
  Listener * listener;
  ListenerItem * next;
 public:
  ListenerItem(Listener*, ListenerItem*);
  virtual ~ListenerItem();
  virtual void tell(Gyoto::Hook::Teller *);
  size_t len() const;
};

Teller::ListenerItem::ListenerItem(Listener*hear, ListenerItem * nxt)
  : listener(hear), next(nxt) {}
Teller::ListenerItem::~ListenerItem() {if (next) delete next;}

void Teller::ListenerItem::tell(Teller * teller) {
  GYOTO_DEBUG<<teller<<" telling to " <<listener<< std::endl;
  listener->tell(teller);
  if (next) next->tell(teller);
}

size_t Teller::ListenerItem::len() const {
  size_t l = 1;
  if (next) l += next->len();
  return l;
}
/////////////////////////////////////


/* LISTENER */
Listener::Listener() {}
Listener::~Listener() {}
void Listener::tell(Teller *) {}

/* TELLER */
Teller::Teller() : listeners_(0) {}
Teller::Teller(const Teller &) : listeners_(0) {}
Teller::~Teller() {if (listeners_) delete listeners_;}
void Teller::hook(Listener * hear)
{
  GYOTO_DEBUG <<"new hook: "<< hear <<" for Teller " << this << std::endl;
  listeners_ = new ListenerItem(hear, listeners_);
}

void Teller::unhook(Listener * hear) {
  GYOTO_DEBUG <<"removing hook: "<< hear <<" from Teller " << this << std::endl;
  ListenerItem * item = listeners_;
  ListenerItem ** parent = &listeners_;
  while (item) {
    if (item->listener==hear) {
      *parent=item->next;
      item->next=0;
      delete item;
      item=*parent;
    } else {
      parent=&(item->next);
      item=item->next;
    }
  }
}

void Teller::tellListeners() { if (listeners_) listeners_->tell(this); }


