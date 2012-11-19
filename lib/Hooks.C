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

using namespace Gyoto;
using namespace Gyoto::Hook;

/* LISTENER */
Listener::Listener() {}
Listener::~Listener() {}
void Listener::tell(Teller *) {};

/* TELLER */
Teller::Teller() : listeners(0) {}
Teller::Teller(const Teller &o) : listeners(0) {}
Teller::~Teller() {if (listeners) delete listeners;}
void Teller::hook(Listener * hear)
{listeners = new ListenerItem(hear, listeners);}

void Teller::unhook(Listener * hear) {
  ListenerItem * item = listeners;
  ListenerItem ** parent = &listeners;
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

/* LISTENER LIST */
ListenerItem::ListenerItem(Listener*hear, ListenerItem * nxt)
  : listener(hear), next(nxt) {}
ListenerItem::~ListenerItem() {if (next) delete next;}

void ListenerItem::tell(Teller * teller) {
  listener->tell(teller);
  if (next) next->tell(teller);
}

size_t ListenerItem::len() const {
  size_t l = 1;
  if (next) l += next->len();
  return l;
}
