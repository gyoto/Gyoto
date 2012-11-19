#ifndef __GyotoHooks_h
#define __GyotoHooks_h

#include <cstddef>

namespace Gyoto {
  namespace Hook {
    class Teller;
    class Listener;
    class ListenerItem;
  };
};

class Gyoto::Hook::Listener {
  friend class Gyoto::Hook::Teller;
 public:
  Listener();
  ~Listener();
  virtual void tell(Gyoto::Hook::Teller *);
};

class Gyoto::Hook::Teller {
  friend class Gyoto::Hook::Listener;
 protected:
  ListenerItem *listeners;
 public:
  Teller();
  Teller(const Teller &);
  ~Teller();
  virtual void hook (Listener *);
  virtual void unhook (Listener *);
};

class Gyoto::Hook::ListenerItem {
  friend  class Gyoto::Hook::Teller;
  friend class Gyoto::Hook::Listener;
 protected:
  Listener * listener;
  ListenerItem * next;
 public:
  ListenerItem(Listener*, ListenerItem*);
  ~ListenerItem();
  virtual void tell(Gyoto::Hook::Teller *);
  size_t len() const;
};

#endif
