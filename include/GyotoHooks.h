#ifndef __GyotoHooks_h
#define __GyotoHooks_h

#include <cstddef>

namespace Gyoto {
  namespace Hook {
    class Teller;
    class Listener;
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
 private:
  class ListenerItem;
  ListenerItem *listeners_;
 public:
  Teller();
  Teller(const Teller &);
  ~Teller();
  virtual void hook (Listener *);
  virtual void unhook (Listener *);
 protected:
  virtual void tellListeners();
};

#endif
