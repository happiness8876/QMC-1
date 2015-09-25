#ifndef WORM_HPP_CMC
#define WROM_HPP_CMC
#include "Node.hpp"

template <class Temp> class GeneralWorm
{
  typedef typename Node<Temp>::Iter Iter;

  public:
    void set (const Iter& it, bool up) { _self = it; _up = up; }
    const Iter& it () const { return _self; }
    bool up () const { return _up; }
    void it (const Iter& it_) { _self = it_; }
    void up (bool up_) { _up = up_; }
    void bounce () { _up = !_up;
#ifdef DEBUG_TRACK_MODE
__DEBUG_TRACK
std::cout << (_up ? "  up\n" : "  down\n");
#endif
    }

  private:
    bool _up;
    Iter _self;
};
#endif
