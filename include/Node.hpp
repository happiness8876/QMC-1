#ifndef Node_HPP
#define Node_HPP
#include <list>
#include <sstream>
#include <cstdlib>
#include "Lattice.hpp"

template <class Temp> class Diagram;

template <class Temp> class Node
{
  public:
    typedef typename std::list<Node>::iterator Iter;
    // Get value:
    int site () const { return _site; }
    double time () const { return _time; }
    bool end () const { return _end; }
    const Iter& assoc (int i) const { return _assoc[i]; }
    Iter ass_bef (int i) const { return i == _conji ? _assoc[i] : --(Iter(_assoc[i])); }
    Iter ass_aft (int i) const { return (i == _conji || _end) ? ++(Iter(_assoc[i])) : _assoc[i]; }
    const Iter& conj () const { return _assoc[_conji]; }
    bool is_hopping () const { return _conji != -1; }
    Temp z;

    const std::string info (const Lattice& latt) const;
    const std::string self_info () const;

  private:
    double _time;
    int _site, _conji;
    Iter _assoc[6]; // link to the neighbor nodes just upper or at the same time
    bool _end; // indicate whether this is an end-node

    Node () {}
    Node (int s, double t, bool end = false) { _site = s; _time = t; _end = end; _conji = -1; }

    friend class Diagram<Temp>;
};

template <class Temp>
const std::string Node<Temp> :: info (const Lattice& latt) const
{
  std::stringstream sstr;
  sstr << "self:\n";
  sstr << "  x, y, z = " << latt(_site).xi() << ", " << latt(_site).yi() << ", " << latt(_site).zi() << "\n";
  sstr << this->self_info();
  for(int i = 0; i < latt(_site).nbs(); i++) {
    sstr << "assoc " << i << ":\n";
    sstr << _assoc[i]->self_info();
  }
  if (_conji != -1) {
    sstr << "conjugate:\n";
    sstr << conj()->self_info();
  }
  return sstr.str();
}

template <class Temp>
const std::string Node<Temp> :: self_info () const
{
  std::stringstream sstr;
  sstr << "  address = " << this << "\n";
  sstr << "  site = " << _site << "\n";
  sstr << "  time = " << _time << "\n";
  sstr << z << "\n";
  return sstr.str();
}
#endif
