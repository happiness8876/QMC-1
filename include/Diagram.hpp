#ifndef DIAGRAM_HPP_CMC
#define DIAGRAM_HPP_CMC
#include <cmath>
#include <fstream>
#include <sstream>
#include <limits>
#include <typeinfo>
#include <vector>
#include "MotherLattice.hpp"
#include "Node.hpp"

template <class Temp> class Diagram
{
  public:
    typedef typename std::list< Node<Temp> >::iterator Iter;
    typedef typename std::list< Node<Temp> >::const_iterator const_Iter;
    static typename std::list< Node<Temp> > NULL_LIST;
    const static typename std::list< Node<Temp> >::iterator NULL_ITER () { return NULL_LIST.end(); }

    Diagram (int Lx, int Ly, int Lz, double Beta, bool BCx=0, bool BCy=0, bool BCz=0);

    void diag_insert (const Iter& it_goal, double time, Iter& it_bef, Iter& it_aft);
    void diag_insert (int site, double time, Iter& it_bef, Iter& it_aft);
    void diag_walk (const Iter& it_worm, double time, bool up);
    void diag_hop (const Iter& it_worm, int nbi, bool up, Iter& bef, Iter& aft);
    void diag_delete_hop (Iter& it_worm, bool up);
    void diag_relink_hop (Iter& it_worm, int nbi, bool up, Iter& bef, Iter& aft);
    void diag_cross_end (const Iter& it_worm, bool up);
    void diag_cross_node (const Iter& it_worm, bool up);
    void pass_neighbor (bool up, const Iter& it_worm, Iter& it_nb, int nbi);
    void diag_halt (const Iter& it_worm, bool up);
    void diag_remove (Iter& it_worm, bool up);

    // Information
    int size () const { return SITEN; }
    const Iter& beg (int site) const { return it_beg[site]; }
    const Iter& last (int site) const { return it_end[site]; }
    int hopdir (const const_Iter& it) const;

    // Debug
    std::string info () const;
    std::string info (const Node<Temp>& node) const;
    void diag_info_plot (std::string filename, std::string ext="diag") const;
    virtual void info_plot_line (const const_Iter& it, std::ofstream& ofs) const;
    void info_plot_hop (const const_Iter& it, std::ofstream& ofs) const;
    void check (int site, std::string funcname="unknown") const;
    void check_neighbors (int site, std::string funcname) const;
    void do_check (int site, std::ofstream& ofs, const std::string& funcname) const;

    virtual void write (std::string name);
    virtual void read (std::string name);
    MotherLattice _latt;

  protected:
    int SITEN;
    double BETA;
    std::vector <std::list< Node<Temp> > > nodeList;
    std::vector<Iter> it_beg, it_end;

    Iter find_aft (int site, double time);
    void insert2 (const Iter& goal, double t, Iter& bef, Iter& aft);
    inline void remove2 (const Iter& it1, const Iter& it2, const Iter& aft);

    void link (int nbi, const Iter& it);
    void unlink (const Iter& it, const Iter& it_aft, int nbi);
    void set_conj (const Iter& it, int nbi, Iter& it_conj); 
};
template <class Temp> typename std::list< Node<Temp> > Diagram<Temp>::NULL_LIST;

template <class Temp>
Diagram<Temp> :: Diagram (int Lx, int Ly, int Lz, double Beta, bool BCx, bool BCy, bool BCz) : SITEN (Lx*Ly*Lz), BETA (Beta), _latt (Lx, Ly, Lz, BCx, BCy, BCz)
{
  nodeList.resize (SITEN);
  it_beg.resize (SITEN);
  it_end.resize (SITEN);
  // Insert terminal nodeList
  for(int site = 0; site < SITEN; ++site) {
    Node<Temp> node_head (site, 0., true);
    Node<Temp> node_tail (site, BETA, true);
    nodeList[site].push_back (node_head);
    nodeList[site].push_back (node_tail);
    it_beg[site] = nodeList[site].begin();
    it_end[site] = it_beg[site];  ++it_end[site];
  }
  // Links the nodeList
  Iter it, it_nb;
  for(int site = 0; site < SITEN; ++site) {
    for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
      // find the neighbor site
      int neighbor = _latt(site).nb(nbi);
      it = nodeList[site].begin();
      it_nb = nodeList[neighbor].begin();
      // link the head nodeList
      it->_assoc[nbi] = it_nb;
      // link the tail nodeList
      (++it)->_assoc[nbi] = (++it_nb);
    }
  }
#ifdef DEBUG_MODE
  for(int site = 0; site < SITEN; ++site)
    check (site, __func__);
#endif
}

template <class Temp>
inline void Diagram<Temp> :: diag_insert (const Iter& it_goal, double time, Iter& it_bef, Iter& it_aft)
{
  insert2 (it_goal, time, it_bef, it_aft);
#ifdef DEBUG_MODE
  check_neighbors (it_goal->site(), __func__);
#endif
}

template <class Temp>
inline void Diagram<Temp> :: diag_insert (int site, double time, Iter& it_bef, Iter& it_aft)
{
  Iter it_goal = find_aft (site, time);
  insert (it_goal, time, it_bef, it_aft);
}

template <class Temp>
void Diagram<Temp> :: diag_remove (Iter& it_worm, bool up)
{
  Iter it_fix, it_aft;
  if (up) {
    it_fix = it_worm;
    ++it_fix;
    it_aft = it_fix;
    ++it_aft;
  }
  else {
    it_fix = it_worm;
    --it_fix;
    it_aft = it_worm;
    ++it_aft;
  }
  remove2 (it_worm, it_fix, it_aft);
#ifdef DEBUG_MODE
  int site = it_worm->site();
  check_neighbors (site, __func__);
#endif
}

template <class Temp>
void Diagram<Temp> :: diag_walk (const Iter& it_worm, double time, bool up)
{
  // Set time
  it_worm->_time = time;

  // Reset the linking
  int site = it_worm->_site;
  if (up) {
    for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
      const int& nb = _latt(site).nb(nbi);
      const int& opp = _latt(site).opp(nbi);
      Iter it_ass = it_worm->_assoc[nbi];
      while (it_ass->_time < it_worm->_time)
        (it_ass++)->_assoc[opp] = it_worm;
      it_worm->_assoc[nbi] = it_ass;
    }
  }
  else {
    Iter it_aft = it_worm;  ++it_aft;
    for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
      const int& nb = _latt(site).nb(nbi);
      const int& opp = _latt(site).opp(nbi);
      Iter it_ass = it_worm->_assoc[nbi];  --it_ass;
      while (it_ass->_time > it_worm->_time)
        (it_ass--)->_assoc[opp] = it_aft;
      it_worm->_assoc[nbi] = ++it_ass;
    }
  }
#ifdef DEBUG_MODE
  check_neighbors (site, __func__);
#endif
}

template <class Temp>
void Diagram<Temp> :: diag_hop (const Iter& it_worm, int nbi, bool up, Iter& it_bef, Iter& it_aft)
{
  // Insert nodes
  Iter it_goal = it_worm->assoc(nbi);
  insert2 (it_goal, it_worm->time(), it_bef, it_aft);

  // Set special associate
  int opp = _latt(it_worm->site()).opp(nbi);
  while (it_bef->assoc(opp) != it_worm) // same time issue
    (it_bef->_assoc[opp]++)->_assoc[nbi] = it_bef;
  if (up) {
    set_conj (it_worm, nbi, it_bef);
    it_aft->_assoc[opp] = ++Iter(it_worm);
  }
  else {
    set_conj (it_worm, nbi, it_aft);
    it_bef->_assoc[opp] = it_worm;
  }
#ifdef DEBUG_MODE
  int site = it_worm->site();
  check_neighbors (site, __func__);
  check_neighbors (it_bef->site(), __func__);
#endif
}

template <class Temp>
void Diagram<Temp> :: diag_delete_hop (Iter& it_worm, bool up)
{
#ifdef DEBUG_MODE
int oldsite = it_worm->site();
#endif
  Iter it_aft, it_halt, it_conj;
  if (up) {
    it_halt = it_worm;
    ++it_halt;
    it_aft = it_halt;
    ++it_aft;
    it_conj = it_halt->conj();
  }
  else {
    it_halt = it_worm;
    --it_halt;
    it_aft = it_worm;
    ++it_aft;
    it_conj = it_halt->conj();
  }

  // Delete nodes
  remove2 (it_worm, it_halt, it_aft);

  // Remove conjugate
  it_conj->_assoc[it_conj->_conji] = it_aft;
  it_conj->_conji = -1;

  it_worm = it_conj;
#ifdef DEBUG_MODE
  int site = it_worm->site();
  check_neighbors (oldsite, __func__);
  check_neighbors (site, __func__);
#endif
}

template <class Temp>
void Diagram<Temp> :: diag_relink_hop (Iter& it_worm, int nbi, bool up, Iter& bef, Iter& aft)
{
#ifdef DEBUG_MODE
int oldsite = it_worm->site();
#endif
  diag_delete_hop (it_worm, up);
  diag_hop (it_worm, nbi, !up, bef, aft);
#ifdef DEBUG_MODE
  int site = it_worm->site();
  check_neighbors (oldsite, __func__);
  check_neighbors (site, __func__);
  check_neighbors (bef->site(), __func__);
#endif
}

template <class Temp>
void Diagram<Temp> :: pass_neighbor (bool up, const Iter& it_worm, Iter& it_nb, int nbi)
{
//std::cout << "ppp\n";
  it_worm->_time = it_nb->time();
  if (up) {
    ++(it_worm->_assoc[nbi]);
    int opp = _latt(it_worm->site()).opp (nbi);
    it_nb->_assoc[opp] = it_worm;
  }
  else {
    --(it_worm->_assoc[nbi]);
    int opp = _latt(it_worm->site()).opp (nbi);
    ++(it_nb->_assoc[opp]);
  }
}

template <class Temp>
void Diagram<Temp> :: diag_cross_end (const Iter& it_worm, bool up)
{
  int site = it_worm->site();
  if (up) {
    Iter it_aft = it_beg[site];
    ++it_aft;
    for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
      // Unlink
      unlink (it_worm, it_end[site], nbi);
      // Link
      Iter it_ass = it_beg[site]->_assoc[nbi];
      it_worm->_assoc[nbi] = ++it_ass;
    }
    // Move node
    nodeList[site].splice (it_aft, nodeList[site], it_worm);
    it_worm->_time = 0.;
  }
  else {
    Iter it_aft = it_worm;  ++it_aft;
    for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
      // Unlink
      unlink (it_worm, it_aft, nbi);
      unlink (it_end[site], it_worm, nbi);
      // Link
      int nbsite = _latt.nb (site, nbi);
      it_worm->_assoc[nbi] = it_end[nbsite];
    }
    // Move node
    nodeList[site].splice (it_end[site], nodeList[site], it_worm);
    it_worm->_time = BETA;
  }
#ifdef DEBUG_MODE
  check_neighbors (site, __func__);
#endif
}

template <class Temp>
void Diagram<Temp> :: diag_cross_node (const Iter& it_worm, bool up)
{
  int site = it_worm->site();
  if (up) {
    Iter it_aft = it_worm;
    ++it_aft;
    for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
      // Unlink
      unlink (it_worm, it_aft, nbi);
      // Link
      it_worm->_assoc[nbi] = it_aft->ass_aft(nbi);
    }
    // Move node
    it_worm->_time = it_aft->_time;
    nodeList[site].splice (++it_aft, nodeList[site], it_worm);
  }
  else {
    Iter it_bef = it_worm, it_aft = it_worm;
    ++it_aft;
    --it_bef;
    Iter it_bef_conj = it_bef->conj();
    for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
      // Unlink
      unlink (it_worm, it_aft, nbi);
      unlink (it_bef, it_worm, nbi);
      // Link
      it_worm->_assoc[nbi] = it_bef->_assoc[nbi];
    }
    // Set conjugate back
    set_conj (it_bef, it_bef->_conji, it_bef_conj);
    // Move node
    nodeList[site].splice (it_bef, nodeList[site], it_worm);
    it_worm->_time = it_bef->_time;
  }
#ifdef DEBUG_MODE
  check_neighbors (site, __func__);
#endif
}

template <class Temp>
void Diagram<Temp> :: diag_halt (const Iter& it_worm, bool up)
{
  Iter it_halt = up ? ++Iter(it_worm) : --Iter(it_worm);
  diag_walk (it_worm, it_halt->time(), up);
#ifdef DEBUG_MODE
  int site = it_worm->site();
  check_neighbors (site, __func__);
#endif
}

template <class Temp>
void Diagram<Temp> :: insert2 (const Iter& goal, double t, Iter& bef, Iter& aft)
{
  int site = goal->site();
  // Insert nodes:
  bef = nodeList[site].insert (goal, Node<Temp>(site, t));
  aft = nodeList[site].insert (goal, Node<Temp>(site, t));
  // Set associate:
  for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
    link (nbi, bef);
    aft->_assoc[nbi] = bef->_assoc[nbi];
  }
}

template <class Temp>
void Diagram<Temp> :: remove2 (const Iter& it1, const Iter& it2, const Iter& aft)
{
  int site = aft->_site;
  // Unlink the node:
  for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
    unlink (it1, aft, nbi);
    unlink (it2, aft, nbi);
  }
  // Remove
  nodeList[site].erase (it1);
  nodeList[site].erase (it2);
}

template <class Temp>
void Diagram<Temp> :: link (int nbi, const Iter& it)
{
  const int& opp = _latt(it->_site).opp(nbi);
  Iter it_bef = it;  --it_bef;
  Iter ass = it_bef->ass_aft(nbi);
  while (ass->_time < it->_time) {
    ass->_assoc[opp] = it;
    ++ass;
  }
  it->_assoc[nbi] = ass;
}

template <class Temp>
inline typename Diagram<Temp>::Iter Diagram<Temp> :: find_aft (int site, double t)
{
  Iter it = nodeList[site].begin();
  while (it->_time <= t) ++it;
  return it;
}

template <class Temp>
void Diagram<Temp> :: set_conj (const Iter& it, int nbi, Iter& it_conj)
{
  // Set _conji
  it->_conji = nbi;
  int opp = _latt(it->_site).opp(nbi);
  it_conj->_conji = opp;
  // Set assoc[_conji]
  it->_assoc[nbi] = it_conj;
  it_conj->_assoc[opp] = it;
}

template <class Temp>
void Diagram<Temp> :: unlink (const Iter& it, const Iter& it_aft, int nbi)
{
  Iter ass = it->ass_bef(nbi);
  int opp = _latt(it->_site).opp(nbi);
  while (ass->_assoc[opp] == it) {
    ass->_assoc[opp] = it_aft;
    --ass;
  }
}

template <class Temp>
int Diagram<Temp> :: hopdir (const const_Iter& it) const
{
  if (it->_conji == -1) return 0;
  int dir, di, dj, L;
  int i = it->site(), j = it->conj()->site();
  if (_latt(i).xi() != _latt(j).xi()) {
    // x direction
    dir = 1;
    di = _latt(i).xi();
    dj = _latt(j).xi();
    L = _latt.Lx();
  }
  else if (_latt(i).yi() != _latt(j).yi()) {
    // y direction
    dir = 2;
    di = _latt(i).yi();
    dj = _latt(j).yi();
    L = _latt.Ly();
  }
  else if (_latt(i).zi() != _latt(j).zi()) {
    // z direction
    dir = 3;
    di = _latt(i).zi();
    dj = _latt(j).zi();
    L = _latt.Lz();
  }
  if (di < dj)
    dir *= -1;
  if ((di == 0 && dj == L) || (di == L && dj == 0))
    // crossing boundary
    dir *= -1;
  return dir;
}

template <class Temp>
std::string Diagram<Temp> :: info (const Node<Temp>& node) const
{
  std::stringstream sstr;
  sstr << "location: " << &node << "\n"
       << "site = " << node._site << "\n"
       << "time = " << node._time << "\n";
       //<< "n bef = " << node.nbef() << "\n";
  for(int i = 0; i < _latt(node._site).nbs(); ++i) {
    sstr << "assoc_memory[" << i << "]: " << &(node.assoc(i)) << "\n";
    sstr << "assoc_site [" << i << "] = " << node.assoc(i)->_site << "\n";
    sstr << "assoc_time [" << i << "] = " << node.assoc(i)->_time << "\n";
  }
  return sstr.str();
}

template <class Temp>
std::string Diagram<Temp> :: info () const
{
  std::stringstream sstr;
  Iter it;
  sstr << " --- Information of Worldline Diagram ---\n\n";
  for(int is = 0; is < SITEN; ++is) {
    sstr << "site " << is << ":\n\n";
    for(it = nodeList[is].begin(); it != nodeList[is].end(); ++it) {
      sstr << info(*it) << "\n";
      }
    }
  sstr << " --- End of Worldline Information ---\n\n";
  return sstr.str();
}

template <class Temp>
void Diagram<Temp> :: diag_info_plot (std::string filename, std::string ext) const
{
  std::ofstream ofs ((filename+".line."+ext).c_str());
  std::ofstream ofs2 ((filename+".hop."+ext).c_str());
  if (!ofs) {
    std::cout << __FILE__ << " :: " << __func__ << " : cannot open file '" << (filename+".line."+ext) << "'\n";
    exit(100);
  }
  if (!ofs2) {
    std::cout << __FILE__ << " :: " << __func__ << " : cannot open file '" << (filename+".hop."+ext) << "'\n";
    exit(100);
  }
  for(int i = 0; i < SITEN; ++i) {
    ofs << _latt(i).xi() << " " << _latt(i).yi() << " " << _latt(i).zi();
    for(const_Iter it = nodeList[i].begin(); it != nodeList[i].end(); ++it) {
      // Write line information
      info_plot_line (it, ofs);
      // Write hopping information
      info_plot_hop (it, ofs2);
    }
    ofs << "\n";
  }
  ofs.close();
  ofs2.close();
}

template <class Temp>
void Diagram<Temp> :: info_plot_line (const const_Iter& it, std::ofstream& ofs) const
{
  ofs << " " << it->time();
}

template <class Temp>
void Diagram<Temp> :: info_plot_hop (const const_Iter& it, std::ofstream& ofs) const
{
  int i = it->site();
  if (!(it->is_hopping())) {
    if (it != nodeList[i].begin() && it != it_end[i])
      // worms
      ofs << 0 << " " << _latt(i).xi() << " " << _latt(i).yi() << " " << _latt(i).zi() << " " << it->time() << "\n";
  }
  else if (it->site() > it->conj()->site()) {
    int hp = hopdir (it);
    int j = it->conj()->site();
    int dsite = hp > 0 ? 1 : -1;
    if (abs(hp) == 1) {
      // x direction
      if ((_latt(i).xi() - _latt(j).xi()) > 1)
        // cross boundary
        ofs << -1 << " " << _latt(i).xi() << " " << _latt(j).xi();
      else if ((_latt(i).xi() - _latt(j).xi()) < -1)
        // cross boundary
        ofs << -1 << " " << _latt(j).xi() << " " << _latt(i).xi();
      else
        ofs << 1 << " " << _latt(i).xi() << " " << (_latt(i).xi() - dsite);
      ofs << " " << _latt(i).yi() << " " << _latt(i).zi() << " " << it->time() << "\n";
    }
    else if (abs(hp) == 2) {
      // y direction
      if ((_latt(i).yi() - _latt(j).yi()) > 1)
        ofs << -2 << " " << _latt(i).yi() << " " << _latt(j).yi();
      else if ((_latt(i).yi() - _latt(j).yi()) < -1)
        ofs << -2 << " " << _latt(j).yi() << " " << _latt(i).yi();
      else
        ofs << 2 << " " << _latt(i).yi() << " " << (_latt(i).yi() - dsite);
      ofs << " " << _latt(i).xi() << " " << _latt(i).zi() << " " << it->time() << "\n";
    }
    else if (abs(hp) == 3) {
      // z direction
      if ((_latt(i).zi() - _latt(j).zi()) > 1)
        ofs << -3 << " " << _latt(i).zi() << " " << _latt(j).zi();
      else if ((_latt(i).zi() - _latt(j).zi()) < -1)
        ofs << -3 << " " << _latt(j).zi() << " " << _latt(i).zi();
      else
        ofs << 3 << " " << _latt(i).zi() << " " << (_latt(i).zi() - dsite);
      ofs << " " << _latt(i).xi() << " " << _latt(i).yi() << " " << it->time() << "\n";
    }
  }
}

template <class Temp>
void Diagram<Temp> :: check_neighbors (int site, std::string funcname) const
{
  check (site, funcname);
  for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
    int nbsite = _latt.nb (site, nbi);
    check (nbsite, funcname);
  }
}

template <class Temp>
void Diagram<Temp> :: write (std::string name)
{
  _latt.write (name);
  // Rename the old bakcup-file if exist
  std::string fullname = name+".diagram";
  std::ifstream ifs (fullname.c_str());
  bool exist = ifs.good();
  ifs.close();
  std::string oldname = fullname + ".old";
  if (exist) rename (fullname.c_str(), oldname.c_str());
  std::ofstream ofs (fullname.c_str(), std::ios::binary);
  if (!ofs) {
    std::cout << "Error open file: DoubleDiagram: write: " << fullname << "\n";
    exit(100);
  }
  // Write parameters
  ofs.write((char*)&SITEN, sizeof(SITEN));
  // Write nodeList
  for(int i = 0; i < SITEN; ++i) {
    int size = nodeList[i].size();
    ofs.write((char*)&size, sizeof(size));
    for(Iter it = nodeList[i].begin(); it != nodeList[i].end(); ++it)
      ofs.write((char*)&(*it),sizeof(*it));
  }
  // Write associate
  for(int i = 0; i < SITEN; ++i)
    for(Iter it = nodeList[i].begin(); it != nodeList[i].end(); ++it) {
      for(int nbi = 0; nbi < _latt(i).nbs(); ++nbi) {
        Iter it_ass = it->assoc(nbi);
        int nbsite = _latt(i).nb(nbi);
        int dist = distance (nodeList[nbsite].begin(), it_ass);
        ofs.write((char*)&dist, sizeof(dist));
      }
    }
  ofs.close();
  // remove old backup file
  if (exist) std::remove (oldname.c_str());
}

template <class Temp>
void Diagram<Temp> :: read (std::string name)
{
  _latt.read (name);
  // Check the old backup file
  std::string fullname = name+".diagram";
  std::string oldname = fullname + ".old";
  std::ifstream ifs_old (oldname.c_str());
  bool old_exist = ifs_old.good();
  ifs_old.close();
  std::ifstream ifs;
  if (old_exist) ifs.open (oldname.c_str(), std::ios::binary);
  else ifs.open (fullname.c_str(), std::ios::binary);
  if (!ifs) {
    std::cout << "Error open file: DoubleDiagram: read: " << fullname << "\n";
    exit(100);
  }
  // Read parameters
  ifs.read((char*)&SITEN, sizeof(SITEN));
  nodeList.resize (SITEN);
  it_beg.resize (SITEN);
  it_end.resize (SITEN);
  // Read nodeList
  for(int i = 0; i < SITEN; ++i) {
    int size;  ifs.read((char*)&size, sizeof(size));
    nodeList[i].resize(size, Node<Temp>());
    for(Iter it = nodeList[i].begin(); it != nodeList[i].end(); ++it)
      ifs.read((char*)&(*it),sizeof(*it));
  }
  // Read associate
  for(int i = 0; i < SITEN; ++i) {
    for(Iter it = nodeList[i].begin(); it != nodeList[i].end(); ++it) {
      for(int nbi = 0; nbi < _latt(i).nbs(); ++nbi) {
        int nbsite = _latt(i).nb(nbi);
        int dist;  ifs.read((char*)&dist, sizeof(dist));
        Iter it_ass = nodeList[nbsite].begin();
        advance (it_ass, dist);
        it->_assoc[nbi] = it_ass;
      }
    }
    // Set it_beg[], it_end[]
    it_beg[i] = nodeList[i].begin();
    it_end[i] = --Iter(nodeList[i].end());
  }
  ifs.close();
}

template <class Temp>
void Diagram<Temp> :: check (int site, std::string funcname) const
{
  std::ofstream ofs ("Check.info");
  if (!ofs) {
    std::cout << __FILE__ << " :: " << __func__ << " : cannot open file 'Check.error'\n";
    exit(100);
  }

  // Check end-nodes
  for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi)
    if (!(nodeList[site].begin()->assoc(nbi)->end()) || !(it_end[site]->assoc(nbi)->end())) {
      std::cout << "User check error\n";
      std::cout << "\n In function  " << funcname << "\n\n";
      std::cout << "\n*** End-node associate error ***\n\n";
      ofs << "\n In function  " << funcname << "\n\n";
      ofs << "\n*** End-node associate error ***\n\n";
      ofs << "self info:\n----\n" << nodeList[site].begin()->info(_latt) << "\n";
      ofs.close();
      exit(100);
    }

  for(const_Iter it = ++(nodeList[site].begin()); it != it_end[site]; ++it) {
    // Check conjugate
    if (it->_conji != -1) {
      if (it->conj()->conj() != it) {
        std::cout << "User check error 1\n";
        std::cout << "\n In function  " << funcname << "\n\n";
        std::cout << "\n*** Conjugate error ***\n\n";
        ofs << "\n In function  " << funcname << "\n\n";
        ofs << "\n*** Conjugate error ***\n\n";
        ofs << "self info:\n----\n" << it->info(_latt) << "\n";
        ofs << "conjugate info:\n----\n" << it->conj()->info(_latt) << "\n";
        ofs.close();
        exit(100);
      }
      if (it->conj()->time() != it->time() ) {
        std::cout << "User check error 2\n";
        std::cout << "\n In function  " << funcname << "\n\n";
        std::cout << "\n*** Conjugate time error ***\n\n";
        ofs << "\n In function  " << funcname << "\n\n";
        ofs << "\n*** Conjugate time error ***\n\n";
        ofs << "self info:\n----\n" << it->info(_latt) << "\n";
        ofs << "conjugate info:\n----\n" << it->conj()->info(_latt) << "\n";
        ofs.close();
        exit(100);
      }
    }
    // Check site
    if (it->site() < 0) {
      std::cout << "User check error 9\n";
      std::cout << "\n In function  " << funcname << "\n\n";
      std::cout << "\n*** Site < 0 error ***\n\n";
      ofs << "\n In function  " << funcname << "\n\n";
      ofs << "\n*** Site < 0 error ***\n\n";
      ofs.close();
      exit(100);
    }
    const_Iter itpre = it;
    --itpre;
    if (it->site() != itpre->site() ) {
      std::cout << "User check error 3\n";
      std::cout << "\n In function  " << funcname << "\n\n";
      std::cout << "\n*** Site error ***\n\n";
      ofs << "\n In function  " << funcname << "\n\n";
      ofs << "\n*** Site error ***\n\n";
      ofs.close();
      exit(100);
    }
    // Check time
    if (it->time() < itpre->time()) {
      std::cout << "User check error 4\n";
      std::cout << "\n In function  " << funcname << "\n\n";
      std::cout << "\n*** Time error ***\n\n";
      ofs << "\n In function  " << funcname << "\n\n";
      ofs << "\n*** Time error ***\n\n";
      ofs << "self info:\n----\n" << it->info(_latt) << "\n";
      ofs << "pre-node info:\n----\n" << itpre->info(_latt) << "\n";
      ofs.close();
      exit(100);
    }
    // For associate:
    for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
      Iter assoc = it->assoc(nbi);
      Iter assbef = assoc;
      --assbef;
      int opp = _latt(site).opp(nbi);
      // Check recursive associate
      if (nbi != it->_conji && assoc->assoc(opp) == it) {
        std::cout << "User check error 5\n";
        std::cout << "\n In function  " << funcname << "\n\n";
        std::cout << "\n*** Recursive associate error ***\n\n";
        ofs << "\n In function  " << funcname << "\n\n";
        ofs << "\n*** Recursive associate error ***\n\n";
        ofs << "self info:\n----\n" << it->info(_latt) << "\n";
        ofs << "associate info:\n----\n" << assoc->info(_latt) << "\n";
        ofs.close();
        exit(100);
      }
      // Check associate time
      if (assoc->time() < it->time() ) {
        std::cout << "User check error 6\n";
        std::cout << "\n In function  " << funcname << "\n\n";
        std::cout << "\n*** Assoc time error ***\n\n";
        ofs << "\n In function  " << funcname << "\n\n";
        ofs << "\n*** Assoc time error ***\n\n";
        ofs << "self info:\n----\n" << it->info(_latt) << "\n";
        ofs << "associate info:\n----\n" << assoc->info(_latt) << "\n";
        ofs.close();
        exit(100);
      }
      // Check associate site
      if (assoc->site() != _latt(site).nb(nbi) ) {
        std::cout << "User check error 7\n";
        std::cout << "\n In function  " << funcname << "\n\n";
        std::cout << "\n*** assoc site error ***\n\n";
        ofs << "\n In function  " << funcname << "\n\n";
        ofs << "\n*** assoc site error ***\n\n";
        ofs << "self info:\n----\n" << it->info(_latt) << "\n";
        ofs << "associate-before info:\n----\n" << assbef->info(_latt) << "\n";
        ofs.close();
        exit(100);
      }
      // Check associate-before time
      if (assbef->time() > it->time() ) {
        std::cout << "User check error 8\n";
        std::cout << "\n In function  " << funcname << "\n\n";
        std::cout << "\n*** assoc between error ***\n\n";
        ofs << "\n In function  " << funcname << "\n\n";
        ofs << "\n*** assoc between error ***\n\n";
        ofs << "self info:\n----\n" << it->info(_latt) << "\n";
        ofs << "ass info:\n----\n" << assoc->info(_latt) << "\n";
        ofs << "assbef info:\n----\n" << assbef->info(_latt) << "\n";
        ofs.close();
        exit(100);
      }
    }
  }
  ofs.close();
}
#endif
