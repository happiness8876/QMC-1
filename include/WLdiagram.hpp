#ifndef WLDIAGRAM_HPP_CMC
#define WLDIAGRAM_HPP_CMC
#include <fstream>
#include "Diagram.hpp"
#include "Worm.hpp"

class Operator
{
  public:
    int nbef () const { return _nbef; }
    int naft () const { return _creat ? _nbef + 1 : _nbef - 1; }
    bool creat () const { return _creat; }

  friend std::ostream& operator<<(std::ostream& os, const Operator& op) {
    os << "  nbef = " << op._nbef << "\n";
    os << "  creat = " << op._creat << "\n";
    return os;
  }

  private:
    int _nbef;
    bool _creat;

  friend class WLdiagram;
};

class WLdiagram : public Diagram<Operator>
{
  public:
    typedef Diagram<Operator>::Iter Iter;
    typedef Diagram<Operator>::const_Iter const_Iter;
    typedef GeneralWorm<Operator> Worm;

    WLdiagram (int Lx, int Ly, int Lz, double Beta, bool BCx=0, bool BCy=0, bool BCz=0, int n_init=0);
    virtual void insert (const Iter& it_goal, double time, bool creat, bool up, Worm& mob, Iter& fix);
    //void insert (int site, double time, bool creat, bool up, Worm& mob, Iter& fix);
    void walk (const Worm& worm, double dtime);
    virtual void hop (Worm& worm, int nbi);
    virtual void delete_hop (Worm& worm);
    virtual void relink_hop (Worm& worm, int nbi);
    void cross_end (const Worm& worm);
    void cross_node (const Worm& worm);
    void halt (const Worm& worm);
    void remove (Worm& worm);

    void info_plot (std::string filename, std::string ext="wld") const;
    void info_plot_line (const const_Iter& it, std::ofstream& ofs) const;
    void check (int site, const std::string& funcname="unknown") const;

  protected:
    typedef Diagram<Operator> Diag;
    void set_operator (const Iter& it_bef, const Iter& it_aft, bool creat);
};

WLdiagram :: WLdiagram (int Lx, int Ly, int Lz, double Beta, bool BCx, bool BCy, bool BCz, int n_init)
: Diagram<Operator> (Lx, Ly, Lz, Beta, BCx, BCy, BCz)
{
#ifdef TWO_DOMAIN_INITIAL
  for(int i = 0; i < Diag::SITEN; ++i) {
    int ni = i < (Diag::SITEN/2) ? n_init * ((_latt(i).xi()+_latt(i).yi()) % 2) : n_init * ((_latt(i).xi()+_latt(i).yi()+1) % 2);
    Diag::it_end[i]->z._nbef = ni;
  }
#else
  for(int i = 0; i < Diag::SITEN; ++i) {
    int ni = n_init < 0 ? abs(n_init) * ((_latt(i).xi() + _latt(i).yi()) % 2) : n_init;
    Diag::it_end[i]->z._nbef = ni;
  }
#endif
#ifdef DEBUG_MODE
  for(int site = 0; site < SITEN; ++site)
    check (site, __func__);
#endif
#ifdef DEBUG_TRACK_MODE
__DEBUG_TRACK
#endif
}

void WLdiagram :: insert (const Iter& it_goal, double time, bool creat, bool up, Worm& worm, Iter& fix)
{
  Iter it_bef, it_aft;
  diag_insert (it_goal, time, it_bef, it_aft);
  set_operator (it_bef, it_aft, creat);
  if (up) {
    worm.set (it_aft, up);
    fix = it_bef;
  }
  else {
    worm.set (it_bef, up);
    fix = it_aft;
  }

#ifdef DEBUG_MODE
  int site = it_goal->site();
  //for(site = 0; site < SITEN; site++)
  check (site, __func__);
#endif
#ifdef DEBUG_TRACK_MODE
__DEBUG_TRACK
#endif
}

void WLdiagram :: walk (const Worm& worm, double dtime)
{
  double time = worm.it()->time() + (worm.up() ? dtime : -dtime);
  diag_walk (worm.it(), time, worm.up());
#ifdef DEBUG_TRACK_MODE
__DEBUG_TRACK
#endif
}

void WLdiagram :: hop (Worm& worm, int nbi)
{
  bool up = worm.up();
  Iter it_bef, it_aft;
  diag_hop (worm.it(), nbi, up, it_bef, it_aft);
  if (up) {
    if (worm.it()->z.creat())
      set_operator (it_bef, it_aft, false);
    else
      set_operator (it_bef, it_aft, true);
    worm.it (it_aft);
  }
  else {
    if (worm.it()->z.creat())
      set_operator (it_bef, it_aft, true);
    else
      set_operator (it_bef, it_aft, false);
    worm.it (it_bef);
  }
#ifdef DEBUG_MODE
  int site = it_bef->site();
  //for(site = 0; site < SITEN; site++)
  check (site, __func__);
#endif
#ifdef DEBUG_TRACK_MODE
__DEBUG_TRACK
#endif
}

void WLdiagram :: delete_hop (Worm& worm)
{
#ifdef DEBUG_MODE
int site = worm.it()->site();
#endif
  Iter it_worm = worm.it();
  diag_delete_hop (it_worm, worm.up());
  worm.it (it_worm);
#ifdef DEBUG_MODE
  //for(site = 0; site < SITEN; site++)
  check (site, __func__);
#endif
#ifdef DEBUG_TRACK_MODE
__DEBUG_TRACK
#endif
}

void WLdiagram :: relink_hop (Worm& worm, int nbi)
{
#ifdef DEBUG_MODE
int oldsite = worm.it()->site();
#endif
  Iter it_bef, it_aft, it_worm = worm.it();
  bool up = worm.up();
  diag_relink_hop (it_worm, nbi, up, it_bef, it_aft);
  if (up) {
    if (it_worm->z.creat())
      set_operator (it_bef, it_aft, true);
    else
      set_operator (it_bef, it_aft, false);
    worm.set (it_bef, false);
  }
  else {
    if (it_worm->z.creat())
      set_operator (it_bef, it_aft, false);
    else
      set_operator (it_bef, it_aft, true);
    worm.set (it_aft, true);
  }
#ifdef DEBUG_MODE
  int site = it_bef->site();
  check (oldsite, __func__);
  for(site = 0; site < SITEN; site++)
  check (site, __func__);
#endif
#ifdef DEBUG_TRACK_MODE
__DEBUG_TRACK
#endif
}

void WLdiagram :: cross_end (const Worm& worm)
{
  bool up = worm.up();
  diag_cross_end (worm.it(), up);
  int site = worm.it()->site();
  if (up)
    Diag::it_end[site]->z._nbef = worm.it()->z.nbef();
  else
    Diag::it_end[site]->z._nbef = worm.it()->z.naft();
#ifdef DEBUG_MODE
  //for(site = 0; site < SITEN; site++)
  check (site, __func__);
#endif
#ifdef DEBUG_TRACK_MODE
__DEBUG_TRACK
#endif
}

void WLdiagram :: cross_node (const Worm& worm)
{
  bool up = worm.up();
  diag_cross_node (worm.it(), up);
  if (up) {
    Iter it_bef = worm.it();
    --it_bef;
    int n = worm.it()->z.nbef();
    worm.it()->z._nbef = it_bef->z.nbef();
    it_bef->z._nbef = n;
  }
  else {
    Iter it_aft = worm.it();
    ++it_aft;
    int n = worm.it()->z.nbef();
    worm.it()->z._nbef = it_aft->z.nbef();
    it_aft->z._nbef = n;
  }
#ifdef DEBUG_MODE
  int site = worm.it()->site();
  //for(site = 0; site < SITEN; site++)
  check (site, __func__);
#endif
#ifdef DEBUG_TRACK_MODE
__DEBUG_TRACK
#endif
}

void WLdiagram :: halt (const Worm& worm)
{
  diag_halt (worm.it(), worm.up());
#ifdef DEBUG_TRACK_MODE
__DEBUG_TRACK
#endif
}

void WLdiagram :: remove (Worm& worm)
{
#ifdef DEBUG_MODE
int site = worm.it()->site();
#endif
  Iter it_worm = worm.it();
  diag_remove (it_worm, worm.up());
#ifdef DEBUG_MODE
  //for(site = 0; site < SITEN; site++)
  check (site, __func__);
#endif
#ifdef DEBUG_TRACK_MODE
__DEBUG_TRACK
#endif
}

void WLdiagram :: set_operator (const Iter& it_bef, const Iter& it_aft, bool creat)
{
  Iter aftaft = it_aft;
  ++aftaft;
  int n = aftaft->z.nbef();
  if (creat) {
    it_aft->z._nbef = n + 1;
    it_bef->z._nbef = n;
    it_aft->z._creat = false;
    it_bef->z._creat = true;
  }
  else {
    it_aft->z._nbef = n - 1;
    it_bef->z._nbef = n;
    it_aft->z._creat = true;
    it_bef->z._creat = false;
  }
}

void WLdiagram :: info_plot (std::string filename, std::string ext) const
{
  diag_info_plot (filename, ext);
}

void WLdiagram :: info_plot_line (const const_Iter& it, std::ofstream& ofs) const
{
  ofs << " " << it->time() << " " << it->z.nbef();
}

void WLdiagram :: check (int site, const std::string& funcname) const
{
  std::ofstream ofs ("Check.info");
  if (!ofs) {
    std::cout << __FILE__ << " :: " << __func__ << " : cannot open file 'Check.error'\n";
    exit(100);
  }

  // Check End-nodes n match
  const_Iter it_1st = ++(nodeList[site].begin());
  if (it_1st->z.nbef() != it_end[site]->z.nbef()) {
    std::cout << "User check error\n";
    std::cout << "\n In function  " << funcname << "\n\n";
    std::cout << "\n*** end-node n not match ***\n\n";
    ofs << "\n In function  " << funcname << "\n\n";
    ofs << "\n*** end-node n not match ***\n\n";
    ofs << "first node info:\n----\n" << it_1st->info(_latt) << "\n";
    ofs << "end node info:\n----\n" << it_end[site]->info(_latt) << "\n";
    ofs.close();
    exit(100);
  }

  for(const_Iter it = ++(nodeList[site].begin()); it != it_end[site]; ++it) {
    const_Iter it_aft = ++const_Iter(it);
    // Check creation/annihilation
    if (it != it_end[site]) {
      if (it->z.creat()) {
        if (it->z.nbef() != (it_aft->z.nbef() - 1)) {
          std::cout << "User check error\n";
          std::cout << "\n In function  " << funcname << "\n\n";
          std::cout << "\n*** creation/annihilation error ***\n\n";
          ofs << "\n In function  " << funcname << "\n\n";
          ofs << "\n*** creation/annihilation error ***\n\n";
          ofs << "self info:\n----\n" << it->info(_latt) << "\n";
          ofs << "aft info:\n----\n" << it_aft->info(_latt) << "\n";
          ofs.close();
          exit(100);
        }
      }
      else {
        if (it->z.nbef() != (it_aft->z.nbef() + 1)) {
          std::cout << "User check error\n";
          std::cout << "\n In function  " << funcname << "\n\n";
          std::cout << "\n*** creation/annihilation error ***\n\n";
          ofs << "\n In function  " << funcname << "\n\n";
          ofs << "\n*** creation/annihilation error ***\n\n";
          ofs << "self info:\n----\n" << it->info(_latt) << "\n";
          ofs << "aft info:\n----\n" << it_aft->info(_latt) << "\n";
          ofs.close();
          exit(100);
        }
      }
    }
    // Check n > 0
    if (it->z.nbef() < 0) {
      std::cout << "User check error\n";
      std::cout << "\n In function  " << funcname << "\n\n";
      std::cout << "\n*** n < 0 error ***\n\n";
      ofs << "\n In function  " << funcname << "\n\n";
      ofs << "\n*** n < 0 error ***\n\n";
      ofs << "self info:\n----\n" << it->info(_latt) << "\n";
      ofs.close();
      exit(100);
    }
  }
  ofs.close();
}
#endif
