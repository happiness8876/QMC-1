#ifndef INFOWLDIAGRAM_HPP_CMC
#define INFOWLDIAGRAM_HPP_CMC
#include <map>
#include "LiveWLdiagram.hpp"
#include "CheckList.hpp"
#include "ParaBox.hpp"
#include "ParaBox.hpp"

class InfoWLdiagram : public LiveWLdiagram
{
  public:
    typedef LiveWLdiagram LWld;
    InfoWLdiagram (int Lx, int Ly, int Lz, double Beta, double bh_uhalf, double bh_t, double bh_mu, double bh_vnn,
                   double vtrapx, double vtrapy, double vtrapz, double eoffset, bool BCx=0, bool BCy=0, bool BCz=0, int ninit=0, int nmax=10000, unsigned long rndseed=777777777)
      : LiveWLdiagram (Lx, Ly, Lz, Beta, bh_uhalf, bh_t, bh_mu, bh_vnn, vtrapx, vtrapy, vtrapz, eoffset, BCx, BCy, BCz, ninit, nmax,
                   rndseed), _updated_sites (Lx*Ly*Lz) { init_MS(); }
    InfoWLdiagram (const ParaBox& para)
      : LiveWLdiagram (para["LX"], para["LY"], para["LZ"],1./double(para["TEMPERATURE"]), 0.5*(double)para["BH_U"], para["BH_t"], para["BH_MU"], para["BH_Vnn"], para["VX_TRAP"], para["VY_TRAP"], para["VZ_TRAP"], para["EOFFSET"], para["BOUNDARY_CONDITION_X"], para["BOUNDARY_CONDITION_Y"], para["BOUNDARY_CONDITION_Z"], para["INIT_n"], para["nMAX"], para["RandSeed"]), _updated_sites (_latt.size()) { init_MS(); }

    // Overload function
    Status updateZe (Worm& worm);
    void updateZ ();
    void read (std::string name) { LWld::read(name); init_MS(); }

    // Member
    int size () const { return _latt.size(); }
    double beta () const { return BETA; }
    const MotherLattice& latt () const { return _latt; }
    const CheckList& updated_sites () const { return _updated_sites; }
    unsigned long seed () const { return rand.seed(); }

    // Information
    int n (int i) const { return it_end[i]->z.nbef(); }
    double nBETA (int i) const;
    double nnBETA (int i, int j) const;
    double eff_mu (int i) const { return EFF_MU[i]; }
    int N () const;
    int NA () const;
    double Ed () const;
    int Nhop () const;
    std::vector<int> W ();
    double end_time () const { return it_end[0]->time(); }
    const_Iter it_beg (int i) const { return nodeList[i].begin(); }
#ifdef TRACING_WINDING_NUMBER
    int W (int i) const { return MS_W[i]; }
    int W2 () const { return MS_W[0]*MS_W[0] + MS_W[1]*MS_W[1] + MS_W[2]*MS_W[2]; }
#endif
    std::map<int,double> NA_weight () const; // Only for 1D system
#ifdef MEASURE_CORRELATION
    const int C_bb (int i) const { return (size()%2==0 && i==halfL) ? MS_bb(i)*extended_shift*2 : MS_bb(i)*extended_shift; }
    int Ze_shift_weight () const { return extended_shift; }
#endif
    void set_subsystem_A (const CheckList& A);
    void info_plot (std::string filename) const { WLdiagram::info_plot (filename); }

  private:
    int MS_N, MS_Nhop, MS_NA;
    std::vector<int> MS_W;
    double MS_Ediag;
    std::vector<double> MS_Ebond, MS_Eon;
    CheckList _updated_sites, _subA, _Aboundary;
#ifdef MEASURE_CORRELATION
    MyVector<int> MS_bb;
    int extended_shift;
#endif

    void init_MS ();
    void hop (Worm& worm, int nbi);
    void delete_hop (Worm& worm);
    void relink_hop (Worm& worm, int nbi);
    std::pair<int,int> winding_element (int site1, bool creat1, int site2) const;
    void add_winding_num (const Iter& it1, const Iter& it2, std::vector<int>& ms_w);
    void del_winding_num (const Iter& it1, const Iter& it2, std::vector<int>& ms_w);
#ifdef MEASURE_CORRELATION
    void insert (const Iter& it_goal, double time, bool creat, bool up, Worm& worm, Iter& fix);
    LiveWLdiagram::Status walk (double dtime, const Worm& worm, Iter& it_block);
#endif
};

class Time_slice
{
  public:
    Time_slice (const InfoWLdiagram& wld, const std::vector<int>& sites) { init(wld,sites); }
    Time_slice (const InfoWLdiagram& wld) { std::vector<int> s; for(int i = 0; i < wld.latt().size(); ++i) s.push_back(i); init(wld,s); }
    void init (const InfoWLdiagram& wld, const std::vector<int>& sites);
    const InfoWLdiagram::const_Iter& it (int i) const { return its[i]; }
    const InfoWLdiagram::const_Iter& first () const { return its_sort.front(); }
    double dt () const { return _dt; }
    bool next ();
  private:
    double END_TIME, _dt, tpre;
    std::vector<InfoWLdiagram::const_Iter> its; // Store the iterators by site
    std::vector<InfoWLdiagram::const_Iter> its_sort; // Store the iterators sorted by time
    static bool compare_time (InfoWLdiagram::const_Iter it1, InfoWLdiagram::const_Iter it2) { return (it1->time() < it2->time()); }
    void sort_first (std::vector<InfoWLdiagram::const_Iter>& a);
};

void InfoWLdiagram :: set_subsystem_A (const CheckList& A)
{
  _subA = A;
  _Aboundary = latt().in_boundary(A); 
  // Initialize MS_NA
  MS_NA = 0;
  for(int i = 0; i < _subA.subsize(); ++i)
    MS_NA += it_end[_subA.ele(i)]->z.nbef();
}

#ifdef MEASURE_CORRELATION
void InfoWLdiagram :: insert (const Iter& it_goal, double time, bool creat, bool up, Worm& worm, Iter& fix)
{
  WLdiagram::insert (it_goal, time, creat, up, worm, fix);
  extended_shift = worm.it()->z.creat() ? fix->z.nbef() : worm.it()->z.nbef();
  //if (!creat)
    ++(MS_bb(0));
}

LiveWLdiagram::Status InfoWLdiagram :: walk (double dtime, const Worm& worm, Iter& it_block)
{
  double tbef = worm.it()->time();
  Status stat = LiveWLdiagram::walk (dtime, worm, it_block);
  double taft = worm.it()->time();
  if (stat == CROSS_END) {
    if (worm.up() == (tbef < fix->time())) {
  #ifdef ONE_DIMENSION
      int d = PBC_abs_dis (worm.it()->site(), fix->site(), size());
      ++MS_bb (d);
  #else
      int dx = PBC_abs_dis (_latt(worm.it()->site()).xi(), _latt(fix->site()).xi(), _latt.Lx());
      int dy = PBC_abs_dis (_latt(worm.it()->site()).yi(), _latt(fix->site()).yi(), _latt.Ly());
      int dz = PBC_abs_dis (_latt(worm.it()->site()).zi(), _latt(fix->site()).zi(), _latt.Lz());
      ++MS_bb (dx,dy,dz);
  #endif
    }
  }
  else if (stat == CRASH) {
      ++MS_bb(0);
  }
  else {
    if ((tbef < fix->time() && taft >= fix->time()) || (tbef > fix->time() && taft <= fix->time())) {
  #ifdef ONE_DIMENSION
      int d = PBC_abs_dis (worm.it()->site(), fix->site(), size());
      ++MS_bb (d);
  #else
      int dx = PBC_abs_dis (_latt(worm.it()->site()).xi(), _latt(fix->site()).xi(), _latt.Lx());
      int dy = PBC_abs_dis (_latt(worm.it()->site()).yi(), _latt(fix->site()).yi(), _latt.Ly());
      int dz = PBC_abs_dis (_latt(worm.it()->site()).zi(), _latt(fix->site()).zi(), _latt.Lz());
      ++MS_bb (dx,dy,dz);
  #endif
    }
  }
  return stat;
}
#endif

std::map<int,double> InfoWLdiagram :: NA_weight () const
{
  std::map<int,double> NA_w;
  if (_subA.size() == 0) return NA_w;
  // Count the particles at t = 0
  int NA = 0;
  for(int i = 0; i < _subA.subsize(); ++i)
    NA += it_end[_subA.ele(i)]->z.nbef();
  Time_slice boundary_state (*this, _Aboundary.contains());
  double tpre = 0.;
  do {
    const_Iter it = boundary_state.first();
    if (it->end())
      NA_w[NA] += it->time() - tpre;
    else if (!_subA.has(it->conj()->site())) {
      NA_w[NA] += it->time() - tpre;
      tpre = it->time();
      NA += it->z.creat() ? 1 : -1;
    }
  }while (boundary_state.next());
  return NA_w;
}

std::vector<int> InfoWLdiagram :: W ()
{
#ifdef TRACING_WINDING_NUMBER
  return MS_W;
#endif
  std::vector<int> tempW (3, 0);
  for(int i = 0; i < SITEN; ++i)
    for(Iter it = ++nodeList[i].begin(); it != it_end[i]; ++it)
      add_winding_num (it, it->conj(), tempW);
  for(int i = 0; i < 3; i++) {
    if (tempW[i] % 2 != 0) {
      std::cout << "InfoWLdiagram :: init_MS: Winding number error\n" << i << ": " << tempW[i] << "\n";
      exit(100);
    }
    tempW[i] /= 2;
  }
  return tempW;
}

int InfoWLdiagram :: Nhop () const
{
#ifdef TRACING_N_HOP
  return MS_Nhop;
#endif
  int Nh = 0;
  for(int i = 0; i < SITEN; ++i)
    Nh += nodeList[i].size();
  Nh -= 2*SITEN;
  if (Nh % 2 != 0) { std::cout << "***Nhop error\n" << Nh << "\n"; exit(100); }
  return Nh / 2;
}

double InfoWLdiagram :: Ed () const
{
#ifdef TRACING_DIAG_E
  return MS_Ediag;
#endif
  double Eon = 0., Enb = 0.;
  for(int i = 0; i < SITEN; ++i) {
    Eon += LWld::onsite_energy (i, it_end[i]->z.nbef());
#ifdef NN_INTERACTION
    Enb += LWld::longrange_energy (it_end[i], it_end[i]->z.nbef());
#endif
  }
  return Eon + 0.5*Enb;
}

int InfoWLdiagram :: N () const
{
#ifdef TRACING_N
  return MS_N;
#endif
  int N = 0;
  for(int i = 0; i < SITEN; ++i)
    N += it_end[i]->z.nbef();
  return N;
}

int InfoWLdiagram :: NA () const
{
  #ifdef TRACING_NA
  return MS_NA;
  #endif
  int N = 0;
  for(int i = 0; i < _subA.subsize(); ++i) {
    int site = _subA.ele(i);
    N += it_end[site]->z.nbef();
  }
  return N;
}

double InfoWLdiagram :: nBETA (int i) const
{
  double nbeta = 0., tpre = 0.;
  for(const_Iter it = ++(nodeList[i].begin()); it != nodeList[i].end(); ++it) {
    double dt = it->time() - tpre;
    nbeta += dt * it->z.nbef();
    tpre = it->time();
  }
  return nbeta;
}

double InfoWLdiagram :: nnBETA (int i, int j) const
{
  const_Iter iti = ++(nodeList[i].begin()), itj = ++(nodeList[j].begin());
  double nnbeta = 0., tpre = 0.;
  do {
    double nn = iti->z.nbef() * itj->z.nbef();
    if (iti->time() < itj->time()) {
      nnbeta += nn * (iti->time() - tpre);
      tpre = iti->time();
      ++iti;
    }
    else if (iti->time() > itj->time()) {
      nnbeta += nn * (itj->time() - tpre);
      tpre = itj->time();
      ++itj;
    }
    else { // ti == tj
      nnbeta += nn * (iti->time() - tpre);
      tpre = iti->time();
      ++iti; ++itj;
    }
  }while (iti != nodeList[i].end() || itj != nodeList[j].end());
  return nnbeta;
}

void InfoWLdiagram :: init_MS ()
{
  // Initialize MS_N
  MS_N = 0;
  for(int i = 0; i < SITEN; ++i)
    MS_N += it_end[i]->z.nbef();
  // Initialize energy: Eon[], Ebond[], Ediag
  MS_Ediag = 0.;
  MS_Eon.resize (SITEN);
  MS_Ebond.resize (_latt.bonds());
  std::vector<bool> unbond (_latt.bonds(), true);
  for(int i = 0; i < SITEN; i++) {
    MS_Eon[i] = onsite_energy (i, it_end[i]->z.nbef());
    MS_Ediag += MS_Eon[i];
#ifdef NN_INTERACTION
    for(int nbi = 0; nbi < _latt(i).nbs(); ++nbi) {
      int nb = _latt(i).nb(nbi);
      int bi = _latt(i).bi(nbi);
      if (unbond[bi]) {
        MS_Ebond[bi] = BH_Vnn * it_end[i]->z.nbef() * it_end[nb]->z.nbef();
        MS_Ediag += MS_Ebond[bi];
        unbond[bi] = false;
      }
    }
#endif
  }
  // Initialize energy: MS_Nhop
  int Nh = 0;
  for(int i = 0; i < SITEN; ++i)
    Nh += nodeList[i].size();
  Nh -= 2*SITEN;
  if (Nh % 2 != 0) { std::cout << "InfoWLdiagram :: init_MS: Nhop error\n" << Nh << "\n"; exit(100); }
  MS_Nhop = Nh / 2;
  // Initialize winding number: MS_W
  MS_W.resize (3, 0);
  for(int i = 0; i < SITEN; ++i)
    for(Iter it = ++nodeList[i].begin(); it != it_end[i]; ++it)
      add_winding_num (it, it->conj(), MS_W);
  for(int i = 0; i < 3; i++) {
    if (MS_W[i] % 2 != 0) {
      std::cout << "InfoWLdiagram :: init_MS: Winding number error\n" << i << ": " << MS_W[i] << "\n";
      exit(100);
    }
    MS_W[i] /= 2;
  }
}

void InfoWLdiagram :: hop (Worm& worm, int nbi)
{
  WLdiagram::hop (worm, nbi);
#ifdef TRACING_N_HOP
  ++MS_Nhop;
#endif
#ifdef TRACING_WINDING_NUMBER
  Iter it = worm.up() ? --Iter(worm.it()) : ++Iter(worm.it());
  add_winding_num (it, it->conj(), MS_W);
#endif
}

void InfoWLdiagram :: delete_hop (Worm& worm)
{
#ifdef TRACING_WINDING_NUMBER
  Iter it = worm.up() ? ++Iter(worm.it()) : --Iter(worm.it());
  del_winding_num (it, it->conj(), MS_W);
#endif
  WLdiagram::delete_hop (worm);
#ifdef TRACING_N_HOP
  --MS_Nhop;
#endif
}

void InfoWLdiagram :: relink_hop (Worm& worm, int nbi)
{
#ifdef TRACING_WINDING_NUMBER
  Iter it = worm.up() ? ++Iter(worm.it()) : --Iter(worm.it());
  del_winding_num (it, it->conj(), MS_W);
#endif
  WLdiagram::relink_hop (worm, nbi);
#ifdef TRACING_WINDING_NUMBER
  it = worm.up() ? --Iter(worm.it()) : ++Iter(worm.it());
  add_winding_num (it, it->conj(), MS_W);
#endif
}

void InfoWLdiagram :: updateZ ()
{
#ifdef MEASURE_CORRELATION
  MS_bb.clear();
#endif
  _updated_sites.clear();
  LWld::updateZ();
}

InfoWLdiagram::LiveWLdiagram::Status InfoWLdiagram :: updateZe (Worm& worm)
{
  LWld::Status stat = LWld::updateZe (worm);
  int site = worm.it()->site();
  if (stat != CRASH) _updated_sites.add (site);

  // Update MS_N
  if (stat == CROSS_END) {
#ifdef TRACING_N
    MS_N += (worm.it()->z.creat() == worm.up() ? -1 : 1);
#endif
#ifdef TRACING_DIAG_E
    // Update MS_Eon[], MS_Ebond[] and MS_Ediag
    MS_Ediag -= MS_Eon[site];
    MS_Eon[site] = onsite_energy (site, it_end[site]->z.nbef());
    MS_Ediag += MS_Eon[site];
  #ifdef NN_INTERACTION
    for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
      int nb = _latt(site).nb(nbi);
      int bi = _latt(site).bi(nbi);
      MS_Ediag -= MS_Ebond[bi];
      MS_Ebond[bi] = BH_Vnn * it_end[site]->z.nbef() * it_end[nb]->z.nbef();
      MS_Ediag += MS_Ebond[bi];
    }
  #endif
#endif
  }
#ifdef TRACING_NA
  // Update MS_NA
  if (stat == CROSS_END && _subA.has(site))
    MS_NA += (worm.it()->z.creat() == worm.up() ? -1 : 1);
#endif
  return stat;
}

void InfoWLdiagram :: add_winding_num (const Iter& it1, const Iter& it2, std::vector<int>& ms_w)
{
  int sitea = it1->site(), siteb = it2->site();
  if (_latt(sitea).absxi() == 0 && _latt(siteb).absxi() == 1)
    ms_w[0] += it1->z.creat() ? -1 : 1;
  else if (_latt(sitea).absxi() == 1 && _latt(siteb).absxi() == 0)
    ms_w[0] += it2->z.creat() ? -1 : 1;
  else if (_latt(sitea).absyi() == 0 && _latt(siteb).absyi() == 1)
    ms_w[1] += it1->z.creat() ? -1 : 1;
  else if (_latt(sitea).absyi() == 1 && _latt(siteb).absyi() == 0)
    ms_w[1] += it2->z.creat() ? -1 : 1;
  else if (_latt(sitea).abszi() == 0 && _latt(siteb).abszi() == 1)
    ms_w[2] += it1->z.creat() ? -1 : 1;
  else if (_latt(sitea).abszi() == 1 && _latt(siteb).abszi() == 0)
    ms_w[2] += it2->z.creat() ? -1 : 1;
}

void InfoWLdiagram :: del_winding_num (const Iter& it1, const Iter& it2, std::vector<int>& ms_w)
{
  int sitea = it1->site(), siteb = it2->site();
  if (_latt(sitea).absxi() == 0 && _latt(siteb).absxi() == 1)
    ms_w[0] -= it1->z.creat() ? -1 : 1;
  else if (_latt(sitea).absxi() == 1 && _latt(siteb).absxi() == 0)
    ms_w[0] -= it2->z.creat() ? -1 : 1;
  else if (_latt(sitea).absyi() == 0 && _latt(siteb).absyi() == 1)
    ms_w[1] -= it1->z.creat() ? -1 : 1;
  else if (_latt(sitea).absyi() == 1 && _latt(siteb).absyi() == 0)
    ms_w[1] -= it2->z.creat() ? -1 : 1;
  else if (_latt(sitea).abszi() == 0 && _latt(siteb).abszi() == 1)
    ms_w[2] -= it1->z.creat() ? -1 : 1;
  else if (_latt(sitea).abszi() == 1 && _latt(siteb).abszi() == 0)
    ms_w[2] -= it2->z.creat() ? -1 : 1;
}

void Time_slice :: init (const InfoWLdiagram& wld, const std::vector<int>& sites)
{
  END_TIME = wld.end_time();
  int siteN = wld.latt().size();
  its.resize (wld.latt().size(), InfoWLdiagram::NULL_ITER());
  its_sort.clear();
  for(int i = 0; i < sites.size(); ++i) {
    InfoWLdiagram::const_Iter it = ++wld.it_beg(sites[i]);
    its[sites[i]] = it; // Iterators for a state of specific time
    its_sort.push_back (it); // Iterators sorted in time
  }
  std::sort (its_sort.begin(), its_sort.end(), compare_time);
  _dt = its_sort.front()->time();
  tpre = _dt;
}

void Time_slice :: sort_first (std::vector<InfoWLdiagram::const_Iter>& a)
{
  InfoWLdiagram::const_Iter it_first = a.front();
  double t = a.front()->time();
  std::vector<InfoWLdiagram::const_Iter>::iterator it = a.begin(), it_run = ++a.begin();
  while (it_run != a.end() && t > (*it_run)->time()) {
    *it = *it_run;
    ++it_run;
    ++it;
  }
  *it = it_first;
}

bool Time_slice :: next ()
{
  // ++ the iterators with smallest time
  double time = its_sort.front()->time();
  if (time == END_TIME) return false;
  std::list<InfoWLdiagram::const_Iter> temp;
  do {
    int site = its_sort.front()->site();
    ++(its[site]);
    ++(*(its_sort.begin()));
    sort_first (its_sort); // Sort 'its_sort'
  }while (its_sort.front()->time() == time); // loop is for iterators with the same time 
  // Calculate the time-distance for this state
  _dt = its_sort.front()->time() - tpre;
  tpre = its_sort.front()->time();
  return true;
}
#endif
