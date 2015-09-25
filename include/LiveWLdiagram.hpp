#ifndef LIVEWLDIAGRAM_HPP_CMC
#define LIVEWLDIAGRAM_HPP_CMC
#include <vector>
#include <cmath>
#include "WLdiagram.hpp"
#include "RandGen.hpp"
#include "ProbTable.hpp"
#include "ParaBox.hpp"

class LiveWLdiagram : public WLdiagram
{
  public:
    typedef WLdiagram::Iter Iter;
    LiveWLdiagram (int Lx, int Ly, int Lz, double Beta, double bh_uhalf, double bh_t, double bh_mu, double bh_vnn,
                   double vtrapx=0., double vtrapy=0., double vtrapz=0., double eoffset=0.1, bool BCx=0, bool BCy=0, bool BCz=0, int ninit=0, int nmax=10000,
                   unsigned long rndseed=777777777);
    virtual ~LiveWLdiagram () {}
    virtual void updateZ ();
    double diagonal_energy (const Iter& it, int n) const;
    double onsite_energy (int site, int n) const;
    double longrange_energy (const Iter& it, int n) const;
    void set_mu (double mu);
    virtual void write (std::string name);
    virtual void read (std::string name);
    void save (const std::string& name) { write (name); rand.save (name); }
    void load (const std::string& name) { read (name); rand.load (name); }

  protected:
    typedef WLdiagram Wld;
    typedef WLdiagram::Worm Worm;
    enum Status { UNKNOWN, FREE, HALTED, CRASH, PASS_NEIGHBOR, CROSS_NODE, CROSS_END };

    Worm mob;
    Iter fix;
    RandGen rand;
    ProbTable table;
    // Parameters
    int nMAX;
    static const int nMIN = 0;
    double BH_UHALF, BH_t, BH_MU, BH_Vnn, VTRAPX, VTRAPY, VTRAPZ, EOFFSET;
    std::vector<double> EFF_MU;

    virtual Status updateZe (Worm& worm);
    virtual Status walk (double dtime, const Worm& worm, Iter& it_block);
    void act (Worm& worm, Status stat, double Ewalk, const Iter& it_block);
    double walk_energy (const Iter& it, bool up);
    double shift_energy (double E, double Esh);
    void shifted_diagonal_energy (const Iter& it, double& Ebef, double& Eaft);
    bool check_nlimit (int n, bool creat, int nmax, int nmin);
};

LiveWLdiagram :: LiveWLdiagram(int Lx, int Ly, int Lz, double Beta, double bh_uhalf, double bh_t, double bh_mu, double bh_vnn, double vtrapx, double vtrapy, double vtrapz, double eoffset, bool BCx, bool BCy, bool BCz, int n_init, int nmax, unsigned long rndseed)
: WLdiagram (Lx, Ly, Lz, Beta, BCx, BCy, BCz, n_init), BH_UHALF (bh_uhalf), BH_t (bh_t), BH_MU (bh_mu), BH_Vnn(bh_vnn), VTRAPX (vtrapx)
, VTRAPY (vtrapy), VTRAPZ (vtrapz), EOFFSET (eoffset), nMAX (nmax), rand (rndseed)
{
  // Calculate effactive mu
  EFF_MU.resize (Wld::Diag::SITEN);
  for(int i = 0; i < Wld::Diag::SITEN; ++i) {
    int x = Wld::Diag::_latt(i).xi();
    int y = Wld::Diag::_latt(i).yi();
    int z = Wld::Diag::_latt(i).zi();
    double mu_shift = VTRAPX * x * x + VTRAPY * y * y + VTRAPZ * z * z;
    EFF_MU[i] = BH_MU - mu_shift;
    #ifdef NN_INTERACTION
    #ifdef SPIN_HALF_XXZ_MODEL
    EFF_MU[i] += 0.5*BH_Vnn*__latt(i).nbs();
    #endif
    #endif
  }
}

void LiveWLdiagram :: set_mu (double mu)
{
  BH_MU = mu;
  for(int i = 0; i < SITEN; ++i) {
    int x = _latt(i).xi();
    int y = _latt(i).yi();
    int z = _latt(i).zi();
    double mu_shift = VTRAPX * x * x + VTRAPY * y * y + VTRAPZ * z * z;
    EFF_MU[i] = BH_MU - mu_shift;
    #ifdef NN_INTERACTION
    #ifdef SPIN_HALF_XXZ_MODEL
    EFF_MU[i] += 0.5*BH_Vnn*__latt(i).nbs();
    #endif
    #endif
  }
}

void LiveWLdiagram :: updateZ ()
{
  // Insert worms to the diagram
  int site = rand.choice (Wld::Diag::SITEN);
  double time = Wld::Diag::BETA * rand.uniform();
  bool up = rand.choice2();
  bool creat = rand.choice2();
  Iter it_goal = Wld::Diag::find_aft (site, time);
  int ninit = it_goal->z.nbef();
  // Check whether cannot creat or annihilate particle anymore
  if (!check_nlimit (ninit, creat, nMAX, nMIN)) return;
  insert (it_goal, time, creat, up, mob, fix);

  // Worm runs until crash
  while ((updateZe (mob)) != CRASH);
}

LiveWLdiagram::Status LiveWLdiagram :: updateZe (Worm& worm)
{
  double Ewalk = walk_energy (worm.it(), !worm.up());
  double dtime = -log (rand.uniform()) / Ewalk;
  double time_pre = worm.it()->time();

  Iter it_block;
  Status stat = walk (dtime, worm, it_block);
  act (worm, stat, Ewalk, it_block);
  return stat;
}

void LiveWLdiagram :: act (Worm& worm, Status stat, double Ewalk, const Iter& it_block)
{
  table.clear();
  int site = worm.it()->site();

  if (stat == FREE) {
    // 1.1 bounce
    table.add (Ewalk);
    // 1.2 insert interaction
    for(int nbi = 0; nbi < Wld::Diag::_latt(site).nbs(); ++nbi) {
      int n = worm.it()->assoc(nbi)->z.nbef();
      bool creat = (worm.up() != worm.it()->z.creat());
      if (check_nlimit (n, creat, nMAX, nMIN)) {
        // Check whether cannot creat or annihilate particle anymore
        double weight = BH_t * (creat ? n+1 : n);
        table.add (weight, nbi);
      }
    }
    // 1.3 choose
   #ifdef LOCAL_OPTIMAL
    table.locally_optimal();
   #endif
    int choice = table.get_choice (rand.uniform());
    if (choice == 0) worm.bounce();
    else
      hop (worm, table.nbi(choice));
  }
  else if (stat == HALTED) {
    // 2.1 bounce
    double weight = BH_t * (worm.it()->z.creat() ? worm.it()->z.naft() : worm.it()->z.nbef());
    table.add (weight);
    // 2.2 remove interaction
    Iter it_to = it_block->conj(); // 'it_to' is the to-node in remove case, or the middle-node in relink case

    if (mob.up()) {
      int site = it_to->site();
      double Ebef = onsite_energy (site, it_to->z.nbef());
      double Eaft = onsite_energy (site, it_to->z.naft());
#ifdef NN_INTERACTION
      // long range energy
      int nbns = 0;
      for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi)
        nbns += it_to->ass_aft(nbi)->z.nbef();
      Ebef += BH_Vnn * (it_to->z.nbef() * nbns);
      Eaft += BH_Vnn * (it_to->z.naft() * nbns);
#endif
      weight = shift_energy (Eaft, Ebef);
    }
    else {
      int site = it_to->site();
      double Ebef = onsite_energy (site, it_to->z.nbef());
      double Eaft = onsite_energy (site, it_to->z.naft());
#ifdef NN_INTERACTION
      // long range energy
      int nbns = 0;
      for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi)
        nbns += it_to->assoc(nbi)->z.nbef();
      Ebef += BH_Vnn * (it_to->z.nbef() * nbns);
      Eaft += BH_Vnn * (it_to->z.naft() * nbns);
#endif
      weight = shift_energy (Ebef, Eaft);
    }
    table.add (weight);
    // 2.3 relink interaction
    for(int nbi = 0; nbi < Wld::Diag::_latt(it_to->site()).nbs(); ++nbi) {
      Iter it_nb = it_to->assoc(nbi);
      if (it_nb->site() != worm.it()->site()) {
        // Check whether relink to the self site
        int n = it_nb->z.nbef();
        bool creat = (worm.up() == worm.it()->z.creat());
        if (check_nlimit (n, creat, nMAX, nMIN)) {
          // Check whether cannot creat or annihilate particle anymore
          weight = BH_t * (creat ? n+1 : n);
          table.add (weight, nbi);
        }
      }
    }
    // 2.4 choose
#ifdef LOCAL_OPTIMAL
    table.locally_optimal();
#endif
    int choice = table.get_choice (rand.uniform());
    if (choice == 0) worm.bounce();
    else if (choice == 1)
      delete_hop (worm);
    else
      relink_hop (worm, table.nbi(choice));
  }

  else if (stat == CRASH)
    Wld::remove (worm);
}

LiveWLdiagram::Status LiveWLdiagram :: walk (double dtime, const Worm& worm, Iter& it_block)
{
  // Get the block and the time-distance to the block
  int site = worm.it()->site();
  it_block = worm.up() ? ++Iter(worm.it()) : --Iter(worm.it());
#ifdef NN_INTERACTION
  // Check the neighbor kinks
  int block_nbi = -1;
  if (worm.up())
    for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
      Iter it_ass = worm.it()->assoc(nbi);
      if (it_ass->time() < it_block->time()) {
        it_block = it_ass;
        block_nbi = nbi;
      }
    }
  else
    for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi) {
      Iter it_ass = --Iter(worm.it()->assoc(nbi));
      if (it_ass->time() > it_block->time()) {
        it_block = it_ass;
        block_nbi = nbi;
      }
    }
#endif
  double dtime_block = fabs (it_block->time() - worm.it()->time());

  // Free
  if (dtime < dtime_block) {
    Wld::walk (worm, dtime);
    return FREE;
  }
#ifdef NN_INTERACTION
  // Pass neighbor
  else if (block_nbi != -1) {
    pass_neighbor (worm.up(), worm.it(), it_block, block_nbi);
    return PASS_NEIGHBOR;
  }
#endif
  // Cross end
  else if (it_block->end()) {
    Wld::cross_end (worm);
    return CROSS_END;
  }
  // Cross node
  else if (it_block->z.creat() == worm.it()->z.creat()) {
    Wld::cross_node (worm);
    return CROSS_NODE;
  }
  // Meet fix-worm
  else if (it_block == fix) {
    return CRASH;
  }
  // Halted
  else {
    Wld::halt (worm);
    return HALTED;
  }
}

inline double LiveWLdiagram :: walk_energy (const Iter& it, bool up)
{
  double Ebef, Eaft;
  shifted_diagonal_energy (it, Ebef, Eaft);
  return up ? Eaft : Ebef;
}

void LiveWLdiagram :: shifted_diagonal_energy (const Iter& it, double& Ebef, double& Eaft)
{
  // On-site energy
  int nbef = it->z.nbef(), naft = it->z.naft(), site = it->site();
  int dn = naft - nbef;
  double dE = -EFF_MU[site] * dn;
  dE += BH_UHALF * (naft*(naft-1) - nbef*(nbef-1));
#ifdef NN_INTERACTION
    int neighbor_n = 0;
    for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi)
      neighbor_n += it->assoc(nbi)->z.nbef();
    dE += BH_Vnn * (dn * neighbor_n);
#endif
  if (dE > 0) {
    Eaft = dE + EOFFSET;
    Ebef = EOFFSET;
  }
  else if (dE < 0) {
    Eaft = EOFFSET;
    Ebef = -dE + EOFFSET;
  }
  else {
    Eaft = EOFFSET;
    Ebef = EOFFSET;
  }
}

inline double LiveWLdiagram :: diagonal_energy (const Iter& it, int n) const
{
  double Ed = onsite_energy (it->site(), n);
#ifdef NN_INTERACTION
  Ed += longrange_energy (it, n);
#endif
  return Ed;
}

double LiveWLdiagram :: longrange_energy (const Iter& it, int n) const
{
  int site = it->site();
  // Nearest neighbors interaction
  int neighbor_n = 0;
  for(int nbi = 0; nbi < _latt(site).nbs(); ++nbi)
    neighbor_n += it->assoc(nbi)->z.nbef();
  return BH_Vnn * (n * neighbor_n);
}

inline double LiveWLdiagram :: onsite_energy (int site, int n) const
{
  return (BH_UHALF * (n-1) - EFF_MU[site]) * n;
}

inline double LiveWLdiagram :: shift_energy (double E, double Esh)
{
  return E <= Esh ? EOFFSET : E - Esh + EOFFSET;
}

inline bool LiveWLdiagram :: check_nlimit (int n, bool creat, int nmax, int nmin)
{
  if (n == nmax && creat) return false;
  else if (n == nmin && !creat) return false;
  return true;
}

void LiveWLdiagram :: write (std::string name)
{
  Wld::write (name);
  std::string fullname = name+".livewld";
  std::ofstream ofs (fullname.c_str(), std::ios::binary);
  if (!ofs) {
    std::cout << "Error open file: LiveWLdiagram: write: " << fullname << "\n";
    exit(100);
  }
  // Parameters
  ofs.write ((char*)&nMAX, sizeof(nMAX));
  ofs.write ((char*)&BH_UHALF, sizeof(BH_UHALF));
  ofs.write ((char*)&BH_t, sizeof(BH_t));
  ofs.write ((char*)&BH_MU, sizeof(BH_MU));
  ofs.write ((char*)&BH_Vnn, sizeof(BH_Vnn));
  ofs.write ((char*)&VTRAPX, sizeof(VTRAPX));
  ofs.write ((char*)&VTRAPY, sizeof(VTRAPY));
  ofs.write ((char*)&VTRAPZ, sizeof(VTRAPZ));
  ofs.write ((char*)&EOFFSET, sizeof(EOFFSET));
  ofs.write ((char*)&BETA, sizeof(BETA));
  ofs.write ((char*)&EFF_MU[0], sizeof(EFF_MU[0])*EFF_MU.size());
  ofs.close();
}

void LiveWLdiagram :: read (std::string name)
{
  Wld::read (name);
  std::string fullname = name+".livewld";
  std::ifstream ifs (fullname.c_str(), std::ios::binary);
  if (!ifs) {
    std::cout << "Error open file: LiveWLdiagram: read: " << fullname << "\n";
    exit(100);
  }
  // Parameters
  ifs.read ((char*)&nMAX, sizeof(nMAX));
  ifs.read ((char*)&BH_UHALF, sizeof(BH_UHALF));
  ifs.read ((char*)&BH_t, sizeof(BH_t));
  ifs.read ((char*)&BH_MU, sizeof(BH_MU));
  ifs.read ((char*)&BH_Vnn, sizeof(BH_Vnn));
  ifs.read ((char*)&VTRAPX, sizeof(VTRAPX));
  ifs.read ((char*)&VTRAPY, sizeof(VTRAPY));
  ifs.read ((char*)&VTRAPZ, sizeof(VTRAPZ));
  ifs.read ((char*)&EOFFSET, sizeof(EOFFSET));
  ifs.read ((char*)&BETA, sizeof(BETA));
  ifs.read ((char*)&EFF_MU[0], sizeof(EFF_MU[0])*EFF_MU.size());
  ifs.close();
}
#endif
