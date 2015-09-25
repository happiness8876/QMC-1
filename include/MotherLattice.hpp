#ifndef MOTHERLATTICE_HPP_CMC
#define MOTHERLATTICE_HPP_CMC
#include <vector>
#include "Lattice.hpp"
#include "CheckList.hpp"

class MotherLattice : public Lattice
{
  public:
    //MotherLattice (int Lx, int Ly, int Lz, bool bc=0, int xc=-1, int yc=-1, int zc=-1) : Lattice (Lx, Ly, Lz, xc, yc, zc, bc) {}
    MotherLattice (int Lx, int Ly, int Lz, bool bcx=0, bool bcy=0, bool bcz=0, int xc=0, int yc=0, int zc=0) : Lattice (Lx, Ly, Lz, bcx, bcy, bcz, xc, yc, zc) {}
    CheckList sub (int Lx, int Ly, int Lz, int Ox=0, int Oy=0, int Oz=0) const;
    CheckList whole () const;
    //CheckList diff (const CheckList& biglatt, const CheckList& smalllatt);
    std::vector<CheckList> autosqr (int L, int step) const;
    std::vector<CheckList> autoladder (int Lx, int Ly, int stepx, int stepy);
    CheckList out_boundary (const CheckList& sub_) const;
    CheckList in_boundary (const CheckList& sub_) const;
};

CheckList MotherLattice :: out_boundary (const CheckList& sub_) const
{
  CheckList boun (size());
  for(int i = 0; i < sub_.subsize(); ++i) {
    int site = sub_.ele(i);
    for(int nbi = 0; nbi < nbs(site); ++nbi) {
      int nbsite = nb (site, nbi);
      if (!sub_.has (nbsite))
        boun.add (nbsite);
    }
  }
  return boun;
}

CheckList MotherLattice :: in_boundary (const CheckList& sub_) const
{
  CheckList boun (size());
  for(int i = 0; i < sub_.subsize(); ++i) {
    int site = sub_.ele(i);
    for(int nbi = 0; nbi < nbs(site); ++nbi) {
      int nbsite = nb (site, nbi);
      if (!sub_.has (nbsite))
        boun.add (site);
    }
  }
  return boun;
}

std::vector<CheckList> MotherLattice :: autoladder (int Lx_, int Ly_, int xstep, int ystep)
{
  if (Lx_ > Lx() || Ly_ > Ly()) {
    std::cout << "Error: MotherLattice: autoladder: sub-lattice length should be smaller than mother-lattice length\n";
    exit(100);
  }
  std::vector<CheckList> subvec;
  CheckList run (size());
  subvec.push_back (run);
  for(int x = 0; x < Lx_; x += xstep)
    for(int y = 0; y < Ly_; y += ystep) {
      // Include a step-lattice
      for(int dx = 0; dx < xstep; ++dx)
        for(int dy = 0; dy < ystep; ++dy) {
          int xi = x + dx, yi = y + dy;
          run.add (get_index (xi, yi, 0));
        }
      subvec.push_back (run);
    }
  return subvec;
}

std::vector<CheckList> MotherLattice :: autosqr (int L, int stage) const
{
  std::vector<CheckList> subvec;
  CheckList run (size());
  int x = 0, y = 0, r = 1, count = 0;
  bool rotup = true;
  subvec.push_back (run);
  run.add (get_index (x, y, 0));
  subvec.push_back (run);
  for(int i = 0; i < L; ++i) {
    if (rotup) {
      ++x;
      run.add (get_index (x, y, 0));
      if (++count >= stage) {
        subvec.push_back (run);
        count = 0;
      }
      //std::cout << x << ", " << y << "\n";
      for(int j = 0; j < r; ++j) {
        ++y;
        run.add (get_index (x, y, 0));
        if (++count >= stage) {
          subvec.push_back (run);
          count = 0;
        }
        //std::cout << x << ", " << y << "\n";
      }
      for(int j = 0; j < r; ++j) {
        --x;
        run.add (get_index (x, y, 0));
        if (++count >= stage) {
          subvec.push_back (run);
          count = 0;
        }
        //std::cout << x << ", " << y << "\n";
      }
    }
    else {
      ++y;
      run.add (get_index (x, y, 0));
      //std::cout << x << ", " << y << "\n";
      for(int j = 0; j < r; ++j) {
        ++x;
        run.add (get_index (x, y, 0));
        if (++count >= stage) {
          subvec.push_back (run);
          count = 0;
        }
        //std::cout << x << ", " << y << "\n";
      }
      for(int j = 0; j < r; ++j) {
        --y;
        run.add (get_index (x, y, 0));
        if (++count >= stage) {
          subvec.push_back (run);
          count = 0;
        }
        //std::cout << x << ", " << y << "\n";
      }
    }
    if (count != 0) {
      subvec.push_back (run);
      count = 0;
    }
    rotup = !rotup;
    ++r;
  }
  return subvec;
}

CheckList MotherLattice :: sub (int Lx, int Ly, int Lz, int Ox, int Oy, int Oz) const
{
  if (Lx > Lattice::Lx() || Ly > Lattice::Ly() || Lz > Lattice::Lz()) {
    std::cout << "error: " << __FILE__ << " :: " << __func__ << ": sub-lattice size lagger than mother-lattice'\n";
    exit(100);
  }
  CheckList sublattice (size());
  for(int x = 0; x < Lx; ++x)
    for(int y = 0; y < Ly; ++y)
      for(int z = 0; z < Lz; ++z) {
        int index = get_index (Ox+x, Oy+y, Oz+z);
        sublattice.add (index);
      }
  return sublattice;
}

CheckList MotherLattice :: whole () const
{
  CheckList sublattice (size());
  for(int x = 0; x < Lx(); ++x)
    for(int y = 0; y < Ly(); ++y)
      for(int z = 0; z < Lz(); ++z) {
        int index = get_index (x, y, z);
        sublattice.add (index);
      }
  return sublattice;
}
#endif
