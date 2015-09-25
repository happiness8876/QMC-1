#ifndef LATTICE_HPP
#define LATTICE_HPP
#include <cstdlib>
#include <vector>
#include <fstream>
#include <cmath>
#include "SingleSite.hpp"

class Lattice
{
  public:
    //Lattice (int xsize, int ysize, int zsize, BoundCond bc=0, int xc=-1, int yc=-1, int zc=-1);
    Lattice (int xsize, int ysize, int zsize, bool bcx=0, bool bcy=0, bool bcz=0, int xc=-1, int yc=-1, int zc=-1);
    virtual ~Lattice () {}

    int Lx () const { return _xsize; }
    int Ly () const { return _ysize; }
    int Lz () const { return _zsize; }
    int size () const { return siteN; }
    int bonds () const { return bondN; }
    int xi (int i) const { return site[i].xi(); }
    int yi (int i) const { return site[i].yi(); }
    int zi (int i) const { return site[i].zi(); }
    int nb  (int isite, int inb) const { return site[isite].nb(inb); } // neighbor site
    int nbs (int i) const { return site[i].nbs(); } // neighbor number for this site
    int opp_nbi (int isite, int inb) const { return site[isite].opp(inb); } // opposite neighbor index
    const SingleSite& operator() (int i) const { return site[i]; }
    // Get the index from the coordiate
    int index (int x, int y, int z) const { return get_index (x+_xcenter, y+_ycenter, z+_zcenter); }
    int xmin () const { return -_xcenter; }
    int xmax () const { return _xsize - _xcenter - 1; }
    int ymin () const { return -_ycenter; }
    int ymax () const { return _ysize - _ycenter - 1; }
    int zmin () const { return -_zcenter; }
    int zmax () const { return _zsize - _zcenter - 1; }
    // Return lines across the original point
    const std::vector<int>& xoline () const { return _xoline; }
    const std::vector<int>& yoline () const { return _yoline; }
    const std::vector<int>& zoline () const { return _zoline; }
    int xoline (int i) const { return _xoline[i]; }
    int yoline (int i) const { return _yoline[i]; }
    int zoline (int i) const { return _zoline[i]; }
    void set_center (int xcenter, int ycenter, int zcenter);
    void write (std::string file) const;
    void read (std::string file);
    void init (int xsize, int ysize, int zsize, int xc, int yc, int zc, bool bcx, bool bcy, bool bcz);
    void clear ();
    int index0 (int x, int y = 0, int z = 0) const;
    double dist (int i, int j) const { int dx = xi(i)-xi(j), dy = yi(i)-yi(j), dz = zi(i)-zi(j); return sqrt(dx*dx+dy*dy+dz*dz); }

  private:
    std::vector<SingleSite> site;
    int _xsize, _ysize, _zsize;
    int _xcenter, _ycenter, _zcenter;
    int wx, wy, wz; // weights of x, y, z from coordinate index to site index
    int siteN, bondN;
    bool bound_cond_x, bound_cond_y, bound_cond_z; // boundary condition: 0 for PERIODIC, 1 for OPEN
    std::vector<int> _xoline, _yoline, _zoline;

  protected:
    void find_neighbors ();
    void set_opp_neighbor_index ();
    void set_bond_indices ();
    int get_index (int xi, int yi, int zi) const;
    int true_coor (int i, int len, bool bc) const;
};

/*Lattice :: Lattice (int xsize, int ysize, int zsize, BoundCond bc, int xc, int yc, int zc)
{
  init (xsize, ysize, zsize, xc, yc, zc, bc, bc, bc);
}*/

Lattice :: Lattice (int xsize, int ysize, int zsize, bool bcx, bool bcy, bool bcz, int xc, int yc, int zc)
{
  init (xsize, ysize, zsize, xc, yc, zc, bcx, bcy, bcz);
}

void Lattice :: clear ()
{
  site.clear(); _xoline.clear(); _yoline.clear(); _zoline.clear();
}

void Lattice :: init (int xsize, int ysize, int zsize, int xc, int yc, int zc, bool bcx, bool bcy, bool bcz)
{
  clear();
  if (xc == -1) xc = xsize/2;
  if (yc == -1) yc = ysize/2;
  if (zc == -1) zc = zsize/2;
  bound_cond_x = bcx, bound_cond_y = bcy, bound_cond_z = bcz;
  _xsize = xsize;
  _ysize = ysize;
  _zsize = zsize;
  wx = 1;
  wy = xsize;
  wz = xsize * ysize;
  siteN = xsize * ysize * zsize;
  site.resize (siteN);
  find_neighbors();
  set_opp_neighbor_index();
  set_center (xc, yc, zc);
  set_bond_indices();
}

void Lattice :: set_bond_indices ()
{
  int bi = 0;
  for(int i = 0; i < siteN; ++i)
    for(int nbi = 0; nbi < site[i].nbs(); ++nbi) {
      int nbsite = site[i].nb(nbi);
      if (nbsite > i)
        site[i].bond[nbi] = bi++;
      else {
        int opp = site[i].opp(nbi);
        site[i].bond[nbi] = site[nbsite].bond[opp];
      }
    }
  bondN = bi;
}

void Lattice :: set_center (int xcenter, int ycenter, int zcenter)
{
  _xcenter = xcenter;
  _ycenter = ycenter;
  _zcenter = zcenter;
  for(int ix = 0; ix < _xsize; ix++)
  for(int iy = 0; iy < _ysize; iy++)
  for(int iz = 0; iz < _zsize; iz++) {
    int index = get_index (ix, iy, iz);
    site[index].xi (site[index].absxi() - xcenter);
    site[index].yi (site[index].absyi() - ycenter);
    site[index].zi (site[index].abszi() - zcenter);
  }
  // Store the sites in the line crossing to original point
  for(int i = 0; i < siteN; i++) {
    // For x line
    if (site[i].yi() == 0 && site[i].zi() == 0)
      _xoline.push_back (i);
    // For y line
    if (site[i].xi() == 0 && site[i].zi() == 0)
      _yoline.push_back (i);
    // For z line
    if (site[i].xi() == 0 && site[i].yi() == 0)
      _zoline.push_back (i);
  }
}

void Lattice :: find_neighbors ()
{
  for(int ix = 0; ix < _xsize; ix++)
  for(int iy = 0; iy < _ysize; iy++)
  for(int iz = 0; iz < _zsize; iz++) {
    int index = get_index (ix, iy, iz);
    site[index]._absxi = ix;
    site[index]._absyi = iy;
    site[index]._abszi = iz;

    site[index].add_neighbor (get_index (ix-1, iy, iz));
    site[index].add_neighbor (get_index (ix+1, iy, iz));
    site[index].add_neighbor (get_index (ix, iy-1, iz));
    site[index].add_neighbor (get_index (ix, iy+1, iz));
    site[index].add_neighbor (get_index (ix, iy, iz-1));
    site[index].add_neighbor (get_index (ix, iy, iz+1));
  }
}

void Lattice :: set_opp_neighbor_index ()
{
  for(int i = 0; i < siteN; i++)
    for(int inb = 0; inb < site[i].nbs(); inb++) {
      int nb_site = site[i].nb(inb);
      int opp_nb = 0;
      while (site[nb_site].nb(opp_nb) != i) {
#ifdef DEBUG_MODE
        if (opp_nb >= site[nb_site].nbs()) {
          std::cout << "*error: Lattice:: set_opp_neighbor_index: cannot find opposite neighbor index\n";
          std::cout << " site = " << i << "\n which neighbor: " << inb << "\n"; exit(100);
        }
#endif
          opp_nb++;
      }
      site[i].set_opp( inb, opp_nb);
    }
}

int Lattice :: true_coor (int i, int len, bool bc) const
{
  if (bc == 1) {
    if (i < 0 || i >= len) return -1;
  }
  else {
    // If the length is less than 3, it reduces to open boundary condition.
    if (len <= 2)
      if (i < 0 || i >= len) return -1;
    // If the coordinate out of lattice, shift back by periodic boundary condition.
    while (i < 0) i += len;
    while (i >= len) i -= len;
  }
  return i;
}

int Lattice :: get_index (int ix, int iy, int iz) const
{
  ix = true_coor (ix, _xsize, bound_cond_x);
  if (ix == -1) return -1;
  iy = true_coor (iy, _ysize, bound_cond_y);
  if (iy == -1) return -1;
  iz = true_coor (iz, _zsize, bound_cond_z);
  if (iz == -1) return -1;
  return ix + iy*wy + iz*wz;
}

int Lattice :: index0 (int x, int y, int z) const
{
  int ind = get_index (x, y, z);
  if (ind == -1) {
    std::cout << "Error: Lattice :: index0: no such coordinate\n";
    exit(100);
  }
  return ind;
}

void Lattice :: write (std::string file) const
{
  std::ofstream ofs ((file+".lattice").c_str());
  if (!ofs) {
    std::cout << "Error open file: Lattice: write: " << (file+".lattice") << "\n";
    exit(100);
  }
  ofs << _xsize << " " << _ysize << " " << _zsize << " " << _xcenter << " " << _ycenter << " " << _zcenter << " " << bound_cond_x << " " << bound_cond_y << " " << bound_cond_z;
  ofs.close();
}

void Lattice :: read (std::string file)
{
  std::ifstream ifs ((file+".lattice").c_str());
  if (!ifs) {
    std::cout << "Error open file: Lattice: read: " << (file+".lattice") << "\n";
    exit(100);
  }
  int xsize, ysize, zsize, xcenter, ycenter, zcenter;
  bool bcx, bcy, bcz;
  ifs >> xsize >> ysize >> zsize >> xcenter >> ycenter >> zcenter >> bcx >> bcy >> bcz;
  ifs.close();
  clear();
  init (xsize, ysize, zsize, xcenter, ycenter, zcenter, bcx, bcy, bcz);
}
#endif
