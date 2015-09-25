#ifndef SINGLESITE_HPP
#define SINGLESITE_HPP
class Lattice;

class SingleSite
{
  private:
    int  _absxi, _absyi, _abszi;
    int  _xi, _yi, _zi; // coordinate index
    bool _si;	 // for honeycomb lattice
    int  _neighbor[6];
    int  neighbor_num;
    int  opp_nbi[6]; // opposite neighbor index
    int  bond[6]; // bond index
  friend class Lattice;

  public:
         SingleSite () { neighbor_num = 0; };

    int  nb (int i) const { return _neighbor[i]; }
    int  nbs () const { return neighbor_num; }
    int  xi () const { return _xi; }
    int  yi () const { return _yi; }
    int  zi () const { return _zi; }
    bool si () const { return _si; }
    int  opp (int i) const { return opp_nbi[i]; }
    int  bi (int i) const { return bond[i]; }
    int  absxi () const { return _absxi; }
    int  absyi () const { return _absyi; }
    int  abszi () const { return _abszi; }

    void xi (int xi_) { _xi = xi_; }
    void yi (int yi_) { _yi = yi_; }
    void zi (int zi_) { _zi = zi_; }
    void si (int si_) { _si = si_; }
    void set_opp (int i, int opp_i) { opp_nbi[i] = opp_i; }
    void add_neighbor (int neighbor_);
};

void SingleSite :: add_neighbor (int neighbor_)
{
  if (neighbor_ != -1) {
    _neighbor [neighbor_num] = neighbor_;
    ++neighbor_num;
  }
}
#endif
