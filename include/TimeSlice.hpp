#ifndef TIME_SLICE_HPP_CMC
#define TIME_SLICE_HPP_CMC
#include <algorithm>
#include <vector>
#include "InfoWLdiagram.hpp"

class Time_slice
{
  public:
    Time_slice (const InfoWLdiagram& wld, const std::vector<int>& sites) { init(wld,sites); }
    Time_slice (const InfoWLdiagram& wld) { std::vector<int> s; for(int i = 0; i < wld.latt().size(); ++i) s.push_back(i); init(wld,s); }
    void init (const InfoWLdiagram& wld, const std::vector<int>& sites);
    const InfoWLdiagram::const_Iter& it (int i) const { return its[i]; }
    const InfoWLdiagram::const_Iter& it_first () const { return its_sort.front(); }
    double dt () const { return _dt; }
    bool next ();
  private:
    double END_TIME, _dt, tpre;
    std::vector<InfoWLdiagram::const_Iter> its; // Store the iterators by site
    std::vector<InfoWLdiagram::const_Iter> its_sort; // Store the iterators sorted by time
    static bool compare_time (InfoWLdiagram::const_Iter it1, InfoWLdiagram::const_Iter it2) { return (it1->time() < it2->time()); }
    void sort_first (std::vector<InfoWLdiagram::const_Iter>& a);
};

void Time_slice :: init (const InfoWLdiagram& wld, const std::vector<int>& sites)
{
  END_TIME = wld.end_time();
  int siteN = wld.latt().size();
  its.clear();
  its_sort.clear();
  for(int i = 0; i < sites.size(); ++i) {
    InfoWLdiagram::const_Iter it = ++wld.it_beg(sites[i]);
    its.push_back (it); // Iterators for a state of specific time
    its_sort.push_back (it); // Iterators sorted in time
  }
  std::sort (its_sort.begin(), its_sort.end(), compare_time);
  _dt = its_sort.front()->time();
  tpre = _dt;
}

void Time_slice :: sort_first (std::vector<InfoWLdiagram::const_Iter>& a)
{
  InfoWLdiagram::const_Iter first = a.front();
  double t = a.front()->time();
  std::vector<InfoWLdiagram::const_Iter>::iterator it = a.begin(), it_run = ++a.begin();
  while (it_run != a.end() && t > (*it_run)->time()) {
    *it = *it_run;
    ++it_run;
    ++it;
  }
  *it = first;
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
