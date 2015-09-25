#ifndef PROBTABLE_HPP
#define PROBTABLE_HPP
#include <sstream>
#include <algorithm>

class ProbTable
{
  private:
    int last;
    double weight [7];
    int _nbi [7];
    //std::list<Node>::iterator _it[neighborN+1];
    std::vector<double*> sort_pt;

  public:
    ProbTable () { last = -1; }
    void clear () { last = -1; sort_pt.clear(); }

    //const std::list<Node>::iterator& it (int i) const { return _it[i]; }
    int nbi (int i) const { return _nbi[i]; }

    void add (double w);
    void add (double w, int nbi) { add(w); _nbi[last] = nbi; }
    //void add (double w, const std::list<Node>::iterator& it_, int nbi_) { add(w); _it[last] = it_; _nbi[last] = nbi_; }
    int get_choice (double rand);

    void locally_optimal ();

    const std::string info () const;
};

inline bool smaller (double *pa, double *pb) { return (*pa < *pb); }

void ProbTable :: locally_optimal ()
{
  if (last == 0) { weight[0] = 1.; return; }

  //sort_pt.sort (smaller);
  std::sort(sort_pt.begin(), sort_pt.end(), smaller);

  double w_bounce = weight[0];
  double w_sum = 0.;
  for(int i = 0; i <= last; ++i) w_sum += weight[i];

  double y = 0., y_sum = 0.;
  w_sum -= *(sort_pt.front());
  // From smallest to bounce weight:
  std::vector<double*>::iterator it = sort_pt.begin();
  while (*it != &weight[0]) {
    y = (1. - y_sum) * (**it / w_sum);
    y_sum += y;
    **it = y;
    ++it;
    w_sum -= **it;
    }
  y = (1. - y_sum) * (**it / w_sum);

  // From bounce to largest weight:
  for(++it; it != sort_pt.end(); ++it)
    **it *= y / w_bounce;
  // Bounce weight:
  //  Bounce weight is the largest one.
  if (sort_pt.back() == &weight[0])
    weight[0] = 1. - y_sum;
  else
    weight[0] = 0.;
}

inline void ProbTable :: add (double w)
{
#ifdef DEBUG_MODE
if (w < 0) { std::cout << "ProbTable: weight cannot be negtive.\n"; exit(100); }
#endif
  ++last;
  weight[last] = w;
  sort_pt.push_back (&weight[last]);
}

int ProbTable :: get_choice (double rand)
{
  for(int i = 1; i <= last; ++i)
    weight[i] += weight[i-1];
  rand *= weight[last];
  for(int i = 0; i <= last; ++i)
    if (weight[i] >= rand) return i;

  // Search nothing
  for(int i = 0; i <= last; ++i)
    std::cout << rand << ": " << weight[i] << "    ";
  std::cout << "\n";
  std::cout << "*error: ProbTalbe:: get_choice: get nothing from table\n"; exit(100);
  return -1;
}

const std::string ProbTable :: info () const 
{
  std::stringstream sstr;
  sstr << std::endl;
  for(int i = 0; i <= last; ++i)
    sstr << "P" << i << " = " << weight[i] << std::endl;
  sstr << std::endl;
  return sstr.str();
}
#endif
