#ifndef CHECKLIST_HPP_CMC
#define CHECKLIST_HPP_CMC
#include <vector>
#include <fstream>
#include <iostream>
#include <iterator>
#include <cstdlib>

class CheckList
{
  private:
    std::vector<int> true_list, _status;
    int SIZE;

  public:
    CheckList (int size = 0, int cap_ = 100);
    CheckList (std::string file) { read (file); }
    void add (int i);
    void clear ();
    int size () const { return SIZE; }
    int subsize () const { return true_list.size(); }
    int ele (int i) { return true_list[i]; }
    int ele (int i) const  { return true_list[i]; }
    bool has (int i) const { return _status[i]; }
    const std::vector<int>& contains () const { return true_list; }
    void write (std::string file) const;
    void read (std::string file);
};

CheckList operator- (const CheckList& biglatt, const CheckList& smalllatt)
{
  if (biglatt.size() != smalllatt.size()) {
    std::cout << "Error: CheckList: diff: The base-sizes do not match\n";
    exit(100);
  }
  CheckList difflatt (biglatt.size());
  for(int i = 0; i < biglatt.subsize(); ++i) {
    int site = biglatt.ele(i);
    if (!smalllatt.has (site))
      difflatt.add (site);
  }
  return difflatt;
}

inline CheckList :: CheckList (int size, int cap_)
{
   SIZE = size;
   _status.resize (SIZE);
   for(int i = 0; i < SIZE; i++)
     _status[i] = false;
   true_list.reserve (cap_);
}

inline void CheckList :: clear()
{
   for(int i = 0; i < true_list.size(); i++)
     _status[true_list[i]] = false;
   true_list.clear();
}

inline void CheckList :: add (int i)
{
  if (_status[i] == false) {
    _status[i] = true;
    true_list.push_back(i);
  }
}

inline void CheckList :: write (std::string file) const
{
  std::ofstream ofs ((file+".checklist").c_str());
  if (!ofs) {
    std::cout << "Error open file: CheckList: write: " << file << "\n";
    exit(100);
  }
  // Write the sizes of vectors
  ofs << SIZE << " " << true_list.size() << "\n";
  // Write the elements of vectors
  for(std::vector<int>::const_iterator it = true_list.begin(); it != true_list.end(); ++it)
    ofs << *it << " ";
  ofs.close();
}

inline void CheckList :: read (std::string file)
{
  std::ifstream ifs ((file+".checklist").c_str());
  if (!ifs) {
    std::cout << "Error open file: CheckList: read: " << file << "\n";
    exit(100);
  }
  // Read the sizes of vectors
  int tsize, ssize;
  ifs >> SIZE >> tsize;
  true_list.resize (tsize);
  _status.clear();
  _status.resize (SIZE, false);
  // Read the elements of vectors
  for(std::vector<int>::iterator it = true_list.begin(); it != true_list.end(); ++it) {
    ifs >> *it;
    _status[*it] = true;
  }
  ifs.close();
}
#endif
