#ifndef OBSERVABLE2_HPP_CMC
#define OBSERVABLE2_HPP_CMC
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>

class Observable2
{
  public:
    typedef std::pair<int,int> Pair;
    Observable2 (const std::string& file, bool append=false);
    void measure (int i, int j, double a) { temp[Pair(i,j)] += a; }//++count; }
    void write (double count);
    void clear () { temp.clear(); }//count = 0.; }

  private:
    std::map<Pair, double> temp;
    std::string FILENAME;
    //double count;
};

Observable2 :: Observable2 (const std::string& file, bool append) : FILENAME (file+".obs3")//, count (0.)
{
  if (!append) remove (FILENAME.c_str());
}

void Observable2 :: write (double count)
{
  // Write the obserables' names
  std::ofstream ofs (FILENAME.c_str(), std::ios::app);
  if (!ofs) {
    std::cout << "Error: Observable2: write: cannot open file '" << FILENAME << "'\n";
    exit(100);
  }
  ofs << std::setprecision(18);
  for(std::map<Pair,double>::iterator it = temp.begin(); it != temp.end(); ++it)
    ofs << it->first.first << " " << it->first.second << " " << ((it->second)/count) << "\n";
  ofs.close();
  temp.clear();
}
#endif
