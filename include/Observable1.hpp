#ifndef OBSERVABLE1_HPP_CMC
#define OBSERVABLE1_HPP_CMC
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <map>

template <class Temp>
class Observable1
{
  public:
    Observable1 (const std::string& file, bool append=false);
    void measure (const Temp& i, double a) { temp[i] += a; }//++count; }
    void write (double count);
    void clear () { temp.clear(); }//count = 0.; }

  private:
    std::map<Temp, double> temp;
    std::string FILENAME;
    //double count;
};

template <class Temp>
Observable1<Temp> :: Observable1 (const std::string& file, bool append) : FILENAME (file+".obs1")//, count (0.)
{
  if (!append) remove (FILENAME.c_str());
}

template <class Temp>
void Observable1<Temp> :: write (double count)
{
  // Write the obserables' names
  std::ofstream ofs (FILENAME.c_str(), std::ios::app);
  if (!ofs) {
    std::cout << "Error: Observable1: write: cannot open file '" << FILENAME << "'\n";
    exit(100);
  }
  ofs << std::setprecision(18);
  for(typename std::map<Temp,double>::iterator it = temp.begin(); it != temp.end(); ++it)
    ofs << it->first << " " << ((it->second)/count) << "\n";
  ofs.close();
  temp.clear();
}
#endif
