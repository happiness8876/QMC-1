#ifndef HISTOGRAM_HPP_CMC
#define HISTOGRAM_HPP_CMC
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iomanip>

class Histogram
{
  public:
    Histogram () : ci(0), count(0.) {}
    //Histogram (std::string name, int bufsize, bool app);
    void add (int n, double times = 1.);
    int size () const { return his.size(); }
    void write (const std::string& file_name);

  private:
    std::vector<double> his;
    double count;
    int ci, buffer_size;
    std::string FILENAME;
};

/*Histogram :: Histogram (std::string name, int bufsize, bool app)
 : buffer_size(bufsize), ci(0), count(0.), FILENAME(name+".his")
{
  his.resize (100, 0.);
  if (app) {
    std::ifstream ifs (FILENAME.c_str());
    if (ifs.good()) {
      int n, times;
      ifs >> n >> times;
      while (!ifs.eof()) {
        add (n, times);
        ifs >> n >> times;
      }
    }
    ifs.close();
  }
}*/

void Histogram :: add (int n, double times)
{
  if (n >= his.size()) his.resize (n+1, 0.);
  his[n] += times;
  ++count;
  /*if (++ci == buffer_size) {
    count += (double)ci;
    //write(FILENAME);
    ci = 0;
  }*/
}

void Histogram :: write (const std::string& file_name)
{
  std::ofstream ofs (file_name.c_str());
  if (!ofs) {
    std::cout << "Error: Histogram: write: cannot open file '" << file_name << "'\n";
    exit(100);
  }
  ofs << std::setprecision(18);
  for(int i = 0; i < his.size(); ++i) {
    if (his[i] != 0.)
      ofs << i << " " << his[i]/count << "\n";
  }
  ofs.close();
}
#endif
