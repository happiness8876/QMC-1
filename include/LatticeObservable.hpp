#ifndef LATTICEOBSERVABLE_HPP_CMC
#define LATTICEOBSERVABLE_HPP_CMC
#include <string>
#include <fstream>
#include <cstdlib>

class Storage
{
  private:
    double _sum, lastVal;
    int lastCount;

  public:
    Storage () { _sum = 0.; lastVal = 0., lastCount = -1; }
    void clear () { _sum = 0.; }
    double sum () const { return _sum; }
    void save (int count, double val);
    void collect (int count);
    void reset (double a) { _sum = 0.; lastVal = a; lastCount = -1; }
};

inline void Storage :: save (int count, double val)
{
  int diff = count - lastCount;
  _sum += (double)diff * lastVal;
  lastCount = count;
  lastVal = val;
}

inline void Storage :: collect (int count)
{
  int diff = count - lastCount;
  _sum += (double)diff * lastVal;
  lastCount = 0;
}

class LatticeObservable
{
  public:
    LatticeObservable (std::string jobname, std::string obsname, int xsize, int ysize, int zsize, int merge_size=1, int buffer_size=1, bool append=false);
    template <class T> void measure (int site, T a) { obs[site].save(count, (double)a); }
    void end ();
    void to_buffer ();
    void write (const std::string& file_name);
    void to_file ();
    void init (int i, double a) { obs[i].reset(a); }

  private:
    std::string FILENAME;
    int SIZE, MERGE_SIZE, BUFFER_SIZE, count, bufcount;
    std::vector<Storage> obs;
    std::vector< std::vector<double> > buffer;
};

LatticeObservable :: LatticeObservable (std::string jobname, std::string obsname, int xsize, int ysize, int zsize, int merge_size, int buffer_size, bool append)
: SIZE(xsize*ysize*zsize), MERGE_SIZE(merge_size), BUFFER_SIZE(buffer_size), count(0), bufcount(0), obs(SIZE), FILENAME(jobname+".lattobs")
{
  if (!append) {
    // Backup the exist file; remove the old backup-file if existed
    std::string backup = FILENAME + ".backup";
    remove (backup.c_str());
    rename (FILENAME.c_str(), backup.c_str());
    // Write observables' names and lattice lengths to the data file
    std::ofstream ofs (FILENAME.c_str());
    if (!ofs) {
      std::cout << "Error: LatticeObservable: constructor: cannot open file '" << FILENAME << "'\n";
      exit(100);
    }
    ofs << obsname << " " << xsize << " " << ysize << " " << zsize << "\n";
    ofs.close();
  }
  // Initialize buffer
  buffer.resize (buffer_size);
  for(int i = 0; i < buffer_size; ++i)
    buffer[i].resize (SIZE);
}

inline void LatticeObservable :: end ()
{
  if (++count >= MERGE_SIZE) {
    count = 0;
    to_buffer();
  }
}

void LatticeObservable :: to_buffer ()
{
  for(int i = 0; i < SIZE; ++i) {
    obs[i].collect (MERGE_SIZE);
    buffer[bufcount][i] = obs[i].sum() / (double)(MERGE_SIZE);
    obs[i].clear();
  }
  if (++bufcount == BUFFER_SIZE) {
    write (FILENAME);
    bufcount = 0;
  }
}

void LatticeObservable :: write (const std::string& file_name)
{
  std::ofstream ofs (file_name.c_str(), std::ios::app);
  ofs << std::setprecision(18);
  for(std::vector<std::vector<double> >::const_iterator it = buffer.begin(); it != buffer.end(); ++it) {
    for(std::vector<double>::const_iterator itt = it->begin(); itt != it->end(); ++itt)
      ofs << *itt << " ";
    ofs << "\n";
  }
  ofs.close();
}
#endif
