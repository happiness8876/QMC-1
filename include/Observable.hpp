#ifndef OBSERVABLE_HPP_CMC
#define OBSERVABLE_HPP_CMC
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>

class Observable
{
  public:
    Observable (std::string jobname, std::string obsnames, int merge_size=1, int buffer_size=1, bool append=false, int data_size=-1);
    Observable& operator>> (double a);
    void operator>> (const std::string& any);
    void clear ();

  private:
    int MERGE_SIZE, BUFFER_SIZE, count, bufcount, obsi;
    std::string FILENAME;
    std::vector<double> temp;
    std::vector< std::vector<double> > buffer;
    std::vector<std::string> names;

    void to_buffer ();
    void write (const std::string& file_name, const std::vector< std::vector<double> >& buf);
};

Observable :: Observable (std::string jobname, std::string obsnames, int merge_size, int buffer_size, bool append, int data_size)
 : MERGE_SIZE(merge_size), BUFFER_SIZE(buffer_size), count(0), bufcount(0), FILENAME(jobname+".obs"), obsi(0)
{
  // Read the observables' names
  std::stringstream sstr (obsnames);
  std::string tempname;
  names.clear();
  while (!sstr.eof()) {
    sstr >> tempname;
    names.push_back (tempname);
  }
  if (!append) {
  // Backup the exist file; remove the backup-file if existed
    std::string backup = FILENAME + ".backup";
    remove (backup.c_str());
    rename (FILENAME.c_str(), backup.c_str());
    // Write observables' names to the data file
    std::ofstream ofs (FILENAME.c_str());
    if (!ofs) {
      std::cout << "Error: Observable: constructor: cannot open file '" << FILENAME << "'\n";
      exit(100);
    }
    for(int i = 0; i < names.size(); ++i)
      ofs << names[i] << " ";
    ofs << "\n";
    ofs.close();
  }
  // Set the number of observables. data_size != -1 means a vector-observable.
  if (data_size == -1) {
    buffer.resize (names.size());
    temp.resize (names.size(), 0.);
  }
  else {
    buffer.resize (data_size);
    temp.resize (data_size, 0.);
  }
}

inline Observable& Observable :: operator>> (double a)
// Measure the obsi-th observable value
{
  temp[obsi++] += a;
  return *this;
}

inline void Observable :: operator>> (const std::string& any)
{
// Finish a round of measurements for all observables
  obsi = 0;
  if (++count >= MERGE_SIZE) {
  // Save to buffer
    count = 0;
    to_buffer();
  }
}

void Observable :: to_buffer ()
{
  for(int i = 0; i < temp.size(); ++i) {
    buffer[i].push_back (temp[i]/(double)MERGE_SIZE);
    temp[i] = 0.;
  }
  if (++bufcount >= BUFFER_SIZE) {
    write (FILENAME, buffer);
    bufcount = 0;
    for(int i = 0; i < buffer.size(); ++i)
      buffer[i].clear();
  }
}

void Observable :: write (const std::string& file_name, const std::vector< std::vector<double> >& buf)
{
  std::ofstream ofs (file_name.c_str(), std::ios::app);
  if (!ofs) {
    std::cout << "Error: Observable: write: cannot open file '" << FILENAME << "'\n";
    exit(100);
  }
  ofs << std::setprecision(18);
  for(int i = 0; i < buf[0].size(); ++i) {
    for(int obsi = 0; obsi < buf.size(); ++obsi)
      ofs << buf[obsi][i] << " ";
    ofs << "\n";
  }
  ofs.close();
}

void Observable :: clear ()
{
  count = 0, bufcount = 0;
  for(int i = 0; i < temp.size(); ++i) {
    temp[i] = 0.;
    buffer[i].clear();
  }
}
#endif
