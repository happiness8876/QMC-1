#ifndef PARABOX_HPP_CMC
#define PARABOX_HPP_CMC
#include <fstream>
#include <map>
#include <cstdlib>
#include <iomanip>
#include <sstream>

class Cast
{
  // Casting a string to other types
  public:
    Cast (std::string t) { temp = t; }
    operator int() const { return atoi (temp.c_str()); }
    operator double() const { return atof (temp.c_str()); }
    operator bool() const { return bool (atoi (temp.c_str())); }
    operator std::string() const { return temp; }
    operator unsigned long() const { return strtoul (temp.c_str(), NULL, 0); }

  private:
    std::string temp;
};

class ParaBox
{
  public:
    ParaBox () {}
    ParaBox (std::string filename);
    Cast operator[] (const std::string& name) const; // Casting function
    bool has (std::string name) const { return _para.find(name) != _para.end(); }
    template <class T> void insert (const std::string& name, const T& val);
    template <class T> void set (const std::string& name, const T& val);

  friend void operator<< (std::ostream& os, ParaBox& p) {
    for(int i = 0; i < p._order.size(); ++i)
      os << std::setw(25) << std::left << p._order[i] << p._para[p._order[i]] << std::endl;
  }

  private:
    std::string _filename;
    std::map<std::string, std::string> _para;
    std::vector<std::string> _order;
};

ParaBox :: ParaBox (std::string filename)
{
  _filename = filename;
  std::ifstream ifs (filename.c_str());
  if (!ifs) {
    std::cout << "*error: Extractor :: constructor: failed to open '" << filename << "\n";
    exit(100);
  }
  std::string name, value;
  ifs >> name >> value;
  while (!ifs.eof()) {
    _order.push_back (name);
    _para[name] = value;
    ifs >> name >> value;
  }
  ifs.close();
}

Cast ParaBox :: operator[] (const std::string& name) const
{
  const std::map<std::string, std::string>::const_iterator it = _para.find(name);
  if (it != _para.end())
    return Cast (it->second);
  else {
    std::cout << "ParaBox : " << _filename << " : cannot find '" << name << "'\n";
    exit(100);
  }
}

template <class T> void ParaBox :: set (const std::string& name, const T& val)
{
  std::stringstream ss;
  ss << val;
  _para[name] = ss.str();
}

template <class T> void ParaBox :: insert (const std::string& name, const T& val)
{
  std::stringstream ss;
  ss << val;
  _para.insert ( std::pair<std::string, std::string>(name, ss.str()) );
}
#endif
