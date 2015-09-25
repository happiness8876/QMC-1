#ifndef SIMULATION_HPP_CMC
#define SIMULATION_HPP_CMC
#include <vector>
#include <fstream>
#include <sstream>
#include "InfoWLdiagram.hpp"
#include "ParaBox.hpp"
#include "Histogram.hpp"
#include "Observable.hpp"
#include "Observable1.hpp"
#include "Observable2.hpp"

class Simulation
{
  public:
    Simulation (const ParaBox& para);
    void thermalization ();
    void start (int bins, int bufs, int runs);
    void test ();
    void run ();
    double find_mu (double mu1, double mu2);
    double get_mu ();
    const InfoWLdiagram& wld () const { return _wld; }
    double quick_estimate_mu (double N, double d);

  private:
    ParaBox _para;
    InfoWLdiagram _wld;
    std::string LOGFILE, OLDLOG;
    void update_histogram (Histogram& his_n_any, std::vector<Histogram>& hiss_n, int N);
    void write_log (int n) const;
    double estimate_N (double mu, int samples=10000);
};

Simulation :: Simulation (const ParaBox& para)
: _wld (para["LX"], para["LY"], para["LZ"],para["BETA"], 0.5*(double)para["BH_U"], para["BH_t"], para["BH_MU"], para["BH_Vnn"], para["VX_TRAP"], para["VY_TRAP"], para["VZ_TRAP"], para["EOFFSET"], para["BOUNDARY_CONDITION_X"], para["BOUNDARY_CONDITION_Y"], para["BOUNDARY_CONDITION_Z"], para["INIT_n"], para["nMAX"], para["RandSeed"])
, _para (para)
{
  LOGFILE = std::string(_para["RAWDATA_DIR"]) + "/" + std::string(_para["JOB_NAME"]) + ".log";
  if (_para["APPEND"]) {
    std::ifstream ifs (LOGFILE.c_str());
    /*if (!ifs) {
      std::cout << "Error: Simulation :: Simulation: cannot open file '" << LOGFILE << "'\n";
      exit(10);
    }*/
    if (ifs.good()) {
      std::stringstream sstr;
      sstr << ifs.rdbuf();
      OLDLOG = sstr.str();
    }
    ifs.close();
  }
  int Lx = _para["LX"], Ly = _para["LY"], Lz = _para["LZ"];
  _wld.set_subsystem_A (_wld.latt().sub(Lx/2,Ly,Lz));
}

void Simulation :: run ()
{
  std::string name = _para["JOB_NAME"], locat = _para["RAWDATA_DIR"], bkdir = _para["BACKUP_DIR"];
  double runs = _para["RUN_SIZE"];
  int subbinsize = _para["SUB_BIN_SIZE"], bufsize = _para["BUFFER_SIZE"], append = _para["APPEND"];
  int Lx = _para["LX"], Ly = _para["LY"], Lz = _para["LZ"];
  std::string bkname = bkdir+"/"+name;
  if (append == 1) _wld.read (bkname);
  else append = 0;
  int writesize = subbinsize * bufsize;
  int canonicalN = _para.has("CANONICAL_N") ? _para["CANONICAL_N"] : -1;
  int NAmean = canonicalN / 2;

  std::cout << "Begin sampling...\n";
  Observable obs (locat+"/"+name, "N N2 Nprob", subbinsize, bufsize, append);
  Observable obs2 (locat+"/"+name+"_2", "NA Ed Nhop W2", subbinsize, bufsize, append);
  Observable obs3 (locat+"/"+name+"_3", "dNA-5 dNA-4 dNA-3 dNA-2 dNA-1 dNA0 dNA1 dNA2 dNA3 dNA4 dNA5", subbinsize, bufsize, append);
#ifdef MEASURE_CORRELATION
  Observable1<double> obs_bb (locat+"/"+name, append);
#endif

  for(int i = 0; i < runs; ++i) {
    int count = 0;
    for(int k = 0; k < writesize; ++k) {
      _wld.updateZ();
      int N = _wld.N();
      obs >> N >> (N*N) >> (N == canonicalN) >> "n";
      if (canonicalN < 0 || N == canonicalN) {
        obs2 >> _wld.NA() >> _wld.Ed() >> _wld.Nhop() >> _wld.W2() >> "n";
        // NA sectors
        std::map<int, double> NA_weight = _wld.NA_weight();
        for(int dNA = -5; dNA <= 5; ++dNA) {
          int NA = NAmean + dNA;
          obs3 >> NA_weight[NA];
        }
        obs3 >> "n";
#ifdef MEASURE_CORRELATION
        // b_i^dagger b_j
        for(std::map<double,int>::const_iterator it = _wld.bb().begin(); it != _wld.bb().end(); ++it)
          obs_bb.measure (it->first, it->second);
#endif
        ++count;
      }
    }
#ifdef MEASURE_CORRELATION
    obs_bb.write(count);
#endif
    _wld.write (bkname);
    write_log (i);
  }
  std::cout << "...complete\n";
  _wld.info_plot ("test");
}

void Simulation :: thermalization ()
{
  std::string name = _para["JOB_NAME"], locat = _para["RAWDATA_DIR"], bkdir = _para["BACKUP_DIR"];
  int thermals = _para["THERMALIZE_STEPS"], subbinsize = _para["SUB_BIN_SIZE"], bufsize = _para["BUFFER_SIZE"], append = _para["APPEND"];
  int writesize = subbinsize * bufsize;
  std::string bkname = bkdir+"/"+name;
  if (append == 1) return;
  else if (append == 2) _wld.read (bkname);
  append = 0;
  double T = 1./(double)_para["BETA"];

  Observable obs (locat+"/"+name+"_n_equi", "N E", subbinsize, bufsize, append);
  std::cout << "Begin thermalization...\n";
  for(int i = 0; i < thermals; ++i) {
    for(int k = 0; k < writesize; ++k) {
      _wld.updateZ();
      obs >> _wld.N() >> (_wld.Ed() - _wld.Nhop()*T) >> "n";
    }
    _wld.write (bkname);
  }
  std::cout << "...complete\n";
}

void Simulation :: write_log (int out_i) const
{
  std::ofstream ofs (LOGFILE.c_str());
  if (!ofs) {
    std::cout << "Error: Simulation: write_log: cannot open file '" << LOGFILE << "'\n";
    exit(100);
  }
  ofs << OLDLOG;
  ofs << "Seed = " << _wld.seed() << "\n";
  ofs << "Times of output = " << out_i << "\n";
  ofs.close();
}

void Simulation :: test ()
{
  for(int i = 0; i < 1000; ++i)
    _wld.updateZ();
  //std::cout << _wld.NA() << "\n";
  _wld.info_plot ("test");
  _wld.NA_weight ();
/*  Time_slice timeslice (_wld);
  do {
    for(int i = 0; i < _wld.latt().size(); ++i)
      std::cout << timeslice.it(i)->z.nbef() << "  ";//":" << timeslice.it(i)->time() << "  ";
    std::cout << ": " << timeslice.dt() << "\n";
  }while (timeslice.next());*/
}

double Simulation :: get_mu ()
{
  //double mu = _para["BH_MU"];
  double mu = _para.has("BH_MU") ? double(_para["BH_MU"]) : quick_estimate_mu (_para["CANONICAL_N"], 0.1);
  std::string locat = _para["RAWDATA_DIR"];
  mu = find_mu (mu, mu+0.1);
  //std::ofstream ofs ((locat+"/U_"+std::string(_para["BH_U"])+"_mu.dat").c_str());
  //std::cout << double(_para["BH_U"]) << " " << mu << "\n";
  //ofs.close();
  _wld.set_mu(mu);
  return mu;
}

double Simulation :: quick_estimate_mu (double N, double d)
{
  double mu = 0., muest = 0.;
  bool signpre, sign;
  int oscillations = 0, oslimit = 100;
  while (oscillations < oslimit) {
    _wld.updateZ();
    if (_wld.N() < N) { mu += d; sign = false; }
    else { mu -= d; sign = true; }
    if (signpre != sign) { ++oscillations; signpre = sign; muest += mu; }
    _wld.set_mu (mu);
  }
  return muest/double(oscillations);
}

double Simulation :: find_mu (double mu1, double mu2)
{
  double N=_para["CANONICAL_N"], precision=0.01;//, mu_init=_para["BH_MU"];
  int samples=_para["FIND_N_STEPS"];
  //bool app=_para["APPEND"];
  std::string name = _para["JOB_NAME"], bkdir = _para["BACKUP_DIR"];
  std::string bkname = bkdir+"/"+name;
  //if (app) _wld.read (bkname);
  double SIZE = int(_para["LX"])*int(_para["LY"])*int(_para["LZ"]);

  for(int i = 0; i < samples; ++i) _wld.updateZ();
  //double mu1 = mu_init;
  double N1 = estimate_N (mu1, samples);
  std::cout << "mu = " << mu1 << ";  f = " << (N1-N) << "; N = " << N1 << "\n";
  if (fabs (N1 - N) < precision) return mu1;
  //double mu2 = mu_init+0.01;
  double N2 = estimate_N (mu2, samples);
  std::cout << "mu = " << mu2 << ";  f = " << (N2-N) << "; N = " << N2 << "\n";
  if (fabs (N2 - N) < precision) return mu2;
  double mu3 = mu2 - (N2-N)*(mu2-mu1)/(N2-N1);
  double N3 = estimate_N (mu3, samples);
  std::cout << "mu = " << mu3 << ";  f = " << (N3-N) << "; N = " << N3 << "\n";
  if (fabs (N3 - N) < precision) return mu3;
  while (fabs (N3 - N) > precision) {
    mu1 = mu2; mu2 = mu3; N1 = N2; N2 = N3;
    mu3 = mu2 - (N2-N)*(mu2-mu1)/(N2-N1);
    N3 = estimate_N (mu3, samples);
    std::cout << "mu = " << mu3 << ";  f = " << (N3-N) << "; N = " << N3 << "\n";
  }
  _wld.write (bkname);
  return mu3;
}

double Simulation :: estimate_N (double mu, int samples)
{
  _wld.set_mu(mu);
  //int skips = samples/2;
  //for(int i = 0; i < skips; ++i) _wld.updateZ();
  double ntemp = 0.;
  for(int i = 0; i < samples; ++i) {
    _wld.updateZ();
    ntemp += _wld.N();
  }
  ntemp /= (double)(samples);
  return ntemp;
}
#endif
