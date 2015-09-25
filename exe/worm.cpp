#include <iostream>
#include <fstream>
//using namespace std;
int dbstart = 35800;
int ii = 0;
//#define DEBUG_MODE
//#define DEBUG_TRACK_MODE
#define __DEBUG_TRACK dbi++; std::cout << dbi << "  " << __func__ << "\n";
#ifdef DEBUG_TRACK_MODE
bool dbbool = true;
std::ofstream dbo ("debug.info");
long dbi = 0;
#endif

#define LOCAL_OPTIMAL
#define NN_INTERACTION
//#define SPIN_HALF_XXZ_MODEL

#define TRACING_N_HOP
#define TRACING_DIAG_E
#define TRACING_N
#define TRACING_NA
#define TRACING_WINDING_NUMBER
//#define MEASURE_CORRELATION

#include "Simulation.hpp"
#include "ParaBox.hpp"

int main (int argc, char *argv[])
{
  std::string parafile = argv[1];
  ParaBox para (parafile);
  if (!para.has ("BETA")) {
    double T = para["TEMPERATURE"];
    para.insert ("BETA", 1./T);
  }
  Simulation qmc (para);
  if (para.has("FIND_N_STEPS") && int(para["FIND_N_STEPS"]) >= 0) {
    para.set ("BH_MU", qmc.get_mu());
    std::ofstream ofs (parafile.c_str());
    ofs << para;
    ofs.close();
  }
  qmc.thermalization ();
  qmc.run ();
  //qmc.get_mu ();
  //qmc.find_mu (para["CANONICAL_N"], para["EOFFSET"], para["BH_MU"], para["SUB_BIN_SIZE"], para["BUFFER_SIZE"]);
}
