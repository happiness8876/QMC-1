#ifndef RANDGEN_HPP
#define RANDGEN_HPP
#include <time.h>
#include <math.h>
//#include "mt19937.cpp"
//#include "randmt.c"
#include "MersenneTwister.h"

class RandGen
{
   private:
      double random;
      unsigned long _seed;
      MTRand gen;

   public:
      RandGen (unsigned long seed_ = time(NULL));

      double myrand () { return gen.rand(); }
      double uniform ();
      bool choice2 () { return (myrand() < 0.5); }
      int choice (int N);
      unsigned long seed () const { return _seed; }
      void save (const std::string& file);
      void load (const std::string& file);
};

void RandGen :: save (const std::string& file)
{
  MTRand::uint32 temp[MTRand::SAVE];
  gen.save (temp);
  std::ofstream ofs ((file+".rand").c_str(), std::ios::binary);
  if (!ofs) {
    std::cout << "Error: RandGen :: save: cannot open file '" << file << "\n";
    exit(10);
  }
  ofs.write ((char*)temp, MTRand::SAVE*sizeof(MTRand::uint32));
  ofs.close();
}

void RandGen :: load (const std::string& file)
{
  MTRand::uint32 temp[MTRand::SAVE];
  std::ifstream ifs ((file+".rand").c_str(), std::ios::binary);
  if (!ifs) {
    std::cout << "Error: RandGen :: save: cannot open file '" << file << "\n";
    exit(10);
  }
  ifs.read ((char*)temp, MTRand::SAVE*sizeof(MTRand::uint32));
  ifs.close();
  gen.load (temp);
}

inline double RandGen :: uniform ()
{
   random = myrand();
   // Exclude the case of random = 0 and 1
   while (random == 0.0 || random == 1.0)
      random = myrand();
   return random;
}

inline int RandGen :: choice (int N)
{
   random = myrand();
   // Exclude the case of random = 1
   while (random == 1.0)
      random = myrand();
   return int(random * N);
}

inline RandGen :: RandGen (unsigned long seed_)
{
   // Set the seed for mt19937.
   //setmt19937 (seed_);
   //init_randmt (seed_);
   if (seed_ == 0) gen.seed ();
   else gen.seed (seed_);
   // Store the seed.
   _seed = gen.get_seed();
}
#endif
