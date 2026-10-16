#include "stats.hpp"


// Seed the generator. A seed of 0 (the programs' default for --seed) means
// "no seed given": one is derived from the clock and the process id, so two
// runs started in the same second still differ. Whichever seed ends up in use
// is printed and returned, so any log is enough to repeat the run with
// --seed <number> and the programs can write it into their outputs.
unsigned setup_rng(gsl_rng* r, unsigned set_seed){
  unsigned seed = set_seed;
  if(seed == 0){
    unsigned long t = (unsigned long)time(NULL);
    unsigned long pid = (unsigned long)getpid();
    // Keep the low 32 bits: --seed is an unsigned, so a wider value could not
    // be passed back, and GSL's mt19937 uses only 32 bits of the seed anyway.
    seed = (unsigned)((t * pid) & 0xffffffffUL);
    if(seed == 0) seed = 1;
  }
  cout << "seed:" << "\t" << seed << endl;
  gsl_rng_set(r, seed);
  return seed;
}


// wrapper for uniform
double runiform(gsl_rng* r, double a, double b){
  double myrandom = a + (b - a) * gsl_rng_uniform(r);

  while(myrandom == 0){
    myrandom = a + (b - a) * gsl_rng_uniform(r);
  }

  return myrandom;
}



// sample an element according to a probability vector
// gsl_ran_multinomial?
int rchoose(gsl_rng* r, const vector<double>& rates){
  //cout << "rchoose, rates:";
  //for(int i = 0; i < rates.size(); ++i) cout << "\t" << rates[i];
  //cout << endl;

  vector<double> p;
  double s = accumulate(rates.begin(), rates.end(), 0.0);

  for(int i = 0; i < rates.size(); ++i){
    p.push_back(rates[i]/s);
  }

  vector<double> psum(p.size(),0.0);
  partial_sum(p.begin(),p.end(),psum.begin());

  double u = gsl_rng_uniform(r);
  int ret = -1;
  for(int i = 0; i < rates.size(); ++i){
    if(u < psum[i]){
      ret = i;
      break;
    }
  }

  //cout << "u=\t" << u << "\t" << ret << endl;
  return ret;
}
