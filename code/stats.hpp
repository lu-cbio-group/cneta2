#ifndef STATS_HPP
#define STATS_HPP


#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>

#include <cstdio>
#include <iostream>

#include <numeric>  // for accumulate, partial_sum
#include <unistd.h>  // for getpid

#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <cmath>

#include <vector>


using namespace std;


// Pairing function: see https://en.wikipedia.org/wiki/Pairing_function
inline int pairInteger(int int1, int int2){
    if(int1 <= int2){
        return ((int1 + int2) * (int1 + int2 + 1)/2 + int2);
    }else{
        return ((int1 + int2) * (int1 + int2 + 1)/2 + int1);
    }
}

// Seed r with set_seed, or with a clock-derived seed when set_seed is 0.
// Returns the seed actually used (printed as "seed:<tab><n>" as well), so the
// caller can record it in its outputs.
unsigned setup_rng(gsl_rng* r, unsigned set_seed);


// A callable over a GSL generator: rng(n) returns a draw in [0, n). This is
// the form portable_shuffle and the tree-generation code take, so any
// generator can drive them, including one reserved for trees (cnets
// --tree_seed). Copying an Rng copies only the pointer; the generator's
// state lives in the gsl_rng it points at.
struct Rng {
    gsl_rng* r;
    explicit Rng(gsl_rng* r_) : r(r_) {}
    unsigned long operator()(unsigned long n) const { return gsl_rng_uniform_int(r, n); }
};


// factorial for choose(n,k)
inline int fact(int n){
  return (n == 1 || n == 0) ? 1 : fact(n - 1) * n;
}


// wrapper for uniform
double runiform(gsl_rng* r, double a, double b);



// sample an element according to a probability vector
// gsl_ran_multinomial?
int rchoose(gsl_rng* r, const vector<double>& rates);


#endif
