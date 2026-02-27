#ifndef COMMON_HPP
#define COMMON_HPP

#include <cstdio>
#include <iostream> // std::cout, std::fixed
#include <iomanip>      // std::setprecision
#include <fstream>
#include <sstream>

#include <string>
#include <cstring>

#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cassert>

#include <algorithm>
#include <random>

#include <map>
#include <vector>
#include <set>
#include <stack>

#include <boost/format.hpp>   // used in newick format
#include <boost/algorithm/string.hpp>

#include "gzstream.h"


using namespace std;

// MK model was used for test, not suitable for CN evolution
// BOUNDT: bounded total copy number model, ignoring hidden haplotype-specific changes, not recommended for real data
// BOUNDA: bounded haplotype-specific copy number model, used in CNETML paper
// INFINITE was used for simulating Poisson process
enum MODEL {MK, BOUNDT, BOUNDA, DECOMP, INFINITE = 3};
enum CN_TYPE {ONLY_SEG, EXCLUDE_SEG, EXCLUDE_CHR, EXCLUDE_WGD, ALL}; // 0-ONLY_SEG: only use segments level CN changes (use only r1/r2); 1-EXCLUDE_SEG: exclude segment level CN changes (use only r3/r4/r5); 2-EXCLUDE_CHR: exclude chromosome level changes (use only r1/r2/r5) ; 3-EXCLUDE_WGD: exclude WGD level changes (use only r1/r2/r3/r4) ; 4-ALL: three types of mutations (use r1 to r5)

const int MAX_AGE = 100;

const int NUM_CHR = 22; // only consider autosome for now, TODD: add ChrX and chrY in the future ??
const int NORM_PLOIDY = 2;
const int NORM_ALLElE_STATE = 4;    // state 4 represents 1/1 in haplotype-specific copy number model BOUNDA
const int NO_CHANGE_WGD = 0; // no change at WGD level
const int NO_CHANGE_HAPLOTYPE = 5;   // no change at chr-level and site-level


const int PRINT_PRECISION = 10;

const double SMALL_VAL = 1.0e-10;   // used to compare floats

const float WGD_CUTOFF = 3.0;    // genome ploidy to determine WGD

// key: chr, seg, copy_number
typedef map<int, map<int, int>> copy_number;
// typedef vector<vector<double>> lnl_table;

const string VERSION = "1.0";


// read-only
struct INPUT_PROPERTY{
  int Ns;
  int cn_max;
  int model;

  int is_total;
  int is_rcn;
  int is_bin;
  int incl_all;
};


inline void print_input_property(const INPUT_PROPERTY& p) {
    std::cout
        << "INPUT_PROPERTY("
        << "Ns=" << p.Ns
        << ", cn_max=" << p.cn_max
        << ", model=" << p.model
        << ", is_total=" << p.is_total
        << ", is_rcn=" << p.is_rcn
        << ", is_bin=" << p.is_bin
        << ", incl_all=" << p.incl_all
        << ")\n";
}


// obtained from input file
struct INPUT_DATA{
  int num_invar_bins;
  int num_total_bins;
  int seg_size;  // Nchar

  vector<int> sample_num_wgd;
  vector<vector<int>> sample_change_chr;
  vector<vector<int>> sample_change_site;
  vector<int> chr_max_change;   
  vector<int> site_max_change;
  vector<int> sample_max_cn;
  vector<double> sample_avg_cn;
  vector<map<int, vector<int>>> sample_chr_cn;    // per sample, map<chr, vector<copy numbers on the chr>>
};


inline std::ostream& operator<<(std::ostream& os, const INPUT_DATA& d) {
    os << "INPUT_DATA {\n";
    os << "  num_invar_bins = " << d.num_invar_bins << "\n";
    os << "  num_total_bins = " << d.num_total_bins << "\n";
    os << "  seg_size (Nchar) = " << d.seg_size << "\n";

    os << "  sample_num_wgd:\n";
    for (size_t i = 0; i < d.sample_num_wgd.size(); ++i) {
        os << "    sample " << i << ": "
           << d.sample_num_wgd[i] << "\n";
    }

    os << "  sample_change_chr:\n";
    for (size_t i = 0; i < d.sample_change_chr.size(); ++i) {
        os << "    sample " << i << ": ";
        for (size_t j = 0; j < d.sample_change_chr[i].size(); ++j) {
            os << d.sample_change_chr[i][j] << " ";
        }
        os << "\n";
    }

    os << "  sample_change_site:\n";
    for (size_t i = 0; i < d.sample_change_site.size(); ++i) {
        os << "    sample " << i << ": ";
        for (size_t j = 0; j < d.sample_change_site[i].size(); ++j) {
            os << d.sample_change_site[i][j] << " ";
        }
        os << "\n";
    }

    os << "  chr_max_change: ";
    for (size_t i = 0; i < d.chr_max_change.size(); ++i) {
        os << d.chr_max_change[i] << " ";
    }
    os << "\n";

    os << "  site_max_change: ";
    for (size_t i = 0; i < d.site_max_change.size(); ++i) {
        os << d.site_max_change[i] << " ";
    }
    os << "\n";

    os << "  sample_max_cn: ";
    for (size_t i = 0; i < d.sample_max_cn.size(); ++i) {
        os << d.sample_max_cn[i] << " ";
    }
    os << "\n";

    os << "  sample_avg_cn: ";
    for (size_t i = 0; i < d.sample_avg_cn.size(); ++i) {
        os << d.sample_avg_cn[i] << " ";
    }
    os << "\n";

    os << "  sample_chr_cn:\n";
    for (size_t i = 0; i < d.sample_chr_cn.size(); ++i) {
        os << "    sample " << i << ":\n";

        std::map<int, std::vector<int> >::const_iterator it;
        for (it = d.sample_chr_cn[i].begin();
             it != d.sample_chr_cn[i].end(); ++it) {

            os << "      chr " << it->first << ": ";
            for (size_t k = 0; k < it->second.size(); ++k) {
                os << it->second[k] << " ";
            }
            os << "\n";
        }
    }

    os << "}";
    return os;
}



// to store observed decomposition information for an observed copy number
// store all values for convenience of access and adaption
struct CN_CHANGE{
    // copy number values
    int cn_state;
    // int cnA;
    // int cnB;

    // decomposed changes
    // TODO: may extend to store haplotype-specific information in the future
    int num_wgd;
    int cn_change_chr;
    int cn_change_site; 


    vector<int> to_vector() const {
        return {
            cn_state,
            num_wgd,
            cn_change_chr,
            cn_change_site
        };
    }    

    bool operator<(const CN_CHANGE& other) const {
        return to_vector() < other.to_vector();
    }
};


inline std::ostream& operator<<(std::ostream& os, const CN_CHANGE& cc) {
    os << "CN_CHANGE("
       << "cn_state=" << cc.cn_state
       << ", num_wgd=" << cc.num_wgd
       << ", cn_change_chr=" << cc.cn_change_chr
       << ", cn_change_site=" << cc.cn_change_site
       << ")";
    return os;
}


// Read strings separated by space into a vector.
// num: the number of element in the string. If not given, assume it is the first number in the string
template <typename T>
void get_vals_from_str(vector<T>& vals, string str_vals, int num = 0){
    assert(str_vals != "");
    stringstream ss(str_vals);
    if(num==0)   ss >> num;
    for(int i = 0; i < num; i++){
        T d1;
        ss >> d1;
        vals.push_back(d1);
    }
    // cout << "vector from " << str_vals << ": ";
    // for(auto v : vals){
    //     cout << "\t" << v;
    // }
    // cout << endl;
}


// used for haplotype-specific copy number model with multiple levels
// when maximum haplotype change is 1, the number of states with each total copy number values is: 1 2 3 2 1, with maximum number is 4
// when maximum haplotype change is 2, the number of states with each total copy number values is: 1 2 3 4 3 2 1, with maximum number is 6
inline vector<int> make_peak_vector(int n) {
    vector<int> v;
    v.reserve(2 * n - 1);  // total length: 1..n..1

    // Increasing part: 1, 2,..., n
    for (int i = 1; i <= n; ++i) {
        v.push_back(i);
    }

    // Decreasing part: n-1, n-2,..., 1
    for (int i = n - 1; i >= 1; --i) {
        v.push_back(i);
    }

    return v;
}

#endif
