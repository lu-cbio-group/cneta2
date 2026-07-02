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
#include <limits>
#include <random>

#include <map>
#include <vector>
#include <set>
#include <stack>
#include <unordered_map>

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
const int NO_CHANGE_WGD = 0; // index of rate matrix for no change at WGD level
const int NO_CHANGE_HAPLOTYPE = 4;   // index of rate matrix for no change at chr-level and site-level
const int CHANGE_CHR = 1000;   // used to represent chromosome gain/loss in the change state, so that site-level changes can be adjusted accordingly to avoid double counting
const int NONE_STATE = -100;   // used to check if state table is initialized for a state

// limits on rate matrix dimension
const int MIN_CHANGE = -2;   // minimum copy number change at a level
const int MIN_CHANGE_HAPLOTYPE = -1;   // minimum haplotype-specific copy number change at a level

// const int MAX_SITE_CHANGE = 4;    // maximum copy number change at site level
// const int MAX_CHR_CHANGE = 2;    // maximum copy number change at chromosome level

const int PRINT_PRECISION = 10;

const double SMALL_VAL = 1.0e-10;   // used to compare floats

const float WGD_CUTOFF = 3.0;    // genome ploidy to determine WGD
const double WGD_WEIGHT = 0.001;   // used to represent very large negative log likelihood

// key: chr, seg, copy_number, shortname for convenience of access and printing
typedef map<int, map<int, int>> copy_number;
// typedef vector<vector<double>> lnl_table;

const string VERSION = "2.0";


// read-only, important input properties that will be used in multiple functions, wrapped in a struct for easier access and passing to functions
struct INPUT_PROPERTY{
    int Ns;
    int cn_max;

    int max_wgd;  // maximum number of WGD events allowed, determining the dimension of rate matrix for WGD
    int max_chr_change_haplotype;
    int max_site_change_haplotype;

    int model;

    int is_total;
    int is_rcn;
    int is_bin;
    int incl_all;
};

// to store decomposition information for an observed copy number, store all values for convenience of access and adaption
struct Transition {
    int from;
    int to;
    double rate;
    string label;
};

inline void print_input_property(const INPUT_PROPERTY& p) {
    std::cout
        << "INPUT_PROPERTY("
        << "Ns=" << p.Ns
        << ", cn_max=" << p.cn_max
        << ", max_chr_change_haplotype=" << p.max_chr_change_haplotype
        << ", max_site_change_haplotype=" << p.max_site_change_haplotype
        << ", model=" << p.model
        << ", is_total=" << p.is_total
        << ", is_rcn=" << p.is_rcn
        << ", is_bin=" << p.is_bin
        << ", incl_all=" << p.incl_all
        << ")\n";
}


// obtained from input copy number file
// wrapped in a struct for easier access and passing to functions during data parsing
struct INPUT_DATA{
    int num_invar_bins;
    int num_total_bins;
    int seg_size;  // Nchar

    // estimations based on averages or read from meta file, used to get maximum changes to switch on/off different rate matrices
    vector<int> sample_num_wgd;
    vector<int> chr_max_change;   
    vector<int> site_max_change;
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

    os << "}";
    return os;
}


// to store decomposition information for an observed copy number
// store all values for convenience of access and adaption
struct CN_CHANGE{
    int cn_state;  // encoded as either total CN (multiple indices in the rate matrix) or index in the haplotype-specific rate matrix

    // decomposed changes
    int num_wgd;
    int cn_change_chr;      // encoded as either total change (multiple indices in the rate matrix) or index in the haplotype-specific rate matrix
    int cn_change_site; 

    vector<int> to_vector() const {
        return {
            cn_state,
            num_wgd,
            cn_change_chr,
            cn_change_site
        };
    }    

    string to_string() const {
        return "CN_CHANGE("
            "cn_state="       + std::to_string(cn_state)       + ", "
            "num_wgd="        + std::to_string(num_wgd)        + ", "
            "cn_change_chr="  + std::to_string(cn_change_chr)  + ", "
            "cn_change_site=" + std::to_string(cn_change_site) + ")";
    }

    bool operator<(const CN_CHANGE& other) const {
        return to_vector() < other.to_vector();
    }
};


inline std::ostream& operator<<(std::ostream& os, const CN_CHANGE& cc) {
    os << cc.to_string();   
    return os;
}

// key: chr, seg, copy_number_change, shortname for convenience of access and printing
typedef map<int, map<int, CN_CHANGE>> copy_number_change;


// Read strings separated by space into a vector.
// num: the number of element in the string. If not given, assume it is the first number in the string
template <typename T>
inline void get_vals_from_str(vector<T>& vals, string str_vals, int num = 0){
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


inline int get_matrix_dim(int x) {
    return (x + 2) * (x + 2);
}


inline bool is_equal_vector(const vector<int>& bin1, const vector<int>& bin2){
    assert(bin1.size() == bin2.size());

    for(int i = 0; i < bin1.size(); i++){
        if(bin1[i] != bin2[i]){
            return false;
        }
    }
    return true;
}


inline string trim(const string& s) {
    size_t start = 0;
    while (start < s.size() && std::isspace(static_cast<unsigned char>(s[start]))) {
        ++start;
    }
    size_t end = s.size();
    while (end > start && std::isspace(static_cast<unsigned char>(s[end - 1]))) {
        --end;
    }
    return s.substr(start, end - start);
}


inline vector<string> split_string(const string& s, char delim) {
    vector<string> parts;
    string item;
    stringstream ss(s);
    while (std::getline(ss, item, delim)) {
        parts.push_back(trim(item));
    }
    return parts;
}


inline bool is_missing_field(const string& s) {
    const string t = trim(s);
    return t.empty() || t == "." || t == "NA" || t == "na" || t == "null";
}

inline int parse_int_strict(const string& s, const string& context) {
    const string t = trim(s);
    if (t.empty()) {
        throw runtime_error("Empty integer field in " + context);
    }

    size_t pos = 0;
    int value = 0;
    try {
        value = std::stoi(t, &pos);
    } catch (...) {
        throw runtime_error("Invalid integer '" + t + "' in " + context);
    }

    if (pos != t.size()) {
        throw runtime_error("Invalid integer '" + t + "' in " + context);
    }
    return value;
}

inline vector<int> parse_int_list(const string& s, char delim, const string& context) {
    vector<int> values;
    if (is_missing_field(s)) {
        return values;
    }

    vector<string> tokens = split_string(s, delim);
    values.reserve(tokens.size());

    for (const string& tok : tokens) {
        if (tok.empty()) {
            throw runtime_error("Empty entry in " + context);
        }
        values.push_back(parse_int_strict(tok, context));
    }
    return values;
}


template <typename T>
inline void print_vector(const vector<T>& data) {
    cout << "vector:" << endl;

    for (size_t j = 0; j < data.size(); ++j) {
        cout << data[j] << " ";
    }
    cout << endl;
    
}


template <typename T>
inline void print_nested_vector(const vector<vector<T>>& data) {
    cout << "nested vector:" << endl;

    for (size_t i = 0; i < data.size(); ++i) {
        cout << "\t[" << i << "]: ";
        print_vector(data[i]);
        cout << endl;
    }
}


template <typename T>
inline void print_nested_vector2(const vector<vector<vector<T>>>& data) {
    cout << "nested vector:" << endl;

    for (size_t i = 0; i < data.size(); ++i) {
        cout << "\t[" << i << "]: ";
        print_nested_vector<T>(data[i]);
        cout << endl;
    }
}


template <typename T>
inline void print_data_map(const map<int, vector<vector<T>>>& data) {
    cout << "Observed copy number data:" << endl;

    for (auto it = data.begin(); it != data.end(); ++it) {
        int key = it->first;
        const vector<vector<T>>& mat = it->second;

        cout << "Key = " << key << endl;

        print_nested_vector<T>(mat);
    }
}


// // Build the ordered list of pairs for a given maximum allowed positive value.
// // Allowed values are always from -1 up to max_allowed.
// // Used to get the state index for haplotype-specific copy number changes in the rate matrices of multiple-level Markov model.
// inline vector<pair<int,int>> build_pair_states(int max_allowed) {
//     if (max_allowed < 1) {
//         throw std::invalid_argument("max_allowed must be >= 1");
//     }

//     // Base order for max_allowed = 1
//     vector<pair<int,int>> states = {
//         {-1,-1}, {-1,0}, {0,-1}, {-1,1}, {1,-1},
//         {0,0}, {0,1}, {1,0}, {1,1}
//     };

//     // Extend from 2, 3, ..., max_allowed
//     for (int k = 2; k <= max_allowed; ++k) {
//         vector<pair<int,int>> new_states;

//         // Add all new states involving k:
//         // (k, x) and (x, k), where x runs from -1 to k
//         // avoiding duplicate (k, k)
//         for (int x = -1; x <= k; ++x) {
//             if (x < k) {
//                 new_states.push_back({x, k});
//                 new_states.push_back({k, x});
//             } else {
//                 new_states.push_back({k, k});
//             }
//         }

//         // Sort new states by:
//         // 1) increasing sum a+b
//         // 2) increasing first component a
//         std::sort(new_states.begin(), new_states.end(),
//                   [](const pair<int,int>& p1, const pair<int,int>& p2) {
//                       int s1 = p1.first + p1.second;
//                       int s2 = p2.first + p2.second;
//                       if (s1 != s2) return s1 < s2;
//                       return p1.first < p2.first;
//                   });

//         // Merge old states and new states by sum,
//         // keeping old states first when sums are equal
//         vector<pair<int,int>> merged;
//         merged.reserve(states.size() + new_states.size());

//         size_t i = 0, j = 0;
//         while (i < states.size() && j < new_states.size()) {
//             int sum_old = states[i].first + states[i].second;
//             int sum_new = new_states[j].first + new_states[j].second;

//             if (sum_old <= sum_new) {
//                 merged.push_back(states[i++]);   // old states first on tie
//             } else {
//                 merged.push_back(new_states[j++]);
//             }
//         }

//         while (i < states.size()) merged.push_back(states[i++]);
//         while (j < new_states.size()) merged.push_back(new_states[j++]);

//         states.swap(merged);
//     }

//     return states;
// }


// states for abosulute haplotype-specific copy number changes: -1/-1   -1/0    0/-1    -1/1    0/0 
inline vector<pair<int,int>> build_pair_states(int max_allowed) {
    if (max_allowed < 1) {
        throw std::invalid_argument("max_allowed must be >= 1");
    }

    vector<pair<int,int>> states;
    int max_hap_cn = NORM_PLOIDY / 2 + max_allowed;  // maximum copy number on one haplotype, used to determine the state space

     for (int total = 0; total <= 2 * max_hap_cn; ++total) {
        for (int a = 0; a <= total; ++a) {
            int b = total - a;

            if (a >= 0 && a <= max_hap_cn &&
                b >= 0 && b <= max_hap_cn) {
                states.push_back({a - NORM_PLOIDY / 2, b - NORM_PLOIDY / 2});
            }
        }
    }

    return states;    
}


// states for absolute haplotype-specific copy numbers: 0/0     0/1     1/0     0/2     1/1
inline vector<pair<int,int>> build_pair_states_cn(int max_allowed) {
    if (max_allowed < 1) {
        throw std::invalid_argument("max_allowed must be >= 1");
    }

    vector<pair<int,int>> states;
    int max_hap_cn = NORM_PLOIDY / 2 + max_allowed;  // maximum copy number on one haplotype, used to determine the state space

     for (int total = 0; total <= 2 * max_hap_cn; ++total) {
        for (int a = 0; a <= total; ++a) {
            int b = total - a;

            if (a >= 0 && a <= max_hap_cn &&
                b >= 0 && b <= max_hap_cn) {
                states.push_back({a, b});
            }
        }
    }

    return states;
}

// Return the index of (a,b) in the generated state list.
// Check the index boundaries to obtain correct state indices and reset when needed
// Returns -1 if not found.
inline int get_pair_index(int a, int b, int max_change_haplotype, const vector<pair<int,int>>& states) {
    if(a < MIN_CHANGE_HAPLOTYPE) a = MIN_CHANGE_HAPLOTYPE;
    if(a > max_change_haplotype) a = max_change_haplotype;
    if(b < MIN_CHANGE_HAPLOTYPE) b = MIN_CHANGE_HAPLOTYPE;
    if(b > max_change_haplotype) b = max_change_haplotype;  

    for (size_t i = 0; i < states.size(); ++i) {
        if (states[i].first == a && states[i].second == b) {
            return static_cast<int>(i);
        }
    }
    return -1;
}


inline pair<int,int> get_pair_from_index(int index, const vector<pair<int,int>>& states) {
    if (index < 0 || index >= static_cast<int>(states.size())) {
        throw std::out_of_range("Pair index out of range");
    }

    return states[index];
}

// Optional helper: build lookup map for faster repeated queries
inline map<pair<int,int>, int> build_pair_index_map(const vector<pair<int,int>>& states) {
    map<pair<int,int>, int> idx;
    for (size_t i = 0; i < states.size(); ++i) {
        idx[states[i]] = static_cast<int>(i);
    }
    return idx;
}


inline int encode_change_pair(int a, int b, int max_allowed) {
    static unordered_map<int, vector<std::pair<int,int>>> cache;
    if (!cache.count(max_allowed)) {
        cache[max_allowed] = build_pair_states(max_allowed);
    }
    const auto& states = cache[max_allowed];

    for (size_t i = 0; i < states.size(); ++i) {
        if (states[i].first == a && states[i].second == b) {
            return static_cast<int>(i);
        }
    }

    throw runtime_error("Unsupported haplotype-specific change pair: " +
                        std::to_string(a) + "/" + std::to_string(b));
}


inline void print_pair_states(const vector<pair<int,int>>& states) {
    cout << "Haplotype-specific change states:" << endl;
    for (size_t i = 0; i < states.size(); ++i) {
        cout << "Index " << i << ": " << states[i].first << "/" << states[i].second << endl;
    }
}

inline int max_abs_change(const vector<int>& vals) {
    int m = 0;
    for (int x : vals) m = max(m, std::abs(x));
    return m;
}


inline double mean_int_vector(const vector<int>& v) {
    if (v.empty()) {
        throw runtime_error("Cannot compute mean of empty vector.");
    }
    return static_cast<double>(accumulate(v.begin(), v.end(), 0)) / v.size();
}

// Function to round a floating point value to a specified precision
// https://www.geeksforgeeks.org/cpp/floating-point-values-as-keys-in-std-map/
inline double round_to_precision(double value, int precision){
    // Calculate the factor to scale the value
    double factor = pow(10, precision);
    // Round the value to the nearest integer after scaling,
    // then scale it back
    return round(value * factor) / factor;
}


#endif
