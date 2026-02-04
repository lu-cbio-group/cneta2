#ifndef LIKELIHOOD_HPP
#define LIKELIHOOD_HPP

//
// Likelihood computation
//

#include "parse_cn.hpp"
#include "evo_tree.hpp"
#include "tree_op.hpp"
#include "model.hpp"


// using namespace std;


// from book Effective STL
struct DeleteObject{
  template<typename T>
  void operator()(const T* ptr) const{
    delete [] ptr;
  }
};


// information derived from input data for DECOMP model
struct OBS_DECOMP{
  int m_max;   // maximum copy of a segment before chr-level events, used in likelihood table initialization

  // used to get dimension of Pmatrix at different levels
  int max_wgd;
  int max_chr_change;
  int max_site_change;

  // used in initialize_lnl_table_decomp, also defined in INPUT_DATA
  vector<int> sample_num_wgd;  // possible number of WGD events, used in likelihood table initialization
  vector<vector<int>> sample_change_chr; // possible number of chr-level events, used in likelihood table initialization
  // vector<int> sample_max_cn;  // not used in likelihood computation

  void print() const {
    cout << "OBS_DECOMP:" << endl;
    cout << "\tm_max: " << m_max << endl;
    cout << "\tmax_wgd: " << max_wgd << endl;
    cout << "\tmax_chr_change: " << max_chr_change << endl;
    cout << "\tmax_site_change: " << max_site_change << endl;

    // cout << "\tsample_num_wgd: ";
    // for (int v : sample_num_wgd) {
    //     cout << v << " ";
    // }
    // cout << endl;

    // cout << "\tsample_change_chr:" << endl;
    // for (size_t i = 0; i < sample_change_chr.size(); ++i) {
    //     cout << "\t\t[" << i << "]: ";
    //     for (int v : sample_change_chr[i]) {
    //         cout << v << " ";
    //     }
    //     cout << endl;
    // }
  }
};


// TODO: merged duplicated data structures
// information derived from input data, used in ancestral state reconstruction
struct MAX_DECOMP{
  int m_max;   // maximum copy of a segment before chr-level events, used in likelihood table initialization

  // used to get dimension of Pmatrix at different levels
  int max_wgd;
  int max_chr_change;
  int max_site_change;
};


// information derived from input data, not used for now
struct MAX_CHANGE{
  // used to get dimension of Pmatrix at different levels
  int max_wgd;
  int max_chr_change;
  int max_site_change;

  // TODO: may be extended in the future to distinguish gain/loss and arm-level changes

  void print() const {
      cout << "MAX_CHANGE:" << endl;
      cout << "\tMax WGD: " << max_wgd << endl;
      cout << "\tMax chromosome change: " << max_chr_change << endl;
      cout << "\tMax site/segment change: " << max_site_change << endl;
  }
};


// dimensions of rate matrices at different levels of copy number changes
struct DIM_DECOMP{
  int dim_wgd;
  int dim_chr;
  int dim_seg;

  void print() const {
      cout << "DIM_DECOMP:" << endl;
      cout << "\tWGD states: " << dim_wgd << endl;
      cout << "\tChromosome states: " << dim_chr << endl;
      cout << "\tSegment/site states: " << dim_seg << endl;
  }  
};


struct QMAT_DECOMP{
  double* qmat_wgd;
  double* qmat_chr;
  double* qmat_seg;

  // Ensure pointers are initialised
  QMAT_DECOMP() noexcept
      : qmat_wgd(nullptr),
        qmat_chr(nullptr),
        qmat_seg(nullptr)
  {}

  // forbid copying (owning raw memory)
  QMAT_DECOMP(const QMAT_DECOMP&) = delete;
  QMAT_DECOMP& operator=(const QMAT_DECOMP&) = delete;

  // move constructor
  QMAT_DECOMP(QMAT_DECOMP&& other) noexcept
      : qmat_wgd(other.qmat_wgd),
        qmat_chr(other.qmat_chr),
        qmat_seg(other.qmat_seg)
  {
      other.qmat_wgd = nullptr;
      other.qmat_chr = nullptr;
      other.qmat_seg = nullptr;
  }

  // move assignment
  QMAT_DECOMP& operator=(QMAT_DECOMP&& other) noexcept {
      if (this != &other) {
          delete [] qmat_wgd;
          delete [] qmat_chr;
          delete [] qmat_seg;

          qmat_wgd = other.qmat_wgd;
          qmat_chr = other.qmat_chr;
          qmat_seg = other.qmat_seg;

          other.qmat_wgd = nullptr;
          other.qmat_chr = nullptr;
          other.qmat_seg = nullptr;
      }
      return *this;
  }

  void free_all() noexcept {
      delete [] qmat_wgd;
      delete [] qmat_chr;
      delete [] qmat_seg;    
  }

  ~QMAT_DECOMP() {
    free_all();
  }
};


struct PMAT_DECOMP{
  map<double, double*> pmats_wgd;
  map<double, double*> pmats_chr;
  map<double, double*> pmats_seg;

  PMAT_DECOMP() = default;
  // forbid copying (owning raw memory)
  PMAT_DECOMP(const PMAT_DECOMP&) = delete;
  PMAT_DECOMP& operator=(const PMAT_DECOMP&) = delete;

  // move constructor (transfer ownership)
  PMAT_DECOMP(PMAT_DECOMP&& other) noexcept
      : pmats_wgd(std::move(other.pmats_wgd)),
        pmats_chr(std::move(other.pmats_chr)),
        pmats_seg(std::move(other.pmats_seg))
  {
      // other maps are left in a valid but unspecified state (usually empty)
      // we do NOT delete pointers in 'other' because we stole them
      other.pmats_wgd.clear();
      other.pmats_chr.clear();
      other.pmats_seg.clear();
  }

  // delete all owned pointers and clear maps
  void free_all() noexcept {
      auto free_map = [](std::map<double, double*>& m) noexcept {
          for (auto& kv : m) {
              delete [] kv.second;   // assumes allocation with new[]
              kv.second = nullptr;
          }
          m.clear();
      };

      free_map(pmats_wgd);
      free_map(pmats_chr);
      free_map(pmats_seg);
  }

  // move assignment
  PMAT_DECOMP& operator=(PMAT_DECOMP&& other) noexcept {
      if (this != &other) {
          free_all();  // free currently owned memory

          pmats_wgd = std::move(other.pmats_wgd);
          pmats_chr = std::move(other.pmats_chr);
          pmats_seg = std::move(other.pmats_seg);

          other.pmats_wgd.clear();
          other.pmats_chr.clear();
          other.pmats_seg.clear();
      }
      return *this;
  }


  ~PMAT_DECOMP() {
      free_all();
  }
};


// used in computing likelihood of children nodes
struct PROB_DECOMP{
  double* pbli_wgd;
  double* pblj_wgd;

  double* pbli_chr;
  double* pblj_chr;

  double* pbli_seg;
  double* pblj_seg;

  
  PROB_DECOMP() = default;
  PROB_DECOMP(const PROB_DECOMP&) = delete;
};

// separated for function refactoring
struct PROB_DECOMP1{
  double* pbl_wgd;
  double* pbl_chr;
  double* pbl_seg;


  PROB_DECOMP1() = default;
  PROB_DECOMP1(const PROB_DECOMP1&) = delete;
};


// variables used in likelihood computation
struct LNL_TYPE{
  int model;
  int cn_max;
  int is_total; // whether or not the input is total copy number

  // used to check tree validity
  int cons;
  double max_tobs;  // maximum sampling time
  int patient_age;   // used in BFGS constraints

  int use_repeat;   // whether or not to use repeated site patterns, used in get_likelihood_chr*
  int correct_bias; // Whether or not to correct acquisition bias, used in get_likelihood_*
  int num_invar_bins; // used in get_likelihood_*

  int cn_type; // Whether or not to only consider segment-level mutations, used in get_likelihood_revised

  int infer_wgd; // whether or not to infer WGD status of a sample, called in initialize_lnl_table_decomp
  int infer_chr; // whether or not to infer chromosome gain/loss status of a sample, called in initialize_lnl_table_decomp

  vector<int> knodes;
};


struct LNL_TABLE{
  vector<vector<double>> lnl_table_wgd;  // likelihood table for one site
  vector<vector<double>> lnl_table_chr;
  vector<vector<double>> lnl_table_seg; 
};


struct LNL_VAL{
  double lnl_wgd;
  double lnl_chr;
  double lnl_seg;
};


const double LARGE_LNL = -1e9;
const double SMALL_LNL = -1e20;
const double MAX_NLNL = 1e20;


/****************** common functions *******************/
void print_lnl_at_tips(const evo_tree& rtree, const vector<int>& obs, const vector<vector<double>>& L_sk_k, int nstate);


void print_tree_lnl(const evo_tree& rtree, vector<vector<double>>& L_sk_k, int nstate);

// From https://stackoverflow.com/questions/17074324/how-can-i-sort-two-vectors-in-the-same-way-with-criteria-that-uses-only-one-of
template <typename T, typename Compare>
std::vector<std::size_t> sort_permutation(
    const std::vector<T>& vec,
    const Compare& compare)
{
    std::vector<std::size_t> p(vec.size());
    std::iota(p.begin(), p.end(), 0);
    std::sort(p.begin(), p.end(),
        [&](std::size_t i, std::size_t j){ return compare(vec[i], vec[j]); });
    return p;
}

template <typename T>
std::vector<T> apply_permutation(
    const std::vector<T>& vec,
    const std::vector<std::size_t>& p)
{
    std::vector<T> sorted_vec(vec.size());
    std::transform(p.begin(), p.end(), sorted_vec.begin(),
        [&](std::size_t i){ return vec[i]; });
    return sorted_vec;
}



/****************** functions for non DECOMP model *******************/
// Create likelihood vectors at the tip node, one table for each site
// obs: a vector of CNs for one site
// L_sk_k has one row for each tree node and one column for each possible state
// only used in get_likelihood_chr
// extracted as a function to avoid duplication in selection statement
void initialize_lnl_table(vector<vector<double>>& L_sk_k, const vector<int>& obs, const evo_tree& rtree, int model, int nstate, int is_total);


/** 
 * @brief Get the probability of children nodes given the state at the parent node 
 * only used in get_likelihood_site
 * extracted as a function to avoid duplication in selection statement 
 * @param L_sk_k likelihood table
 * @param rtree the evolutionary tree
 * @param pbli transition probability matrix for branch i from parent state si to child state sk, a matrix of nstate * nstate, stored as a 1D array
 * @param pblj transition probability matrix for branch j from parent state sj to child state sk
 * @param nsk state at the parent node 
 * @param ni the first child node id
 * @param nj the second child node id
 * @param bli branch length for branch i
 * @param blj branch length for branch j
 * @param model evolutionary model
 * @param nstate number of possible states, also the dimension of transition probability matrix
 * nstate = cn_max + 1 when the input is total copy number under BOUNDT model
 * @return the probability of children nodes given the state at the parent node
 */
inline double get_prob_children(vector<vector<double>>& L_sk_k, const evo_tree& rtree, double* pbli, double* pblj, int nsk, int ni, int nj, int bli, int blj, int model, int nstate){
    double Li = 0.0;
    for(int si = 0; si < nstate; ++si){   // sum over all possible states at the parent node
      if(L_sk_k[ni][si] > 0){
        if(model == MK){
          Li += get_transition_prob(rtree.mu, bli, nsk, si) * L_sk_k[ni][si];
        }else{
          // pbli[nsk + si * nstate] is the transition probability from state nsk to si along branch of length bli
          Li += pbli[nsk + si * nstate] * L_sk_k[ni][si];
        }
      }
    }

    double Lj = 0.0;
    for(int sj = 0; sj < nstate; ++sj){
      if(L_sk_k[nj][sj] > 0){
        if(model == MK){
          Lj += get_transition_prob(rtree.mu, blj, nsk, sj) * L_sk_k[nj][sj];
        }else{
          Lj += pblj[nsk + sj * nstate] * L_sk_k[nj][sj];
        }
      }
    }

    return Li * Lj;
}



// Get the likelihood on one site of a chromosome (assuming higher level events on nodes)
// z: possible changes in copy number caused by chromosome gain/loss
// only used in get_likelihood_chr
// extracted as a function to avoid duplication in selection statement
void get_likelihood_site(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const vector<int>& knodes, const vector<double>& blens, const vector<double*>& pmat_per_blen, const int& has_wgd, const int& z, const int& model, const int& nstate);


// Get the likelihood on a set of chromosomes
// only used in get_likelihood_revised
// extracted as a function to avoid duplication in selection statement
double get_likelihood_chr(const map<int, vector<vector<int>>>& vobs, const evo_tree& rtree, const vector<int>& knodes, const vector<double>& blens, const vector<double*>& pmat_per_blen, const int& has_wgd, const int& cn_type, const int& use_repeat, const int& model, const int& nstate, const int& is_total);



// Incorporate chromosome gain/loss and WGD
// Model 2: Treat total copy number as the observed data and the haplotype-specific information is missing
// additional parameters moved to lnl_type to facilitate calling in optimization
double get_likelihood_revised(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, LNL_TYPE& lnl_type);


/* Compute the likelihood without grouping sites by chromosome, only considering site duplication/deletion (not used)
Precondition: the tree is valid
Ns: number of samples
Nchar: number of characters for each sample
vobs: the observed data matrix
rtree: the given tree
model: model of evolution
*/
double get_likelihood(const vector<vector<int>>& vobs, evo_tree& rtree, const vector<int>& knodes, int model, int cn_max, int is_total);


// Get the likelihood of the tree from likelihood table
double extract_tree_lnl(vector<vector<double>>& L_sk_k, int Ns, int model);


/************** functions for model DECOMP (copy number change) **************/
void initialize_lnl_table_change(LNL_TABLE& L_sk_k, const evo_tree& rtree, const vector<CN_CHANGE>& obs_change, const DIM_DECOMP& dim_decomp, const OBS_DECOMP& obs_decomp, int debug);

// not used for now due to pointer issues
void build_rate_matrices(QMAT_DECOMP& qmat_decomp, const evo_tree& rtree, const OBS_DECOMP& obs_decomp, const DIM_DECOMP& dim_decomp, int debug);
void build_transition_matrices(PMAT_DECOMP& pmat_decomp, const evo_tree& rtree, const vector<int>& knodes, const QMAT_DECOMP& qmat_decomp, const DIM_DECOMP& dim_decomp, const OBS_DECOMP& obs_decomp, int debug); 


LNL_VAL compute_child_likelihood(int node, const CN_CHANGE& sk, const LNL_TABLE& L_sk_k, const PROB_DECOMP1& prob_decomp, const OBS_DECOMP& obs_decomp, const DIM_DECOMP& dim_decomp, int debug);

LNL_VAL get_prob_children_change(LNL_TABLE& L_sk_k, const evo_tree& rtree, const CN_CHANGE& sk, PROB_DECOMP1& prob_decompi, PROB_DECOMP1& prob_decompj, const OBS_DECOMP& obs_decomp, const DIM_DECOMP& dim_decomp, int ni, int nj, int is_total, int debug);

void get_likelihood_site_change(LNL_TABLE& L_sk_k, const evo_tree& rtree, const vector<int>& knodes, const PMAT_DECOMP& pmat_decomp, const DIM_DECOMP& dim_decomp, const OBS_DECOMP& obs_decomp, int is_total, int debug);

double get_likelihood_chr_change(const evo_tree& rtree, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const vector<int>& knodes, const PMAT_DECOMP& pmat_decomp, const DIM_DECOMP& dim_decomp, const OBS_DECOMP& obs_decomp, const LNL_TYPE& lnl_type, int is_total, int debug);
double get_likelihood_change(evo_tree& rtree, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const LNL_TYPE& lnl_type, int debug);

double extract_tree_lnl_change(const LNL_TABLE& L_sk_k, const OBS_DECOMP& obs_decomp, int Ns, int debug);  




/************** functions for model DECOMP (copy number state decomposition) **************/

// L_sk_k has one row for each tree node and one column for each possible state; chr starting from 1
// This function is critical in obtaining correct likelihood. If one tip is not initialized, the final likelihood will be 0.
// nstate = comps.size();
vector<vector<double>> initialize_lnl_table_decomp(vector<int>& obs, const OBS_DECOMP& obs_decomp, int chr, const evo_tree& rtree, const set<vector<int>>& comps, int infer_wgd, int infer_chr, int cn_max, int is_total = 1);


// Assume the likelihood table is for each total copy number (no WGD order considered, deprecated)
double get_prob_children_decomp(vector<vector<double>>& L_sk_k, const evo_tree& rtree, map<int, set<vector<int>>>& decomp_table, int sk, int cn_max, int nstate, PROB_DECOMP& prob_decomp, DIM_DECOMP& dim_decomp, int ni, int nj, int bli, int blj, int is_total);


// Assume the likelihood table is for each combination of states
double get_prob_children_decomp2(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const set<vector<int>>& comps, int sk, PROB_DECOMP& prob_decomp, DIM_DECOMP& dim_decomp, int ni, int nj, int bli, int blj, int cn_max, int nstate, int is_total);


// Get the likelihood on one site of a chromosome
// Assuming each observed copy number is composed of three type of events.
// Sum over all possible states for initial and final nodes
// Allow at most one WGD event along a branch
void get_likelihood_site_decomp(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const set<vector<int>>& comps, const vector<int>& knodes, PMAT_DECOMP& pmat_decomp, DIM_DECOMP& dim_decomp, int cn_max, int is_total);

// Used when WGD is considered, dealing with mutations of different types at different levels (no WGD order considered, deprecated)
double get_likelihood_chr_decomp(const map<int, vector<vector<int>>>& vobs, const OBS_DECOMP& obs_decomp, const evo_tree& rtree, const set<vector<int>>& comps, const vector<int>& knodes, PMAT_DECOMP& pmat_decomp, DIM_DECOMP& dim_decomp, int infer_wgd, int infer_chr, int use_repeat, int cn_max, int is_total);

// Computing likelihood when WGD and chr gain/loss are incorporated
// Assume likelihood is for haplotype-specific information
double get_likelihood_decomp(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, int debug = 0);

// Get the likelihood of the tree from likelihood table of state combinations
double extract_tree_lnl_decomp(vector<vector<double>>& L_sk_k, const set<vector<int>>& comps, int Ns);

#endif
