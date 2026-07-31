#ifndef OPTIMIZATION_HPP
#define OPTIMIZATION_HPP

//
// Maximum likelihood and optimisation functions
//

#include <gsl/gsl_multimin.h>

#include "lbfgsb/lbfgsb_new.h"

#include "common.hpp"
#include "likelihood.hpp"
// #include "nni.hpp"

// using namespace std;

const double ERROR_X = 1.0e-4;
// The minimum mutation rates allowed
const double MIN_MRATE = 1.0e-20;
// The maximum mutation rates allowed
const double MAX_MRATE = 1;

const double M_MIN = 1.0e-3;  // minimum multiplier for mutation rates
const double M_MAX = 100;   // maximum multiplier for mutation rates

// Bounds for mode-3 RLC multipliers (not absolute mutation rates — rates can be multiplied up or down)
const double RLC_MULT_MIN = 1e-3;
const double RLC_MULT_MAX = 10;


// [2026-07-14 moved] was previously declared further down, right before
// get_bsr_rate_slots(); moved up here so OPT_TYPE (below) can cache a
// vector<BsrRateSlot> instead of every caller recomputing it. See
// get_bsr_rate_slots() later in this file for how slots are populated.
// One active CNA rate type for BSR optimisation.
// tree_rate: pointer-to-member on evo_tree (e.g. &evo_tree::dup_rate) — the global reference.
// edge_field: pointer-to-member on RateSet (e.g. &RateSet::dup) — the per-branch field.
// Pairing is enforced by construction in get_bsr_rate_slots, not by a raw double snapshot.
struct BsrRateSlot {
    const char*         name;
    double evo_tree::*  tree_rate;   // global reference rate on rtree
    double RateSet::*   edge_field;  // corresponding field in edge_rates
};

// bundle of variables used in optimization (values from input)
struct OPT_TYPE{
  int estmu;   // used to determine optimization function
  int bsr_mode;  // 0: constant rate, 1: shared multiplier per branch, 2: independent per-event rates, 3: random local clock
  double tolerance;
  int miter;
  int opt_one_branch;
  vector<double> tobs;  // sample times, used to get constaints in optimization
  double scale_tobs;   // scale factor to get lower limit of root age when doing constrained optimization (BFGS) based on maximimum sample time difference
  // bsr_mode=3 (ML-RLC) fields
  vector<int> rlc_shift_eids;  // current δ=1 shift edges selected by outer search
  double rlc_raw_logL = std::numeric_limits<double>::quiet_NaN();
  double rlc_penalized_score = std::numeric_limits<double>::quiet_NaN();
  // which outer search drives shift-edge selection under bsr_mode=3:
  // 0: exhaustive search, 1: stepwise greedy search (default), 2: genetic algorithm search
  int rlc_search_method = 1;
  // [2026-07-21 added] model-selection criterion driving shift-edge search comparisons
  // (see compute_rlc_ic): 0 = AIC (default), 1 = BIC.
  int rlc_criterion = 0;
  // [2026-07-14 added] cache of get_bsr_rate_slots(rtree, cn_type), refreshed once per
  // max_likelihood_BFGS call (bsr_mode>0) instead of being recomputed on every BFGS
  // objective/gradient evaluation by update_variables_transformed/
  // update_edge_rates_rlc_from_x — those used to call get_bsr_rate_slots() themselves
  // every call, which is wasted work since it's invariant for the whole optimization
  // (the 5 global reference rates are frozen for bsr_mode>0). See docs/flowcharts.md.
  vector<BsrRateSlot> bsr_slots_cache;
  // [2026-07-14 added] set to true once optimize_tree_by_bsr_mode's auto-calibration
  // sub-optimization has *successfully* estimated the 5 global reference rates (calib_nlnl
  // came back finite, not MAX_NLNL). maximize_tree_likelihood's outer retry loop calls
  // optimize_tree_by_bsr_mode again if an attempt fails, which previously re-ran the whole
  // calibration sub-optimization from scratch every retry even when calibration itself had
  // already succeeded and only the later shift-edge search failed. Checked (and left false
  // on calibration failure, so a retry still retries calibration too) instead of assuming
  // "called before" means "succeeded before". See docs/flowcharts.md.
  bool bsr_calibrated = false;
};

struct GSL_PARAM{
  evo_tree rtree;
  map<int, vector<vector<int>>> vobs;
  LNL_TYPE lnl_type;
};

/*****************************************************
    optimization with GSL simplex methods (deprecated)
*****************************************************/

double my_f(const gsl_vector *v, void *params);
double my_f_mu(const gsl_vector *v, void *params);
double my_f_cons(const gsl_vector *v, void *params);
double my_f_cons_mu(const gsl_vector *v, void *params);


// Output error information without aborting
inline void my_err_handler(const char * reason, const char * file, int line, int gsl_errno){
    fprintf (stderr, "failed, %s, reason: %s, file %s line %d \n", gsl_strerror(gsl_errno), reason, file, line);
    return;
}

// given a tree, maximise the branch lengths (and optionally mu) assuming branch lengths are independent or constrained in time
// use GSL simplex optimization
void max_likelihood(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double& min_nlnl, const double ssize);


/*****************************************************
    One dimensional optimization with Brent method
*****************************************************/

// Compute the likelihood of the tree with one parameter (one branch or one mutation rate)
// Used in Brent optimization
// return negative likelihood function for minimalization
double computeFunction(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, double value, int type = -1);



/* Brents method in one dimension */
/**
    This function calculate f(value) of the f() function, used by other general optimization method to minimize it.
    Please always override this function to adapt to likelihood or parsimony score.
    The default is for function f(x)=x.
    @param value x-value of the function
    @return f(value) of function f you want to minimize
*/
double brent_opt(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, int type, double ax, double bx, double cx, double tol,
                          double *foptx, double *f2optx, double fax, double fbx, double fcx);


/**
    the brent method to find the value that minimizes the computeFunction().
    @return the x-value that minimize the function
    @param xmin lower bound
    @param xmax upper bound
    @param xguess first guess
    @param tolerance tolerance
    @param fx (OUT) function value at the minimum x found
    @param ferror (OUT) Dont know
*/
double minimizeOneDimen(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, int type, double xmin, double xguess, double xmax, double tolerance, double *fx, double *f2x);


// Using Brent method to optimize the likelihood of mutation rate
double optimize_mutation_rates(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, double tolerance);

// Optimizing the branch length of (node1, node2) with Brent method
// Optimize mutation rates if necessary
double optimize_one_branch(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, double tolerance, Node* node1, Node* node2);


// used in optimize_all_branches
void compute_best_traversal(evo_tree& rtree, NodeVector &nodes, NodeVector &nodes2);


// optimize each branch
double optimize_all_branches(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, int my_iterations, double tolerance);



/*****************************************************
    L-BFGS-B method
*****************************************************/

// The number of parameters to estimate, different when the mutation rates are estimated
// n_bsr_edges: number of edges receiving a multiplier (active_eids.size()).
// Defaults to -1, which falls back to nparams_est (correct for cons=0; wrong for cons=1 — see TODO).
inline int get_ndim(int estmu, int nparams_est, int model, int cn_type, int nrates = 0, int n_types_per_edge = 0, int n_bsr_edges = -1){
    int ndim = 0;

    // variable rate: estimate n_types_per_edge multipliers per active branch
    // bsr_mode=1: n_types_per_edge=1, one shared multiplier per branch
    // bsr_mode=2: n_types_per_edge = number of active rate types for cn_type (1-5), one multiplier per type per branch
    // bsr_mode=3: 
    if(n_types_per_edge > 0){
        int n_edges = (n_bsr_edges >= 0) ? n_bsr_edges : nparams_est;
        // + bsr_slots.size() for global rates
        return nparams_est + n_types_per_edge * n_edges + nrates;  // +1 for the global rate (estmu=1) if applicable
    }

    // bsr_mode=0: n_types_per_edge=0, constant rate, handled by estmu logic below
    if(!estmu){
      ndim = nparams_est;
    }else{
        if(model == MK){
            ndim = nparams_est + 1;
        }else{
            switch(cn_type){
                case ALL:{
                    ndim = nparams_est + 5;
                    break;
                }
                case ONLY_SEG:{
                    ndim = nparams_est + 2;
                    break;
                }
                case EXCLUDE_SEG:{
                    ndim = nparams_est + 3;         
                    break;
                }
                case EXCLUDE_CHR:{
                    ndim = nparams_est + 3;
                    break;
                }
                case EXCLUDE_WGD:{
                    ndim = nparams_est + 4;
                    break;
                }
                default:
                    cerr << "### ERROR: Unknown copy number alteration type for optimization!" << endl;
                    exit(1);
            }
        }
    }
    return ndim;
}


// Update the tree after each iteration in the BFGS optimization (deprecated)
void update_variables(evo_tree& rtree, int model, int cons, int estmu, double *x);


// Update the tree after each iteration in the BFGS optimization
// Estimate time intervals rather than branch length in order to avoid negative terminal branch lengths
// Sort node times in increasing order and take the first Ns intervals
void update_variables_transformed(evo_tree &rtree, double *x, LNL_TYPE &lnl_type, OPT_TYPE &opt_type);

// [2026-07-14] bsr_mode param added so this can assert it's never called for bsr_mode>0
// (global reference rates must stay frozen there) — see the assert in the definition
// and docs/flowcharts.md for why this guard exists.
void update_tree_rates(LNL_TYPE &lnl_type, evo_tree &rtree, double *x, int nparams_est, int bsr_mode);

/**
    the target function which needs to be optimized
    @param x the input vector x
    @return the function value at x
*/
double targetFunk(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, double x[]);


/**
	the approximated derivative function
	@param x the input vector x
	@param dfx the derivative at x
	@return the function value at x
*/
double derivativeFunk(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int ndim, double x[], double dfx[]);

double optimFunc(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int nvar, double *vars);

double optimGradient(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int nvar, double *x, double *dfx);

void lbfgsb(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int n, int m, double *x, double *l, double *u, int *nbd,
		double *Fmin, int *fail,
		double factr, double pgtol,
		int *fncount, int *grcount, int maxit, char *msg,
		int trace, int nREPORT);


/**
 Function to access the L-BFGS-B function, taken from IQ-TREE package which is further taken from HAL_HAS software package
 1. int nvar : The number of the variables
 2. double* vars : initial values of the variables
 3. double* lower : lower bounds of the variables
 4. double* upper : upper bounds of the variables
 5. double pgtol: gradient tolerance
 5. int maxit : max # of iterations
 @return minimized function value
 After the function is invoked, the values of x will be updated
*/
double L_BFGS_B(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, int n, double* x, double* l, double* u);

// [2026-07-14] BsrRateSlot moved up next to OPT_TYPE (see near the top of this file) so
// OPT_TYPE can cache a vector<BsrRateSlot>.


// Returns the active rate slots for the given cn_type.
// bsr_mode=1: use slots to know which rates are active, share one multiplier across all.
// bsr_mode=2: each slot gets its own independent multiplier.
// [2026-07-14 added] Shared index formula for bsr_mode=2/3's per-(edge,type) BFGS
// variable layout: active-edge k's rate-type slot t lives at this offset into
// x[]/variables[]. Previously hand-written identically in three places (the
// population loop in max_likelihood_BFGS, the bsr_mode==2 read-back in
// update_variables_transformed, and update_edge_rates_rlc_from_x for bsr_mode==3),
// risking silent drift if only some were updated. See docs/flowcharts.md.
inline int bsr_var_index(int nparams_est, int k, int n_types_per_edge, int t){
    return nparams_est + k * n_types_per_edge + t + 1;
}

inline vector<BsrRateSlot> get_bsr_rate_slots(const evo_tree& rtree, int cn_type){
    vector<BsrRateSlot> slots;
    auto add = [&](const char* name, double evo_tree::* tr, double RateSet::* ef){
        if(rtree.*tr > 0) slots.push_back({name, tr, ef});
    };
    switch(cn_type){
        case ONLY_SEG:
            add("dup", &evo_tree::dup_rate, &RateSet::dup);
            add("del", &evo_tree::del_rate, &RateSet::del);
            break;
        case EXCLUDE_SEG:
            add("chr_gain", &evo_tree::chr_gain_rate, &RateSet::chr_gain);
            add("chr_loss", &evo_tree::chr_loss_rate, &RateSet::chr_loss);
            add("wgd",      &evo_tree::wgd_rate,      &RateSet::wgd);
            break;
        case EXCLUDE_CHR:
            add("dup", &evo_tree::dup_rate, &RateSet::dup);
            add("del", &evo_tree::del_rate, &RateSet::del);
            add("wgd", &evo_tree::wgd_rate, &RateSet::wgd);
            break;
        case EXCLUDE_WGD:
            add("dup",      &evo_tree::dup_rate,      &RateSet::dup);
            add("del",      &evo_tree::del_rate,      &RateSet::del);
            add("chr_gain", &evo_tree::chr_gain_rate, &RateSet::chr_gain);
            add("chr_loss", &evo_tree::chr_loss_rate, &RateSet::chr_loss);
            break;
        default: // ALL
            add("dup",      &evo_tree::dup_rate,      &RateSet::dup);
            add("del",      &evo_tree::del_rate,      &RateSet::del);
            add("chr_gain", &evo_tree::chr_gain_rate, &RateSet::chr_gain);
            add("chr_loss", &evo_tree::chr_loss_rate, &RateSet::chr_loss);
            add("wgd",      &evo_tree::wgd_rate,      &RateSet::wgd);
    }
    return slots;
}

// Returns edge IDs that receive per-branch rate multipliers (excludes the normal-sample edge).
inline vector<int> get_active_bsr_eids(const evo_tree& rtree){
    vector<int> eids;
    for(int eid = 0; eid < (int)rtree.edges.size(); ++eid){
        // [2026-07-14 disabled] see evo_tree::is_normal_sample_edge in evo_tree.hpp
        // if(rtree.edges[eid].start == rtree.nleaf - 1 ||
        //    rtree.edges[eid].end   == rtree.nleaf - 1) continue;
        if(rtree.is_normal_sample_edge(eid)) continue;
        eids.push_back(eid);
    }
    return eids;
}

// Using BFGS method to get the maximum likelihood with lower and upper bounds (minimalize negative likelihood function)
// Note: the topology of rtree is fixed, yet its branch lengths may be updated in the optimization process
void max_likelihood_BFGS(evo_tree &rtree, const map<int, vector<vector<int>>> &vobs, const map<int, vector<vector<CN_CHANGE>>> &vobs_change, const OBS_DECOMP &obs_decomp, const set<vector<int>> &comps, LNL_TYPE &lnl_type, OPT_TYPE &opt_type, double &min_nlnl, int debug = 0);

// [2026-07-14] bsr_mode param added so this can assert it's never called for bsr_mode>0
// (global reference rates must stay frozen there) — see the assert in the definition
// and docs/flowcharts.md for why this guard exists.
void set_tree_rates(int cn_type, evo_tree &rtree, int nparams_est, double* variables, double* lower_bound, double* upper_bound, int bsr_mode);

/****** bsr_mode=3 ML-RLC functions ******/

// Propagate local rates from root to tips given shift edges and their per-type multipliers.
// Non-shift edges inherit parent's local RateSet exactly.
// Shift edges scale parent's local RateSet element-wise by multipliers[k]
// (one independent multiplier per rate type: dup, del, chr_gain, chr_loss, wgd — same
// per-edge structure as bsr_mode=2, just applied sparsely with inheritance).
// Writes rtree.edge_rates[eid] for every edge.
void update_edge_rates_rlc(
    evo_tree& rtree,
    const vector<int>& shift_eids,
    const vector<RateSet>& multipliers
);

// Convenience wrapper: reads per-type multipliers from optimizer vector x[].
// For shift_eids[k], the multiplier for active rate slot t is at
// x[nparams_est + k*n_types_per_edge + t + 1], matching the layout written in
// max_likelihood_BFGS's bsr_mode==2||3 branch. bsr_slots (order/count) should be
// max_likelihood_BFGS's cached opt_type.bsr_slots_cache — [2026-07-14] changed from
// taking cn_type and recomputing get_bsr_rate_slots() on every call.
void update_edge_rates_rlc_from_x(
    evo_tree& rtree,
    const vector<int>& shift_eids,
    double* x,
    int nparams_est,
    const vector<BsrRateSlot>& bsr_slots
);

// [2026-07-21 added] AIC/BIC model-selection values for ML-RLC (bsr_mode=3), replacing
// rlc_penalized_score above. K*n_types_per_edge is the number of free parameters unlocked
// by the accepted shift edges (see rlc_penalized_score's original comment for why
// n_types_per_edge must multiply K). n_sites (lnl_type.n_sites_for_ic) is the effective
// number of independent observations backing logL, used as BIC's sample size. `score` is
// "higher is better" (negated IC) for whichever of AIC/BIC opt_type.rlc_criterion selects,
// so existing comparison code (pick_best_improving, sorting) keeps working unchanged; aic
// and bic are both always computed so callers can log/report the non-selected criterion
// too (search decisions still follow `score`/rlc_criterion only, not the other one).
struct RlcIC {
    double raw_logL;
    double aic;
    double bic;
    double score;
};

inline RlcIC compute_rlc_ic(double nlnl, int K, int n_types_per_edge, int n_sites, int criterion){
    RlcIC ic;
    ic.raw_logL = -nlnl;
    double k_params = (double)K * (double)n_types_per_edge;
    ic.aic = -2.0 * ic.raw_logL + 2.0 * k_params;
    ic.bic = -2.0 * ic.raw_logL + log((double)n_sites) * k_params;
    ic.score = -(criterion == 0 ? ic.aic : ic.bic);
    return ic;
}

// [2026-07-14 added] Scans indices [0, n), skipping any where evaluated[i] is false, and
// returns the index of the best-scoring one that beats cur_score+eps (ties broken by
// smaller id). Returns -1 if none qualify. Previously this ~25-line selection logic was
// hand-copied for stepwise_search_shift_edges's forward candidate selection and backward
// removal selection; centralized here so both (and the future genetic_search_shift_edges)
// share one implementation. See docs/flowcharts.md.
inline int pick_best_improving(int n, const vector<char>& evaluated, const vector<double>& scores,
                                const vector<int>& ids, double cur_score, double eps){
    int best = -1;
    for(int i = 0; i < n; ++i){
        if(!evaluated[i]) continue;
        double score = scores[i];
        if(score <= cur_score + eps) continue;
        if(best == -1){ best = i; continue; }
        bool better = score > scores[best] + eps;
        bool tie_but_smaller = fabs(score - scores[best]) <= eps && ids[i] < ids[best];
        if(better || tie_but_smaller) best = i;
    }
    return best;
}

// Outer stepwise greedy search for shift edges under ML-RLC.
// At each step adds/removes the shift edge that most improves compute_rlc_ic's score
// (AIC or BIC per opt_type.rlc_criterion).
// Updates rtree, opt_type.rlc_shift_eids, and min_nlnl in place.
void stepwise_search_shift_edges(
    evo_tree& rtree,
    const map<int, vector<vector<int>>>& vobs,
    const map<int, vector<vector<CN_CHANGE>>>& vobs_change,
    const OBS_DECOMP& obs_decomp,
    const set<vector<int>>& comps,
    LNL_TYPE& lnl_type,
    OPT_TYPE& opt_type,
    double& min_nlnl,
    int debug = 0
);

// Outer genetic-algorithm search for shift edges under ML-RLC (bsr_mode=3).
// Alternative to stepwise_search_shift_edges, selected via opt_type.rlc_search_method=2.
// Not yet implemented.
void genetic_search_shift_edges(
    gsl_rng* r,
    evo_tree& rtree,
    const map<int, vector<vector<int>>>& vobs,
    const map<int, vector<vector<CN_CHANGE>>>& vobs_change,
    const OBS_DECOMP& obs_decomp,
    const set<vector<int>>& comps,
    LNL_TYPE& lnl_type,
    OPT_TYPE& opt_type,
    double& min_nlnl,
    int debug = 0
);

// Brute-force search over all 2^ncand subsets of get_active_bsr_eids
// (ncand = candidate shift edges), used as a baseline to check what the best achievable
// model fit is. Only feasible for small ncand (e.g. 9 candidates -> 512 subsets on a 5-sample tree); prints a
// runtime warning for large ncand instead of refusing to run. Reports both the AIC-best and
// BIC-best subsets found (see compute_rlc_ic) since exhaustive enumeration makes both free
// from one run, but writes back rtree/opt_type using only the subset picked by
// opt_type.rlc_criterion. Selected via opt_type.rlc_search_method=0.
void exhaustive_search_shift_edges(
    evo_tree& rtree,
    const map<int, vector<vector<int>>>& vobs,
    const map<int, vector<vector<CN_CHANGE>>>& vobs_change,
    const OBS_DECOMP& obs_decomp,
    const set<vector<int>>& comps,
    LNL_TYPE& lnl_type,
    OPT_TYPE& opt_type,
    double& min_nlnl,
    int debug = 0
);

// Wrapper: dispatches to stepwise_search_shift_edges, genetic_search_shift_edges, or
// exhaustive_search_shift_edges for bsr_mode=3 (per opt_type.rlc_search_method), or
// max_likelihood_BFGS for all other modes.
void optimize_tree_by_bsr_mode(
    gsl_rng* r,
    evo_tree& rtree,
    const map<int, vector<vector<int>>>& vobs,
    const map<int, vector<vector<CN_CHANGE>>>& vobs_change,
    const OBS_DECOMP& obs_decomp,
    const set<vector<int>>& comps,
    LNL_TYPE& lnl_type,
    OPT_TYPE& opt_type,
    double& min_nlnl,
    int debug = 0
);

// Optimizing the branch length of (node1, node2) with BFGS to incorporate constraints imposed by patient age and tip timings.
// Optimize mutation rates if necessary
double optimize_one_branch_BFGS(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, OPT_TYPE& opt_type, Node* node1, Node* node2);



#endif
