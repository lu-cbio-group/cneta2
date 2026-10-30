#ifndef STATE_HPP
#define STATE_HPP


#include "likelihood.hpp"

// using namespace std;


struct LNL_TABLE{
    vector<vector<double>> lnl_table_wgd;
    vector<vector<double>> lnl_table_chr;
    vector<vector<double>> lnl_table_seg;

    LNL_TABLE() = default;
    LNL_TABLE(const LNL_TABLE&) = delete;  // disable copy constructor
};


struct STATE_TABLE{
    //  record copy number changes at multiple levels
    vector<vector<int>> state_table_wgd;
    vector<vector<int>> state_table_chr;
    vector<vector<int>> state_table_seg;

    STATE_TABLE() = default;
    STATE_TABLE(const STATE_TABLE&) = delete;  // disable copy constructor
};


// Haplotype-specific changes of each chain reconstructed at one node and site under model 3.
// Copy number is not computed here because it depends on the order of events,
// which cannot be uniquely determined from independent chain posteriors.
struct CN_CALL{
    int num_wgd = 0;
    pair<int, int> chr_change{0, 0};
    pair<int, int> site_change{0, 0};
    double prob = 0.0;          // posterior probability of the most likely chain state (marginal only)
};


/******************* functions for bounded model on site duplication/deletion *************/

void set_pmat(const evo_tree &rtree, int Ns, int nstate, int model, int cn_max, const vector<int> &knodes, vector<double> &blens, vector<double *> &pmat_per_blen, ofstream &fout);

void get_prob_line_orig(const vector<vector<double>> &L_sk_k, int nid, string &line);

void get_prob_line_orig(const vector<vector<double>> &L_sk_k, int nid, string &line);

void print_tree_state(const evo_tree& rtree, const vector<vector<int>>& S_sk_k, int nstate);

// Get the states of the tree from likelihood table at one site, starting from MRCA
// Step 5 of algorithm of Pupko (2000)
void extract_tree_ancestral_state(const evo_tree& rtree, const vector<int>& knodes, const vector<vector<int>>& S_sk_k, int model, map<int, int>& asr_states);

void set_tip_states(int model, int is_total, int obs_val, vector<int>& tip_states);

// Create likelihood vectors and state vectors at the tip node for reconstructing joint ancestral state
// L_sk_k (S_sk_k) has one row for each tree node and one column for each possible state
// Step 1 of algorithm of Pupko (2000)
void initialize_asr_table(const vector<int>& obs, const evo_tree& rtree, map<double, double*>& pmats, vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, int model, int nstate, int is_total);

// Find the most likely state for a node under each possible state, assuming current node has state nsk and parent node (connected by branch of length blen) has state np
// Step 2 of algorithm of Pupko (2000)
double get_max_prob_children(const vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree&rtree, double* pblen, int k, int nstate, int sp, int ni, int nj, double blen, int model, int debug);

void set_tip_states_tcn_decomp(int state_chr, int peak_sum_haplotype, int& si, int& ei, vector<int>& tip_states);

// Get the most likely state on one site of a chromosome (assuming higher level events on nodes)
// L_sk_k: the likelihood of the best reconstruction of the subtree at node sk provided that the father of sk has state k
// S_sk_k: the state of node sk provided that the father of sk has state k in the optimal conditional reconstruction
// Step 3 of algorithm of Pupko (2000)
void get_ancestral_states_site(vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, const vector<int>& knodes, const set<vector<int>>& comps, map<double, double*>& pmats, int nstate, int model);

// Get the probability of all possible copy numbers for one node at a site, for marginal reconstruction
string get_prob_line(const vector<vector<double>>& L_sk_k, int nid, int nchr, int nc, int is_total, int cn_max);

// Get the most likely copy number for each site at a specific node, for marginal reconstruction
void get_site_cnp(const vector<vector<double>>& L_sk_k, int nid, int nchr, int nc, int is_total, int cn_max, copy_number& cnp);

void print_node_cnp(ofstream& fout, const copy_number& cnp, int nid, int cn_max, int is_total);

// Infer the copy number of the MRCA given a tree at a site, assuming only site duplication/deletion
double reconstruct_marginal_ancestral_state(const evo_tree &rtree, const map<int, vector<vector<int>>> &vobs, const vector<int> &knodes, int model, int cn_max, int use_repeat, int is_total, string ofile);

void set_prob_line_header(int nstate, string& header);

// Infer the copy number of all internal nodes given a tree at a site, assuming only site duplication/deletion
void reconstruct_joint_ancestral_state(const evo_tree &rtree, const map<int, vector<vector<int>>> &vobs, vector<int> &knodes, int model, int cn_max, int use_repeat, int is_total, int m_max, string ofile);

void get_inode_cnp(int max_id, int Ns, map<int, int> &asr_states, int is_total, int cn_max, map<int, copy_number> &cnps, int nchr, int nc);

void write_joint_state_line(int max_id, int Ns, map<int, int>& asr_states, int nchr, int nc, int is_total, int cn_max, ofstream& fout);

void write_pattern_uniq(map<vector<int>, vector<vector<double>>>& sites_lnl_map);

void write_pattern_all(const map<int, vector<vector<int>>>& vobs);

/******************* functions for independent chain model on multi-level CNAs *************/

// Create likelihood vectors and state vectors at the tip node for reconstructing joint ancestral state, for independent chain model
// L_sk_k (S_sk_k) has one row for each tree node and one column for each possible state
void initialize_asr_table_decomp(const vector<CN_CHANGE>& obs, const evo_tree& rtree, const PMAT_DECOMP& pmat_decomp, const DIM_DECOMP& dim_decomp, LNL_TABLE& lnl_table, STATE_TABLE& state_table, const LNL_TYPE& lnl_type, bool use_edge_rates, int debug);

// Rate-matrix states compatible with an observed copy number change at a tip
void set_tip_states_decomp(const CN_CHANGE& obs_val, vector<int>& tip_states_chr, vector<int>& tip_states_site, const LNL_TYPE& lnl_type, int debug);

// Rate-matrix states compatible with an observed copy number change at a tip
void set_tip_states_decomp(const CN_CHANGE& obs_val, vector<int>& tip_states_chr, vector<int>& tip_states_site, const LNL_TYPE& lnl_type, int debug);

// Assume the likelihood table is for each combination of states
double get_max_children_decomp2(const vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, const set<vector<int>>& comps, int k, int nstate, double* pbli_wgd, double* pbli_chr, double* pbli_seg, DIM_DECOMP& dim_decomp, int sp, int ni, int nj, int blen);


// Get the ancestral state on one site of a chromosome
// Assuming each observed copy number is composed of three type of events.
// Sum over all possible states for initial and final nodes
// Allow at most one WGD event along a branch
void get_ancestral_states_site_decomp(LNL_TABLE& lnl_table, STATE_TABLE& state_table, const evo_tree& rtree, const vector<int>& knodes, const PMAT_DECOMP& pmat_decomp, const DIM_DECOMP& dim_decomp, bool use_edge_rates, int debug);


// Infer the copy number of the MRCA given a tree at a site, assuming independent Markov chains
// use_edge_rates: build one P matrix per edge from rtree.edge_rates instead of using the global rates
double reconstruct_marginal_ancestral_state_decomp(const evo_tree& rtree, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const vector<int>& knodes, const OBS_DECOMP& obs_decomp, const LNL_TYPE& lnl_type_orig, bool use_edge_rates, string ofile, int debug);


// Infer the copy number of all internal nodes given a tree at a site, assuming independent Markov chains
// use_edge_rates: as in reconstruct_marginal_ancestral_state_decomp
void reconstruct_joint_ancestral_state_decomp(const evo_tree &rtree, const map<int, vector<vector<CN_CHANGE>>> &vobs_change, vector<int> &knodes, const OBS_DECOMP &obs_decomp, const LNL_TYPE &lnl_type_orig, bool use_edge_rates, string ofile, int debug);

void initialize_lnl_state_table(LNL_TABLE& lnl_table, STATE_TABLE& state_table, int ntotn, DIM_DECOMP& dim_decomp);

// Copy number of a node from its WGD, chromosome-level and site-level states:
// total copy number if lnl_type.is_total, else the haplotype-specific state index
int get_cn_from_state_decomp(int state_wgd, int state_chr, int state_site, const LNL_TYPE& lnl_type);

// Posterior of each (cnA, cnB) at one site, from the posterior of each chain's
// states at the node. The chains are independent, so a combination of chain
// states has the product of their probabilities. An empty vector means the
// chain is not modelled and stays unchanged.
map<pair<int, int>, double> get_allele_cn_posterior_decomp(const vector<double>& post_wgd, const vector<double>& post_chr, const vector<double>& post_site, const LNL_TYPE& lnl_type);

// Most likely copy number at one site (total if lnl_type.is_total, else cnA and
// cnB), with its probability summed over all chain combinations giving it, and
// the most likely of those combinations
CN_CALL get_most_likely_cn_call_decomp(const vector<double>& post_wgd, const vector<double>& post_chr, const vector<double>& post_site, const LNL_TYPE& lnl_type);

// Header and lines of <ofile>.mrca.cn / <ofile>.joint.cn under model 3
string get_cn_call_header(int is_total, bool with_prob);
string get_cn_call_line(int nid, int nchr, int nc, const CN_CALL& call, int is_total, bool with_prob);

void set_change_wgd_chr(const vector<CN_CHANGE>& obs, const evo_tree& rtree, vector<int>& sample_num_wgd, vector<int>& change_chr);

void print_cn_mrca(copy_number_change& cn_mrca);

#endif
