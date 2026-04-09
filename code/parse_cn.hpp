#ifndef PARSE_CN_HPP
#define PARSE_CN_HPP


#include "common.hpp"   // for enum model
#include "stats.hpp"


// collection of functions for parsing input files


// https://github-pages.ucl.ac.uk/research-computing-with-cpp/04cpp3/sec03Templates.html 
// Following the first approach in the link to define template function in header file


// const int num_total_bins = 4401;

// template <typename T>
// vector<size_t> sort_indexes(const vector<T> &v){
//   // initialize original index locations
//   vector<size_t> idx(v.size());
//   iota(idx.begin(), idx.end(), 0);
//
//   // sort indexes based on comparing values in v
//   sort(idx.begin(), idx.end(),
//        [&v](size_t i1, size_t i2){return v[i1] < v[i2];});
//
//   return idx;
// }


/**************** function used in reconstructing ancestral states ****************/
// Convert the haplotype-specific state (ordered by 0/0	0/1	1/0	0/2	 1/1	2/0	0/3	 1/2 2/1	3/0	 0/4 1/3 2/2 3/1	4/0) to total copy number
int state_to_total_cn(int state, int cn_max);
// Convert the haplotype-specific state (ordered by 0/0	0/1	1/0	0/2	 1/1	2/0	0/3	 1/2 2/1	3/0	 0/4 1/3 2/2 3/1	4/0) to haplotype-specific copy number
void state_to_allele_cn(int state, int cn_max, int& cnA, int& cnB);

void set_nstates(const int max_change_haplotype, int& cn_max, vector<int>& n_states, vector<int>& sums);

int change_state_to_total_cn(int state, int max_change_haplotype);

/*****************************************/
// Change the allele specific copy number to the state used in substitution rate matrix
int allele_cn_to_state(int cnA, int cnB);

void rcn_to_decomposition(const vector<vector<vector<int>>>& s_info, vector<vector<CN_CHANGE>>& s_info_change, vector<vector<int>>& sample_change_site, vector<int>& site_max_change, int debug);

// Convert the copy number matrix to decomposition format
void cn_to_decomposition(const vector<vector<vector<int>>>& s_info, 
                        vector<vector<CN_CHANGE>>& s_info_change, 
                        const vector<vector<int>>&  sample_change_chr, 
                        const vector<vector<int>>& sample_change_site, 
                        const vector<int>& sample_num_wgd, 
                        const INPUT_PROPERTY& input_property, 
                        int debug);


bool adjust_chr_change_for_wgd(int& rounded_num_change, double avg_cn, int nwgd);

bool adjust_site_change_for_wgd(int& rounded_num_change, double avg_cn, int nwgd);
 
// Find the potential number of WGDs for each sample
void get_num_wgd(const vector<double>& sample_avg_cn, vector<int>& sample_num_wgd, int debug = 0);

// Find the potential number of chromosome changes for each sample
void get_chr_change(const vector<int>& sample_num_wgd, 
                    const vector<double>& sample_avg_cn, 
                    const vector<map<int, vector<int>>>& sample_chr_cn, 
                    vector<vector<int>>& sample_change_chr, 
                    vector<int>& chr_max_change, 
                    int max_chr_change_haplotype, 
                    int debug = 0); 

void get_chr_change_haplotype(const vector<int>& sample_num_wgd, 
                                const vector<double>& sample_avg_cnA, 
                                const vector<double>& sample_avg_cnB, 
                                const vector<map<int, vector<int>>>& sample_chr_cnA, 
                                const vector<map<int, vector<int>>>& sample_chr_cnB, 
                                vector<vector<int>>& sample_change_chr, 
                                vector<int>& chr_max_change, 
                                int max_chr_change_haplotype, 
                                int debug = 0);

// Find the potential number of site changes for each sample
void get_site_change(const vector<int>& sample_num_wgd, 
                    const vector<vector<vector<int>>>& s_info, 
                    const vector<double>& sample_avg_cn, 
                    vector<vector<int>>& sample_change_site, 
                    vector<int>& site_max_change, 
                    int max_site_change_haplotype,  
                    int debug = 0);

void get_site_change_haplotype(const vector<int>& sample_num_wgd, 
                                const vector<vector<vector<int>>>& s_info, 
                                vector<double>& sample_avg_cnA, 
                                vector<double>& sample_avg_cnB, 
                                vector<vector<int>>& sample_change_site, 
                                vector<int>& site_max_change, 
                                int cn_max, 
                                int max_site_change_haplotype, 
                                int debug = 0);

// Compute the average copy number of a sample
double compute_sample_avg_cn(const vector<vector<int>>& s_cn, map<int, vector<int>>& chr_cn, int cn_max, int is_total);

// Compute the average copy number of all samples
void get_sample_ploidy(const vector<vector<vector<int>>>& s_info, vector<double>& sample_avg_cn, vector<map<int, vector<int>>>& sample_chr_cn, int cn_max, int is_total, int debug);

// Compute the average copy number of a sample for haplotype-specific copy number
pair<double, double> compute_sample_avg_cn_haplotype(const vector<vector<int>>& s_cn, map<int, vector<int>>& chr_cnA, map<int, vector<int>>& chr_cnB, int cn_max);  

// Compute the average copy number of all samples for haplotype-specific copy number
void get_sample_ploidy_haplotype(const vector<vector<vector<int>>>& s_info, vector<double>& sample_avg_cnA, vector<map<int, vector<int>>>& sample_chr_cnA, vector<double>& sample_avg_cnB, vector<map<int, vector<int>>>& sample_chr_cnB, int cn_max, int is_total, int debug);


// Find the largest copy number in a sample
void get_sample_mcn(const vector<vector<vector<int>>>& s_info, vector<int>& sample_max_cn, int cn_max, int is_total);

// Find the number of invariable sites for each character (state)
// Loop over and output only the regions that have varied, which provide information for tree building
void get_var_bins(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, vector<int>& var_bins, int is_total, int debug);

// Distinguish invariable and variable sites;
// Combine adjacent invariable sites
// TODO: Distinguish different types of invariant sites: 2 vs 4, which have different likelihoods
vector<vector<int>> get_invar_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int is_total = 1, int debug = 0);

vector<vector<int>> get_all_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int incl_all, int is_total = 1, int debug = 0);

// read segments and copy numbers into different data structures for further processing
int get_segs_cn(vector<vector<vector<int>>>& s_info, vector<vector<int>>& segs, const string& filename, const INPUT_PROPERTY& input_property, INPUT_DATA& input_data, int debug);


// Compute the average copy number of a segment

void compute_segment_cn_state(vector<int>& seg_cn, const vector<int>& segs_curr, const vector<vector<vector<int>>>& s_info, int Ns, int cn_max, int debug);

void compute_segment_cn_change(vector<CN_CHANGE>& seg_cn, const vector<int>& segs_curr, const vector<vector<CN_CHANGE>>& s_info_change, int Ns, int Nseg, int debug);


void group_segs_by_chr_state(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, map<int, vector<vector<int>>>& data, int Ns, int cn_max, const string& seg_file, int debug);

void group_segs_by_chr_change(const vector<vector<int>>& segs, const vector<vector<CN_CHANGE>>& s_info_change, map<int, vector<vector<int>>>& data_change, int Ns, int cn_max, const string& seg_file, int debug);


vector<vector<int>> group_segs(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, int Ns, int cn_max, int debug = 0);

// Get the input matrix of copy numbers by chromosome

void get_obs_vector_by_chr_state(const map<int, vector<vector<int>>>& data,  map<int, vector<vector<int>>>& vobs, int Ns);

void get_obs_vector_by_chr_change(const map<int, vector<vector<int>>>& data_change,  map<int, vector<vector<CN_CHANGE>>>& vobs_change, int Ns);


/** 
 * @brief Generate bootstrap resampled copy number observations by chromosome for tree inference
 * @param data Input data map where each key is a chromosome number and the value is a vector of vectors containing copy number information.
 * @param vobs Output parameter to store bootstrap resampled copy number observations.
 * @param r GSL random number generator.
 * @tparam T Type of the copy number data (e.g., int, CN_CHANGE), i.e. vector<vector<CN_CHANGE>> or vector<vector<int> 
 */
template <typename T>
void get_bootstrap_vector_by_chr(map<int, vector<vector<int>>>& data, map<int, vector<vector<T>>>& vobs, gsl_rng* r){
    // create a copy of vobs to resample
    map<int, vector<vector<T>>> vobs_copy = vobs;
    vobs.clear();

    // cout << "Total number of chromosomes " << total_chr << endl;
    for(auto dcn : data){
        int nchr = dcn.first;
        // cout << "Chr " << nchr << "\t";
        vector<vector<T>> obs_chr;
        for(size_t nc = 0; nc < data[nchr].size(); ++nc){
            // randomly select a site
            size_t i = gsl_rng_uniform_int(r, data[nchr].size());
            // cout << i << ":" << vobs_copy[nchr][i].size() << "\t" ;
            obs_chr.push_back(vobs_copy[nchr][i]);
        }
        // cout << obs_chr.size() << endl;
        vobs[nchr] = obs_chr;
    }
}

/******************* read input file ***********************/
// Read the sampling time and patient age of each sample
// Ns: number of samples
// age: age at earliest sample
// Assume samples are ordered by ID from 1 to Ns
vector<double> read_time_info(const string& filename, const int Ns, int& age, int debug = 0);


// parse meta information file, including known number of WGDs, chromosome gain/loss, and site duplication/deletion inferred from observed copy numbers
// meta file format: sample_id num_wgd chr_change site_change
// sample_id: 1..Ns
// num_wgd: integer, number of WGD events inferred from observed copy numbers, required when meta file is provided
// chr_change: optional, comma-separated list of integer copy number changes for each chromosome, e.g. "0,1,-1" for no change, one gain, one loss. If haplotype-specific changes are provided, the changes should be ordered by haplotype, e.g. "0/0,1/-1" for no change, one gain on haplotype A and one loss on haplotype B. The chromosome changes should have been adjusted according to WGD events to avoid double counting when WGD events are provided.
// site_change: optional, semicolon-separated list of chromosome-specific site changes, optional, e.g. "1:0,1,-1;2:0,0,2;22:-1,1" for changes of three sites on chr 1, three sites on chr 2, and two sites on chr 22 and "1:0/1,1/-1;2:0/0,2/1" for haplotype-specific changes. The chromosome id should be between 1 and 22, and the site changes should be comma-separated integers for each chromosome. The site changes should have been adjusted according to chromosome gain/loss to avoid double counting when chromosome changes are provided.  
// TODO: change it to one line 
void read_meta_info(const string& meta_file,
                    int Ns,
                    vector<int>& sample_num_wgd,
                    vector<vector<int>>& sample_change_chr,
                    vector<vector<int>>& sample_change_site,
                    int is_total,
                    bool debug);


void compute_max_changes(const vector<vector<int>>& sample_change_chr,
                        const vector<vector<int>>& sample_change_site,
                        vector<int>& chr_max_change,
                        vector<int>& site_max_change,
                        const vector<pair<int,int>>& states_chr,
                        const vector<pair<int,int>>& states_site,
                        int is_total,
                        int debug);


// Format of input total copy number: sample, chr, seg, copy_number
// Format of input allele specific copy number: sample, chr, seg, copy_number A, copy_number B -> converted into a specific state by decomposition
// s_info: a vector for each sample, which is composed of a vector of vector<int> vcn{chr, sid, cn}
void read_cn(vector<vector<vector<int>>>& s_info, const string& filename, int Ns, int &num_total_bins, int cn_max, int is_total = 1, int is_rcn = 0, int debug = 0);


/******************* read input all in one ***********************/

// Read the input copy numbers
vector<vector<int>> read_data_var_regions(const string& filename, const INPUT_PROPERTY& input_property, INPUT_DATA& input_data, int debug = 0);


void read_data_var_regions_by_chr_state(map<int, vector<vector<int>>>& data, const string& cn_file, const string& seg_file, const INPUT_PROPERTY& input_property, INPUT_DATA& input_data, int debug);

void read_data_var_regions_by_chr_change(map<int, vector<vector<int>>>& data_change, const string& cn_file, const string& seg_file, const string& meta_file, const INPUT_PROPERTY& input_property, INPUT_DATA& input_data, int debug);


#endif
