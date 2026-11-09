#ifndef PARSE_CN_HPP
#define PARSE_CN_HPP

#include "common.hpp"   // for enum model
#include "stats.hpp"

// collection of functions for parsing input files

// https://github-pages.ucl.ac.uk/research-computing-with-cpp/04cpp3/sec03Templates.html
// Following the first approach in the link to define template function in header file

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
/**
 * @brief Convert a state (0-based index of copy number combinations used for rate matrix) to total copy number.
 * @param state The state index representing a specific copy number combination.
 * @param cn_max The maximum copy number allowed by the program.
 * @return The total copy number corresponding to the given state.
 * 0	1	2	3	4	5	6	7	8	9	10	11	12	13	14
 * 0	1	1	2	2	2	3	3	3	3	4	4	4	4	4
 * 0/0	0/1	1/0	0/2	 1/1	2/0	0/3	 1/2	 2/1	3/0	0/4	 1/3	  2/2	 3/1	4/0
 */
int state_to_total_cn(int state, int cn_max);
/**
 * @brief Convert a haplotype-specific state to its copy numbers (cnA, cnB); the inverse of allele_cn_to_state().
 * States are ordered by total copy number, then by cnA: 0/0 0/1 1/0 0/2 1/1 2/0 0/3 1/2 2/1 3/0 ...
 * Exits if the state is outside 0 .. (cn_max + 1)(cn_max + 2)/2 - 1.
 * @param state The state index.
 * @param cn_max The maximum total copy number, which sets the number of states.
 * @param cnA Output: copy number of haplotype A.
 * @param cnB Output: copy number of haplotype B.
 */
void state_to_allele_cn(int state, int cn_max, int& cnA, int& cnB);

/**
 * @brief Count the haplotype-specific change states for each total copy number change.
 * @param max_change_haplotype Largest change allowed on one haplotype; must be at least 1.
 * @param cn_max Output: the largest total copy number, 2 * (max_change_haplotype + 1).
 * @param n_states Output: number of states for each total copy number, from 0 to cn_max.
 * @param sums Output: cumulative sums of n_states.
 */
void set_nstates(int max_change_haplotype, int& cn_max, vector<int>& n_states, vector<int>& sums);

/**
 * @brief Convert a change state (0-based index of copy number combinations used for rate matrix) to total copy number.
 * @param state The state index representing a specific copy number combination.
 * @param max_change_haplotype Largest change allowed on one haplotype, which sets the number of states.
 * @return The total copy number corresponding to the given state.
 * When max_change_haplotype = 1, 9 states
 * When max_change_haplotype = 2, 16 states
 */
int change_state_to_total_cn(int state, int max_change_haplotype);


/*****************************************/
/**
 * @brief Convert haplotype-specific copy numbers to the state index used in the substitution rate matrix.
 * States are ordered by total copy number, then by cnA: 0/0 0/1 1/0 0/2 1/1 2/0 0/3 ...
 * @param cnA Copy number of haplotype A.
 * @param cnB Copy number of haplotype B.
 * @return The state index.
 */
int allele_cn_to_state(int cnA, int cnB);

/**
 * @brief Find the range of haplotype-specific change states that share an observed total copy number change.
 * Values outside the valid range are clamped to it.
 * @param cn_total Total copy number change plus 2, so that it is non-negative.
 * @param peak_sum_haplotype max_change_haplotype + 2, which sets the number of states.
 * @param si Output: index of the first matching state. Must be 0 on entry, as it is accumulated.
 * @param ei Output: index of the last matching state.
 */
void get_tcn_state_index(int cn_total, int peak_sum_haplotype, int& si, int& ei);

/**
 * @brief Decompose relative total copy numbers (--is_rcn) into multi-level changes. No WGD or chromosome-level change is assumed: the whole difference from NORM_PLOIDY is a site-level change.
 * @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 * @param s_info_change Output: decomposed copy number changes for each sample, same dimensions as s_info.
 * @param sample_change_site Output: site-level change of every site, appended per sample.
 * @param site_max_change Output: largest absolute site-level change in each sample.
 * @param debug Debug flag for verbose output.
 */
void rcn_to_decomposition(const vector<vector<vector<int>>>& s_info, vector<vector<CN_CHANGE>>& s_info_change, vector<vector<int>>& sample_change_site, vector<int>& site_max_change, int debug);

/**
 * @brief Decompose total copy number into multi-level changes (WGD, chromosome, segment) given original data and computed changes.
 * @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 * @param s_info_change A vector of CN_CHANGE structs to store the decomposed copy number changes for each sample, same dimensions as s_info.
 * @param sample_change_chr A vector of vectors containing chromosome-level change states for each sample, indexed by chromosome id (0-based).
 * @param sample_change_site A vector of vectors containing site-level change states for each sample, indexed by segment id.
 * @param sample_num_wgd A vector containing the number of WGD events for each sample.
 * @param input_property struct containing sample-level information needed for decomposition.
 * @param debug Debug flag for verbose output.
 */
void cn_to_decomposition(const vector<vector<vector<int>>>& s_info,
                        vector<vector<CN_CHANGE>>& s_info_change,
                        const vector<vector<int>>&  sample_change_chr,
                        const vector<vector<int>>& sample_change_site,
                        const vector<int>& sample_num_wgd,
                        const INPUT_PROPERTY& input_property,
                        int debug);

/**
 * @brief Reinterpret a chromosome-level change as possibly before a WGD.
 * If the sample probably has a WGD and the change is even, it is halved and multiplied by CHANGE_CHR
 * to mark it as ambiguous (before or after the WGD); otherwise it is left unchanged.
 * @param rounded_num_change The rounded chromosome-level change; updated in place.
 * @param avg_cn Average copy number (ploidy) of the sample.
 * @param nwgd Number of WGDs in the sample.
 * @return true if the change was adjusted.
 */
bool adjust_chr_change_for_wgd(int& rounded_num_change, double avg_cn, int nwgd);

/**
 * @brief Reinterpret a site-level change as before a WGD.
 * If the sample probably has a WGD and the change is even, it is halved; otherwise it is left unchanged.
 * @param rounded_num_change The rounded site-level change; updated in place.
 * @param avg_cn Average copy number (ploidy) of the sample.
 * @param nwgd Number of WGDs in the sample.
 * @return true if the change was halved.
 */
bool adjust_site_change_for_wgd(int& rounded_num_change, double avg_cn, int nwgd);

/**
 *  @brief Estimate the potential number of WGDs for each sample
 *  Ploidy above WGD_CUTOFF gives one WGD, above twice WGD_CUTOFF gives two.
 *  @param sample_avg_cn Average copy number (ploidy) of each sample.
 *  @param sample_num_wgd A vector to store the estimated number of WGDs for each sample in the input.
 *  @param debug Debug flag for verbose output.
 */
void get_num_wgd(const vector<double>& sample_avg_cn, vector<int>& sample_num_wgd, int debug = 0);

/**
 *  @brief Estimate the potential number of total chromosome changes for each sample based on known number of WGD or estimated average copy number for each sample
 *  @param sample_num_wgd Number of WGDs in each sample.
 *  @param sample_avg_cn A vector to store the average copy number for each sample in the input, required for estimating chr-level copy number changes.
 *  @param sample_chr_cn A vector of maps to store copy numbers for each chromosome in each sample.
 *  @param sample_change_chr A vector to store the estimated number of chromosome changes for each sample in the input, indexed by [sample][chromosome].
 *  @param chr_max_change A vector to store the maximum absolute chromosome change for each sample in the input.
 *  @param max_chr_change_haplotype Maximum copy number change for haplotype-specific chr-level events.
 *  @param debug Debug flag for verbose output.
 */
void get_chr_change(const vector<int>& sample_num_wgd,
                    const vector<double>& sample_avg_cn,
                    const vector<map<int, vector<int>>>& sample_chr_cn,
                    vector<vector<int>>& sample_change_chr,
                    vector<int>& chr_max_change,
                    int max_chr_change_haplotype,
                    int debug = 0);

/**
 *  @brief Estimate the potential number of haplotype-specific chromosome changes for each sample
 *  @param sample_num_wgd Number of WGDs in each sample.
 *  @param sample_avg_cnA Average copy number of haplotype A in each sample.
 *  @param sample_avg_cnB Average copy number of haplotype B in each sample.
 *  @param sample_chr_cnA Copy numbers of haplotype A grouped by chromosome, for each sample.
 *  @param sample_chr_cnB Copy numbers of haplotype B grouped by chromosome, for each sample.
 *  @param sample_change_chr A vector to store the estimated number of chromosome changes for each sample in the input, indexed by [sample][chromosome].
 *  @param chr_max_change A vector to store the maximum absolute chromosome change for each sample in the input.
 *  @param max_chr_change_haplotype Maximum copy number change for haplotype-specific chr-level events.
 *  @param debug Debug flag for verbose output.
 */
void get_chr_change_haplotype(const vector<int>& sample_num_wgd,
                                const vector<double>& sample_avg_cnA,
                                const vector<double>& sample_avg_cnB,
                                const vector<map<int, vector<int>>>& sample_chr_cnA,
                                const vector<map<int, vector<int>>>& sample_chr_cnB,
                                vector<vector<int>>& sample_change_chr,
                                vector<int>& chr_max_change,
                                int max_chr_change_haplotype,
                                int debug = 0);

/**
 *  @brief Compute per-sample site-level max deviation from avg_cn (baseline) using total copy numbers, then take ceiling as an integer bound.
 *  @param sample_num_wgd Number of WGDs in each sample.
 *  @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 *  @param sample_avg_cn Average copy number (ploidy) of each sample.
 *  @param sample_change_site Output: site-level change of every site, appended per sample.
 *  @param site_max_change A vector to store the estimated number of site changes for each sample in the input.
 *  @param max_site_change_haplotype Maximum site-level copy number change for haplotype-specific model, not used actually for now
 *  @param debug Debug flag for verbose output.
 */
void get_site_change(const vector<int>& sample_num_wgd,
                    const vector<vector<vector<int>>>& s_info,
                    const vector<double>& sample_avg_cn,
                    vector<vector<int>>& sample_change_site,
                    vector<int>& site_max_change,
                    int max_site_change_haplotype,
                    int debug = 0);

/**
 *  @brief Compute per-sample site-level max deviation from avg_cn (baseline) using haplotype-specific copy numbers, then take ceiling as an integer bound.
 *  @param sample_num_wgd Number of WGDs in each sample.
 *  @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 *  @param sample_avg_cnA Average copy number of haplotype A in each sample.
 *  @param sample_avg_cnB Average copy number of haplotype B in each sample.
 *  @param sample_change_site Output: site-level change state of every site, appended per sample.
 *  @param site_max_change A vector to store the estimated number of site changes for each sample in the input.
 *  @param cn_max Maximum copy number, used for converting state back to copy number for haplotype-specific input.
 *  @param max_site_change_haplotype Maximum copy number change for haplotype-specific site-level events, used for determining the state space of site-level changes and the dimension of rate matrix.
 *  @param debug Debug flag for verbose output.
 */
void get_site_change_haplotype(const vector<int>& sample_num_wgd,
                                const vector<vector<vector<int>>>& s_info,
                                vector<double>& sample_avg_cnA,
                                vector<double>& sample_avg_cnB,
                                vector<vector<int>>& sample_change_site,
                                vector<int>& site_max_change,
                                int cn_max,
                                int max_site_change_haplotype,
                                int debug = 0);

/**
 *  @brief Calculate average total copy number or ploidy for a single sample.
 *  @param s_cn copy number information for a sample in the format [chr, sid, cn].
 *  @param chr_cn A map to store copy numbers for each chromosome.
 *  @param cn_max Maximum copy number.
 *  @param is_total Indicator whether the copy numbers are total (1) or haplotype-specific (0).
 *  @return The average copy number or ploidy for the sample.
 */
double compute_sample_avg_cn(const vector<vector<int>>& s_cn, map<int, vector<int>>& chr_cn, int cn_max, int is_total);

/**
 * @brief compute ploidy for either total copy number or haplotype-specific copy number, depending on the input format
 * @param s_info: real CN for total copy number, state for haplotype-specific copy number
 * @param sample_avg_cn: average totalcopy number for each sample, used for WGD estimation and chromosome change estimation
 * @param sample_chr_cn: raw copy numbers for each chromosome in each sample, used for chromosome change estimation, total CN for total copy number, state for haplotype-specific copy number
 * @param cn_max: maximum copy number
 * @param is_total: indicator whether the copy numbers are total (1) or haplotype-specific (0)
 * @param debug: debug flag for verbose output
 */
void get_sample_ploidy(const vector<vector<vector<int>>>& s_info, vector<double>& sample_avg_cn, vector<map<int, vector<int>>>& sample_chr_cn, int cn_max, int is_total, int debug);

/**
 * @brief Calculate the average copy number of each haplotype for a single sample.
 * @param s_cn Copy number information for a sample in the format [chr, sid, state].
 * @param chr_cnA Output: copy numbers of haplotype A, grouped by chromosome.
 * @param chr_cnB Output: copy numbers of haplotype B, grouped by chromosome.
 * @param cn_max Maximum copy number, used to decode the states.
 * @return The average copy numbers of haplotypes A and B.
 */
pair<double, double> compute_sample_avg_cn_haplotype(const vector<vector<int>>& s_cn, map<int, vector<int>>& chr_cnA, map<int, vector<int>>& chr_cnB, int cn_max);

/**
 * @brief Compute the average copy number of each haplotype for every sample (haplotype-specific input only).
 * @param s_info Haplotype-specific states for each sample, in the format [chr, sid, state].
 * @param sample_avg_cnA Output: average copy number of haplotype A in each sample.
 * @param sample_chr_cnA Output: copy numbers of haplotype A grouped by chromosome, for each sample.
 * @param sample_avg_cnB Output: average copy number of haplotype B in each sample.
 * @param sample_chr_cnB Output: copy numbers of haplotype B grouped by chromosome, for each sample.
 * @param cn_max Maximum copy number, used to decode the states.
 */
void get_sample_ploidy_haplotype(const vector<vector<vector<int>>>& s_info, vector<double>& sample_avg_cnA, vector<map<int, vector<int>>>& sample_chr_cnA, vector<double>& sample_avg_cnB, vector<map<int, vector<int>>>& sample_chr_cnB, int cn_max);

/**
 * @brief Find the largest total copy number in each sample, for diagnostics.
 * @param s_info Copy numbers for each sample, in the format [chr, sid, cn]; states for haplotype-specific input.
 * @param sample_max_cn Output: largest total copy number of each sample.
 * @param cn_max Maximum copy number, used to decode haplotype-specific states.
 * @param is_total 1 for total copy number, 0 for haplotype-specific states.
 */
void get_sample_mcn(const vector<vector<vector<int>>>& s_info, vector<int>& sample_max_cn, int cn_max, int is_total);

/**
 * @brief Identify variable bins across all samples for segment merging.
 * @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 * @param Ns The number of samples.
 * @param num_total_bins The total number of bins.
 * @param num_invar_bins A reference to an integer that will store the number of invariable bins found.
 * @param var_bins A reference to a vector that will store flags indicating whether each bin is variable (1) or invariable (0).
 * @param is_total An integer flag indicating whether the copy numbers are total (1) or haplotype-specific (0).
 * @param debug Debug flag for verbose output.
 */
void get_var_bins(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, vector<int>& var_bins, int is_total, int debug);

/**
 * @brief Merge consecutive variable bins with the same copy number profiles across all samples into segments.
 * TODO: Distinguish different types of invariant sites: 2 vs 4, which have different likelihoods
 * @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample.
 * @param Ns The number of samples.
 * @param num_total_bins The total number of bins.
 * @param num_invar_bins A reference to an integer that will store the number of invariable bins found.
 * @param is_total An integer flag indicating whether the copy numbers are total (1) or haplotype-specific (0).
 * @param debug Debug flag for verbose output.
 * @return A vector of vectors representing the merged segments, where each seg is represented by: {chr, id_start, id_end, seg_start, seg_end}
 */
vector<vector<int>> get_invar_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int is_total = 1, int debug = 0);

/**
 * @brief Make a one-bin segment for every bin, or only for variable bins, without merging.
 * @param s_info Copy numbers for each sample, in the format [chr, sid, cn].
 * @param Ns Number of samples.
 * @param num_total_bins Total number of bins.
 * @param num_invar_bins Output: number of invariable bins.
 * @param incl_all 1 to keep every bin, 0 to keep only variable bins (invariable ones are then accounted for by their count).
 * @param is_total 1 for total copy number, 0 for haplotype-specific states.
 * @param debug Debug flag for verbose output.
 * @return The segments, each represented by: {chr, id_start, id_end, seg_start, seg_end}
 */
vector<vector<int>> get_all_segs(const vector<vector<vector<int>>>& s_info, int Ns, int num_total_bins, int& num_invar_bins, int incl_all, int is_total = 1, int debug = 0);

/**
 * @brief Read the input copy number file and merge bins into segments if specified
 * @param s_info Output parameter to store copy number information for each sample.
 * @param segs Output parameter to store segment information.
 * @param filename Input file name containing copy number data.
 * @param input_property Structure containing properties of the input data.
 * @param input_data Structure to store processed input data.
 * @param debug Debug flag for verbose output.
 * @return The maximum copy number state used in the analysis of bounded models. Different from cn_max in the input_property when haplotype-specific CNs are used.
 */
int get_segs_cn(vector<vector<vector<int>>>& s_info, vector<vector<int>>& segs, const string& filename, const INPUT_PROPERTY& input_property, INPUT_DATA& input_data, int debug);

/**
 * @brief Get the copy number (state) of a segment in each sample; exits if its bins disagree or exceed cn_max.
 * @param seg_cn Output: copy number of the segment in each sample, appended in sample order
 * @param segs_curr The segment, represented by: {chr, id_start, id_end, seg_start, seg_end}
 * @param s_info Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample
 * @param Ns Number of samples.
 * @param cn_max Maximum copy number state, different from maximum copy number value in the input when the input is haplotype-specific copy number.
 * @param debug Debug flag for verbose output
 */
void compute_segment_cn_state(vector<int>& seg_cn, const vector<int>& segs_curr, const vector<vector<vector<int>>>& s_info, int Ns, int cn_max, int debug);

/**
 * @brief Get the copy number changes of a segment in each sample, taken from its first bin
 * @param seg_cn Output: copy number changes of the segment in each sample, appended in sample order
 * @param segs_curr Vector of segments represented by: {chr, id_start, id_end, seg_start, seg_end}, indicating segment locations
 * @param s_info_change Sample information where each innermost vector contains copy number for a segment in the format [chr, sid, cn] for each sample
 * @param Ns Number of samples.
 * @param Nsite Number of copy number sites.
 * @param debug Debug flag for verbose output
 */
void compute_segment_cn_change(vector<CN_CHANGE>& seg_cn, const vector<int>& segs_curr, const vector<vector<CN_CHANGE>>& s_info_change, int Ns, int Nsite, int debug);

/**
 * @brief Merge consecutive bins with the same copy numbers across all samples
 * @param segs Input segments represented by: {chr, id_start, id_end, seg_start, seg_end}
 * @param s_info Sample information where each innermost vector contains copy number state (int) or change (CN_CHANGE) for a segment in the format [chr, sid, cn] for each sample.
 * @param data Output: a map where each key is a chromosome number and the value is a vector of segments with their copy numbers for each sample.
 * @param Ns Number of samples.
 * @param cn_max Maximum copy number.
 * @param seg_file Output file name to write the segments with copy numbers for each sample.
 * @param debug Debug flag for verbose output.
 */
void group_segs_by_chr_state(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, map<int, vector<vector<int>>>& data, int Ns, int cn_max, const string& seg_file, int debug);

/**
 * @brief Merge consecutive bins with the same copy numbers across all samples
 * @param segs Input segments represented by: {chr, id_start, id_end, seg_start, seg_end}
 * @param s_info_change Sample information where each innermost vector contains copy number state (int) or change (CN_CHANGE) for a segment in the format [chr, sid, cn] for each sample.
 * @param data_change A map where each key is a chromosome number and the value is a vector of segments with their copy numbers for each sample.
 * @param Ns Number of samples.
 * @param Nsite Number of copy number sites.
 * @param seg_file Output file name to write the segments with copy numbers for each sample.
 * @param debug Debug flag for verbose output.
 */
void group_segs_by_chr_change(const vector<vector<int>>& segs, const vector<vector<CN_CHANGE>>& s_info_change, map<int, vector<vector<int>>>& data_change, int Ns, int Nsite, const string& seg_file, int debug);

/**
 * @brief Build the copy number matrix: one row per segment with its copy number in each sample.
 * @param segs Segments, each represented by: {chr, id_start, id_end, seg_start, seg_end}
 * @param s_info Copy numbers for each sample, in the format [chr, sid, cn].
 * @param Ns Number of samples.
 * @param cn_max Maximum total copy number, or maximum state for haplotype-specific copy number; used to validate input CNs.
 * @param debug Debug flag for verbose output.
 * @return One row per segment: {chr, id_start, id_end, cn_sample1, ..., cn_sampleNs}
 */
vector<vector<int>> group_segs(const vector<vector<int>>& segs, const vector<vector<vector<int>>>& s_info, int Ns, int cn_max, int debug = 0);


/*****************  Get the input matrix of copy numbers by chromosome ********/
/**
 *  @brief Group observed copy numbers by chromosome to construct the CN matrix for tree inference
 *  @param data Input data map where each key is a chromosome number and the value is a vector of vectors containing copy number information.
 *  @param vobs A map where each key is a chromosome number and the value is a vector of vectors containing copy number observations for each sample.
 *  vobs[chr] = { {cn_sample1, cn_sample2, ...}, ...} }
 *  @param Ns Number of samples.
 */
void get_obs_vector_by_chr_state(const map<int, vector<vector<int>>>& data,  map<int, vector<vector<int>>>& vobs, int Ns);

/**
 *  @brief Group observed copy number changes by chromosome for tree inference
 *  @param data_change Input data map where each key is a chromosome number and the value is a vector of vectors containing copy number information.
 *  @param vobs_change A map where each key is a chromosome number and the value is a vector of vectors containing copy number observations for each sample.
 *  vobs[chr] = { {cn_change_sample1, cn_change_sample2, ...}, ...} }
 *  @param Ns Number of samples.
 */
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
/**
 * @brief Read the sampling time of each sample and the patient age at the earliest sample.
 * Assumes samples are ordered by ID from 1 to Ns.
 * @param filename The name of the input file containing time data.
 * @param Ns The number of samples to read from the file.
 * @param age A reference to an integer that will store the minimum age read from the file.
 * @param debug An integer flag for enabling debug output.
 */
vector<double> read_time_info(const string& filename, const int Ns, int& age, int debug = 0);

/**
 * @brief Parse the meta information file: known number of WGDs, and optionally chromosome gain/loss and site duplication/deletion inferred from observed copy numbers, for total or haplotype-specific copy number changes.
 * TODO: move the format description to the File format section of the documentation
 * Columns (tab-separated): sample_id num_wgd chr_change site_change
 * sample_id: 1..Ns (required)
 * num_wgd: integer, number of WGD events inferred from observed copy numbers, required when meta file is provided
 * chr_change: optional, comma-separated list of integer copy number changes for each chromosome, e.g. "0,1,-1" for no change, one gain, one loss. If haplotype-specific changes are provided, the changes should be ordered by haplotype, e.g. "0/0,1/-1" for no change, one gain on haplotype A and one loss on haplotype B. The chromosome changes should have been adjusted according to WGD events to avoid double counting when WGD events are provided.
 * site_change: optional, semicolon-separated list of chromosome-specific site changes, e.g. "1:0,1,-1;2:0,0,2;22:-1,1" for changes of three sites on chr 1, three sites on chr 2, and two sites on chr 22 and "1:0/1,1/-1;2:0/0,2/1" for haplotype-specific changes. The chromosome id should be between 1 and 22, and the site changes should be comma-separated integers for each chromosome. The site changes should have been adjusted according to chromosome gain/loss to avoid double counting when chromosome changes are provided.
 * @param meta_file Path to the meta file; an empty path leaves every sample with 0 WGDs and no changes.
 * @param Ns Number of samples.
 * @param sample_num_wgd Output: number of WGDs in each sample.
 * @param sample_change_chr Output: chromosome-level changes of each sample (empty when not given).
 * @param sample_change_site Output: site-level changes of each sample (empty when not given).
 * @param input_property Input properties; is_total and the max_*_change_haplotype limits decide how changes are parsed.
 * @param debug Debug flag for verbose output.
 */
void read_meta_info(const string& meta_file,
                    int Ns,
                    vector<int>& sample_num_wgd,
                    vector<vector<int>>& sample_change_chr,
                    vector<vector<int>>& sample_change_site,
                    const INPUT_PROPERTY& input_property,
                    bool debug);

// Helpers for parsing the chr_change and site_change fields of the meta file

/**
 * @brief Parse one change token from the meta file.
 * Total CN: an integer, e.g. "-1". Haplotype-specific CN: "a/b", e.g. "1/-1", encoded as its build_pair_states() index.
 * Throws if the token is empty, malformed, or out of range.
 * @param token The token to parse.
 * @param is_total 1 for total copy number changes, 0 for haplotype-specific changes.
 * @param max_haplotype_change Largest change allowed on one haplotype (the smallest is -1).
 * @param context Description of where the token comes from, used only in error messages.
 * @return The change value or state index.
 */
int parse_change_token(const string& token,
                       int is_total,
                       int max_haplotype_change,
                       const string& context);

/**
 * @brief Parse a comma-separated list of copy number changes from the meta file into integer change values.
 * Used by read_meta_info() for the chromosome-change column (one value per chromosome, in chromosome order)
 * and by parse_site_change_field() for the site changes of each chromosome block.
 * Total CN: "0,1,-1" -> {0, 1, -1}.
 * Haplotype-specific CN: "0/0,1/-1" -> state indices in build_pair_states() order, e.g. "0/0" -> 4.
 * @param s The field to parse; a missing field gives an empty vector.
 * @param is_total 1 for total copy number changes, 0 for haplotype-specific changes written as a/b.
 * @param max_haplotype_change Largest change allowed on one haplotype (the smallest is -1); used to validate and encode a/b tokens.
 * @param context Description of where the field comes from, used only in error messages.
 * @return The parsed change values, one per token.
 */
vector<int> parse_change_list(const string& s,
                              int is_total,
                              int max_haplotype_change,
                              const string& context);

/**
 * @brief Parse the site_change field of the meta file into one flat vector of changes.
 * Format: 1:0,1,-1;2:0,0,2;22:-1,1 (total CN) or 1:0/1,1/-1;2:0/0 (haplotype-specific CN).
 * Chromosome ids (1-22) are validated but not stored; the values of all blocks are concatenated.
 * @param s The field to parse; a missing field gives an empty vector.
 * @param is_total true for total copy number changes, false for haplotype-specific changes.
 * @param max_haplotype_change Largest change allowed on one haplotype (the smallest is -1).
 * @param context Description of where the field comes from, used only in error messages.
 * @return The site changes of all chromosomes, in input order.
 */
vector<int> parse_site_change_field(const string& s,
                                    bool is_total,
                                    int max_haplotype_change,
                                    const string& context);

/**
 * @brief Find the largest absolute chromosome-level and site-level change of each sample, e.g. from a meta file.
 * For haplotype-specific changes, this is the largest absolute change on either haplotype.
 * @param sample_change_chr Chromosome-level changes of each sample.
 * @param sample_change_site Site-level changes of each sample; must have as many samples as sample_change_chr.
 * @param chr_max_change Output: largest absolute chromosome-level change of each sample.
 * @param site_max_change Output: largest absolute site-level change of each sample.
 * @param states_chr Chromosome-level change states, used to decode haplotype-specific indices.
 * @param states_site Site-level change states, used to decode haplotype-specific indices.
 * @param is_total 1 for total copy number changes, 0 for haplotype-specific state indices.
 * @param debug Debug flag for verbose output.
 */
void compute_max_changes(const vector<vector<int>>& sample_change_chr,
                        const vector<vector<int>>& sample_change_site,
                        vector<int>& chr_max_change,
                        vector<int>& site_max_change,
                        const vector<pair<int,int>>& states_chr,
                        const vector<pair<int,int>>& states_site,
                        int is_total,
                        int debug);

/**
 *  @brief Read the input copy number file
 *  Total copy number rows: sample, chr, site, cn.
 *  Haplotype-specific rows: sample, chr, site, cnA, cnB, converted into a state index while reading.
 *  @param s_info A reference to a vector of vectors of vectors that will store the copy number information for each sample in the format s_info[sample - 1] = {vcn = {chr, sid, cn}}.
 *  @param filename The name of the input file containing copy number data.
 *  @param Ns The number of samples to read from the file.
 *  @param num_total_bins A reference to an integer that will store the total number of bins read from the file.
 *  @param cn_max The maximum allowed copy number.
 *  @param is_total An integer flag indicating whether the copy numbers are total (1) or haplotype-specific (0).
 *  @param is_rcn An integer flag indicating whether the copy numbers are relative copy numbers (1) or absolute copy numbers (0).
 *  @param debug An integer flag for enabling debug output.
 */
void read_cn(vector<vector<vector<int>>>& s_info, const string& filename, int Ns, int &num_total_bins, int cn_max, int is_total = 1, int is_rcn = 0, int debug = 0);


/******************* read input all in one ***********************/

/**
 * @brief Read the input copy numbers into a matrix with one row per segment; only used in testing ML tree building.
 * @param filename Input file name containing copy number data.
 * @param input_property Structure containing properties of the input data.
 * @param input_data Structure to store processed input data.
 * @param debug Debug flag for verbose output.
 * @return One row per segment: {chr, id_start, id_end, cn_sample1, ..., cn_sampleNs}
 */
vector<vector<int>> read_data_var_regions(const string& filename, const INPUT_PROPERTY& input_property, INPUT_DATA& input_data, int debug = 0);

/**
 *  @brief Read data and group regions by chromosome, used in tree building for more comprehensive data processing
 *  @param data: A map where each key is a chromosome number and the value is a vector of vectors containing copy number information for each segment on that chromosome.
 *  @param cn_file Input file name containing copy number data.
 *  @param seg_file Optional output file name to save segment information.
 *  @param input_property Structure containing properties of the input data.
 *  @param input_data Structure to store processed input data.
 *  @param debug Debug flag for verbose output.
*/
void read_data_var_regions_by_chr_state(map<int, vector<vector<int>>>& data, const string& cn_file, const string& seg_file, const INPUT_PROPERTY& input_property, INPUT_DATA& input_data, int debug);

/**
 *  @brief Read data and group regions by chromosome, used in tree building for more comprehensive data processing
 * @param data_change: A map where each key is a chromosome number and the value is a vector of vectors containing copy number change information for each segment on that chromosome.
 *  @param cn_file Input file name containing copy number data.
 *  @param seg_file Optional output file name to save segment information.
 *  @param meta_file Optional meta file with known WGDs and chromosome/site changes (see read_meta_info); empty to estimate them from the copy numbers.
 *  @param input_property Structure containing properties of the input data.
 *  @param input_data Structure to store processed input data, will be changed during processing.
 *  @param debug Debug flag for verbose output.
*/
void read_data_var_regions_by_chr_change(map<int, vector<vector<int>>>& data_change, const string& cn_file, const string& seg_file, const string& meta_file, const INPUT_PROPERTY& input_property, INPUT_DATA& input_data, int debug);

#endif
