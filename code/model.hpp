#ifndef MODEL_HPP
#define MODEL_HPP

// This file implements the model of copy number evolution

#include "common.hpp"   // for model enum

// using namespace std;  // must before r8lib.hpp


#include "matexp/matrix_exponential.hpp"
#include "matexp/r8lib.hpp"



typedef vector<double>::iterator DBIter;  // convenience typedefs
typedef pair<DBIter, DBIter> DBIterPair;



inline int compute_haplotype_change_dim(int n) {
    if (n % 2 != 0) {
        // n is odd: (2n + 1)^2
        int x = 2 * n + 1;
        return x * x;
    } else {
        // n is even: (2n)^2
        int x = 2 * n;
        return x * x;
    }
}

// to validate the rate matrix (row sum should be 0)
bool check_matrix_row_sum(double *mat, int nstate);

void check_pmats_blen2(int nstate, const vector<double>& blens, const vector<double*> pmat_per_blen);

// n = cn_max + 1 for model 1 (total copy number)
void get_transition_matrix_bounded(double* q, double* p, double t, int n);

// not used in practice to save effeorts in function call
double get_transition_prob_bounded(double* p, int sk, int sj, int n);


/*************** MK model *****************/
// Mk model based on JC model
inline double get_transition_prob(double mu, double blength, int sk, int sj){
  //if(debug) cout << "\tget_transition_prob" << endl;

  // assume a basic Markov model here cf JC
  // Qij = u/5 when i!=j
  // we have five states 0, 1, 2, 3, 4
  // Pij(t) = exp(-ut) k_ij + (1 - exp(-ut))*pi_j
  // pi_0 = pi_1 = pi_2 = pi_3 = pi_4 = 1/5

  double prob = 0.0;

  if(sk == sj){
    prob = exp(-mu * blength) + (1 - exp(-mu * blength)) * 0.2;
  }else{
    prob = (1 - exp(-mu * blength)) * 0.2;
  }

  return prob;
}

/*************** site-level BOUND model *****************/

// model when copy number is bounded by 0 and cn_max
void get_rate_matrix_bounded(double* m, double dup_rate, double del_rate, int cn_max);

// A matrix with dimension (cn_max + 1) * (cn_max + 2) / 2
// Suppose copy number configuration is in specific order such as:  0/0	0/1	1/0	0/2	 1/1	2/0	0/3	 1/2	 2/1	3/0	0/4	 1/3	  2/2	 3/1	4/0
void get_rate_matrix_haplotype_specific(double* m, double dup_rate, double del_rate, int cn_max);



/*************** multi-scale model *****************/

// model decomposing observed copy numbers into different levels

// rate matrix for segment-level CNAs
void get_rate_matrix_site_change(double* m, double dup_rate, double del_rate, int site_change_max);

// rate matrix for chr-level CNAs
void get_rate_matrix_chr_change(double* m, double chr_gain_rate, double chr_loss_rate, int chr_change_max);

// void get_rate_matrix_site_change_haplotype(double* m, double dup_rate, double del_rate,  int max_site_change_haplotype = 2);

// void get_rate_matrix_chr_change_haplotype(double* m, double chr_gain_rate, double chr_loss_rate,  int max_chr_change_haplotype = 1);

vector<Transition> build_transitions(const vector<pair<int,int>>& states, const map<pair<int,int>, int>& idx, int max_allowed_cn_change, double cn_increase_rate, double cn_decrease_rate);

void get_rate_matrix_change_haplotype(double* m, double cn_increase_rate, double cn_decrease_rate,  int max_change_haplotype = 1);


// Maximum allowed number of WGD events over a time interval: 2
void get_rate_matrix_wgd(double* m, double wgd_rate, int max_wgd = 2);



#endif
