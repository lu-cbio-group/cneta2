#include "model.hpp"


// This file implements models of copy number evolution


bool check_matrix_row_sum(double *mat, int nstate){
    for(int i = 0; i < nstate; i++){
        double sum = 0;
        for(int j = 0; j < nstate; j++){    // jth column
            sum += mat[i + j * nstate];
        }
        cout << "Total probability of changing from " << i << " is " << sum << endl;
        if(sum != 0){
          return false;
        }
    }

    return true;
}


void check_pmats_blen2(int nstate, const vector<double>& blens, const vector<double*> pmat_per_blen){
  cout << "check P matrix" << endl;
  // find index of bli and blj in blens
  double bli = blens[1];
  double blj = blens[2];
  cout << bli << "\t" << blj << "\n";

  auto pi = equal_range(blens.begin(), blens.end(), bli);
  assert(distance(pi.first, pi.second) == 1);
  int idx_bli = std::distance(blens.begin(), pi.first);

  auto pj = std::equal_range(blens.begin(), blens.end(), blj);
  assert(distance(pj.first, pj.second) == 1);
  int idx_blj = std::distance(blens.begin(), pj.first);

  cout << idx_bli << "\t" << blens[idx_bli] << "\t" << bli << "\n";
  cout << idx_blj << "\t" << blens[idx_blj] << "\t" << blj << "\n";

  double* pbli = pmat_per_blen[idx_bli];
  double* pblj = pmat_per_blen[idx_blj];

  cout << "Print P matrix for branch length " << bli << "\t" << blens[idx_bli] << endl;
  r8mat_print(nstate, nstate, pmat_per_blen[idx_bli], "  P matrix:");
  cout << "Print P matrix for branch length " << blj << "\t" << blens[idx_blj] << endl;
  r8mat_print(nstate, nstate, pmat_per_blen[idx_blj], "  P matrix:");

}



/**  
 * @brief Compute transition probability matrix from rate matrix          
 * @param q: rate matrix
 * @param p: transition probability matrix
 * @param t: branch length
 * @param n: dimension of the matrix
*/
void get_transition_matrix_bounded(double* q, double* p, double t, int n){
    int debug = 0;

    double *tmp = new double[n*n];
    memset(tmp, 0.0, n*n*sizeof(double));
    for(int i = 0; i < n*n; i++){
        tmp[i] = q[i] * t;
    }

    double* res = r8mat_expm1(n, tmp);
    for(int i = 0; i < n*n; i++){
        if(res[i] < 0){
            p[i] = 0.0;
        }else{
            p[i] = res[i];
        }
    }

    if(debug){
        cout << "t: " << t << endl;
        r8mat_print(n, n, q, "  Q matrix:");
        r8mat_print(n, n, tmp, "  TMP matrix:");
        r8mat_print(n, n, p, "  P matrix:");
    }

    delete [] tmp;
    delete [] res;
}


// not used in practice to save efforts in function call
double get_transition_prob_bounded(double* p, int sk, int sj, int n){
    int debug = 0;

    int i = sk  + sj * n;
    double v = p[i];

    if(debug){
        r8mat_print(n, n, p, "  P matrix before access:");
        cout << "Prob at " << i << "("<< sk << ", " << sj << ") is: " << v << endl;
    }

    return v;
}


// boost::numeric::ublas::matrix<double> get_rate_matrix_bounded(double dup_rate = 0.01, double del_rate = 0.01){
//     boost::numeric::ublas::matrix<double> m(cn_max + 1, cn_max + 1);
//     for (unsigned i = 0; i < m.size1 (); ++ i){
//         for (unsigned j = 0; j < m.size2 (); ++ j){
//             m (i, j) = 0;
//         }
//     }
//     for( unsigned i = 1; i < cn_max; i++){
//         m (i,i-1) = 2 * i * del_rate;
//         m (i,i + 1) = 2 * i * dup_rate;
//         m(i, i) = 0 - m (i,i-1) - m (i,i + 1);
//     }
//     m (cn_max, cn_max - 1) = 2 * cn_max * del_rate;
//     m (cn_max, cn_max) = 0 - m (cn_max, cn_max - 1);
//
//     // std::cout << m << std::endl;
//     return m;
// }

/************************ functions for site-level bounded models **********************/

// rate matrix for total copy number 
void get_rate_matrix_bounded(double* m, double dup_rate, double del_rate, int cn_max){
    int debug = 0;   // for debugging internally
    int ncol = cn_max + 1;

    for (unsigned i = 0; i < ncol; ++ i){
        for (unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }
    for( unsigned i = 1; i < cn_max; i++){
        m[i+(i-1)*ncol] = 2 * i * del_rate;
        m[i+(i + 1)*ncol] = 2 * i * dup_rate;
        m[i+i*ncol] = 0 - m[i+(i-1)*ncol] - m[i+(i + 1)*ncol];
    }
    m[cn_max + (cn_max - 1)*ncol] = 2 * cn_max * del_rate;
    m[cn_max + cn_max*ncol] = 0 - m[cn_max + (cn_max - 1)*ncol];

    if(debug){
        std::cout << m << std::endl;
        r8mat_print(ncol, ncol, m, "  A:");
    }
}


// rate matrix for haplotype-specific CNAs
void get_rate_matrix_haplotype_specific(double* m, double dup_rate, double del_rate, int cn_max){
    int debug = 0;
    int ncol = (cn_max + 1) * (cn_max + 2) / 2;
    if(debug){
     cout << "Total number of states for copy number " << cn_max << " is " << ncol << endl;
    }

    for (unsigned i = 0; i < ncol; ++ i){
        for (unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }

    int s = cn_max * (cn_max - 1) / 2; // start index for the last group (total copy number 4)
    int r, c;   // row index and column index
    // 1st row for the last group (only deletion is allowed)
    r = ncol - (cn_max + 1);
    if(debug) cout << "Filling row " << r << endl;
    c = s;
    m[r + c*ncol] = del_rate;
    m[r + r*ncol] = 0 - del_rate;
    // last row for the last group
    r = ncol - 1 ;
    if(debug) cout << "Filling row " << r << endl;
    c = r - (cn_max + 1);
    m[r + c*ncol] = del_rate;
    m[r + r*ncol] = 0 - del_rate;
    // middle rows with two possibilities
    for(int j = cn_max; j > 1; j--){
        r = ncol - j;
        if(debug) cout << "Filling row " << r << endl;
        c = s;
        m[r + c*ncol] = del_rate;
        m[r + (c + 1)*ncol] = del_rate;
        m[r + r*ncol] = 0 - 2 * del_rate;

        s++;
    }

    int tcn = 0;    // total copy number
    // For entries with at least one zero, filling bottom up
    for(int i = 1; i < cn_max; i++){
        tcn = cn_max - i;
        int start_r = tcn * (tcn + 1)/2;

        // for zero at the front
        r = start_r;
        if(debug) cout << "Filling row " << r << endl;
        c = r - tcn;
        m[r + c*ncol] = del_rate;
        c = r + (tcn + 1);
        m[r + c*ncol] = dup_rate;
        m[r + r*ncol] = 0 - (del_rate + dup_rate);

        // for zero at the end
        r = start_r + tcn;
        if(debug) cout << "Filling row " << r << endl;
        c = r - (tcn + 1);
        m[r + c*ncol] = del_rate;
        c = r + (tcn + 2);
        m[r + c*ncol] = dup_rate;
        m[r + r*ncol] = 0 - (del_rate + dup_rate);
    }

    // For entries with none zero (from copy number 2 to cn_max-1, filling bottom up (excluding last group)
    for(int i = 1; i < cn_max - 1; i++){
        tcn = cn_max - i;
        int start_r = tcn * (tcn + 1)/2 + 1;
        for(int j = 1; j < tcn; j++){   // tcn - 2 rows to fill
            r = start_r + j - 1;
            if(debug) cout << "Filling row " << r << endl;

            c = r - (tcn + 1);
            m[r + c*ncol] = del_rate;
            c = c + 1;
            m[r + c*ncol] = del_rate;

            c = r + (tcn + 1);
            m[r + c*ncol] = dup_rate;
            m[r + (c + 1)*ncol] = dup_rate;

            m[r + r*ncol] = 0 - 2 * (del_rate + dup_rate);
        }
    }

    if(debug){
        r8mat_print(ncol, ncol, m, "  Haplotype-specific P matrix:");
    }
}


/************************ rate matrices for multi-scale models **********************/
// rate matrix for segment-level CNAs using independent Markov chain model
void get_rate_matrix_site_change(double* m, double dup_rate, double del_rate, int site_change_max){
    int debug = 0;
    int ncol = 2 * site_change_max + 1;
    int lcol = ncol - 1;    // last column

    for(unsigned i = 0; i < ncol; ++ i){
        for (unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }

    // start from 2nd row, index by row
    for(unsigned i = 1; i < ncol - 1; i++){
        m[i+(i - 1)*ncol] = del_rate;
        m[i+(i + 1)*ncol] = dup_rate;
        m[i+i*ncol] = 0 - m[i+(i-1)*ncol] - m[i+(i + 1)*ncol];
    }
    
    // last row
    m[lcol + (lcol - 1) * ncol] = del_rate;  // before the last element
    m[lcol + lcol*ncol] = 0 - del_rate;   // last element

    if(debug){
        std::cout << m << std::endl;
        r8mat_print(ncol, ncol, m, "  Site-level P matrix:");
    }
}



// rate matrix for chr-level CNAs using independent Markov chain model
void get_rate_matrix_chr_change(double* m, double chr_gain_rate, double chr_loss_rate, int chr_change_max){
    int debug = 0;
    int ncol = 2 * chr_change_max + 1;
    int lcol = ncol - 1;

    for(unsigned i = 0; i < ncol; ++ i){
        for (unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }

    // once lost all copies, cannot be regained
    for(unsigned i = 1; i < ncol - 1; i++){
        m[i+(i - 1)*ncol] = chr_loss_rate;
        m[i+(i + 1)*ncol] = chr_gain_rate;
        m[i+i*ncol] = 0 - chr_loss_rate - chr_gain_rate;
    }

    m[lcol + (lcol - 1)*ncol] = chr_loss_rate;
    m[lcol + lcol*ncol] = 0 - chr_loss_rate;

    if(debug){
        std::cout << m << std::endl;
        r8mat_print(ncol, ncol, m, "  Chr-level P matrix:");
    }
}



// TODO: rate matrix for haplotype-specific segment-level CNAs using independent Markov chain model
// 0	1	2	3	4	5	6	7	8	9	10	11	12	13	14	15
// 0	1	1	2	2	2	3	3	3	3	4	4	4	5	5	6
// 0/0	0/1	1/0	0/2	 1/1	2/0	0/3	 1/2	 2/1	3/0	 1/3	  2/2	 3/1	2/3	3/2	3/3
// -1/-1	-1/0	0/-1	-1/+1	0/0	+1/-1	-1/+2	0/+1	+1/0	+1/+1	0/+2	+1/+1	+2/0	+1/+2	+2/+1	+2/+2
// 16 states when at most one change allowed on one chr haplotype
// void get_rate_matrix_site_change_haplotype(double* m, double dup_rate, double del_rate, int max_site_change_haplotype){
//     int debug = 0;
//     int ncol = compute_haplotype_change_dim(max_site_change_haplotype);  
//     assert(ncol == 16);     // hard coded for now, can be extended to higher copy number if needed
//     int lcol = ncol - 1;    // last column

//     for(unsigned i = 0; i < ncol; ++ i){
//         for (unsigned j = 0; j < ncol; ++ j){
//             m[i + j*ncol] = 0;
//         }
//     }

//     // the first two rows, after the first row with all 0s
//     for(unsigned i = 1; i <= 2; i++){
//         m[i] = del_rate;
//         m[i + (i + 2) * ncol] = dup_rate;
//         m[i + i * ncol] = 0 - del_rate - dup_rate;
//     }
//     for(unsigned i = 3; i <= 4; i++){
//         m[i + (i - 2) * ncol] = del_rate;
//         m[i + (i + 5) * ncol] = dup_rate;
//         m[i + i * ncol] = 0 - del_rate - dup_rate;
//     }

//     for(unsigned i = 5; i <= 7; i++){
//         m[i + i * ncol] = 0 - 2 * del_rate - 2 * dup_rate;
//     }

//     // for normal 1/1 state
//     unsigned i = 5;
//     m[i + (i - 3) * ncol] = del_rate;
//     m[i + (i - 4) * ncol] = del_rate;
//     m[i + (i + 1) * ncol] = dup_rate;
//     m[i + (i + 2) * ncol] = dup_rate;  

//     // for 1/2 and 2/1
//     for(unsigned i = 6; i <= 7; i++){
//         m[i + (i - 3) * ncol] = del_rate;
//         m[i + (i + 5) * ncol] = dup_rate;
//     }

//     i = 6;
//     m[i + (i - 3) * ncol] = del_rate;
//     m[i + (i + 4) * ncol] = dup_rate;

//     i = 7;
//     m[i + (i - 2) * ncol] = del_rate;
//     m[i + (i + 3) * ncol] = dup_rate;

//     for(unsigned i = 8; i <= 9; i++){
//         m[i + (i - 5) * ncol] = del_rate;
//         m[i + i * ncol] = 0 - del_rate;
//     }  

//     // 2/2
//     i = 10;
//     m[i + (i - 3) * ncol] = del_rate;
//     m[i + (i - 4) * ncol] = del_rate;
//     m[i + (i + 3) * ncol] = dup_rate;
//     m[i + (i + 4) * ncol] = dup_rate;  
//     m[i + i * ncol] = 0 - 2 * del_rate - 2 * dup_rate;

//     for(unsigned i = 11; i <= 12; i++){
//         m[i + (i - 3) * ncol] = del_rate;
//         m[i + (i - 5) * ncol] = del_rate;
//         m[i + (i + 2) * ncol] = dup_rate;
//         m[i + i * ncol] = 0 - 2 * del_rate - dup_rate;
//     }   

//     unsigned i_wgd = 10; // 2/2
//     for(unsigned i = 13; i < lcol; i++){
//         m[i + (i - 2) * ncol] = del_rate;
//         m[i + i_wgd * ncol] = del_rate;
//         m[i + lcol * ncol] = dup_rate;
//         m[i + i * ncol] = 0 - 2 * del_rate - dup_rate;
//     }   

//     m[lcol + (lcol - 2)*ncol] = del_rate;
//     m[lcol + (lcol - 1)*ncol] = del_rate;
//     m[lcol + lcol*ncol] = 0 - 2 * del_rate;

//     if(debug){
//         std::cout << m << std::endl;
//         r8mat_print(ncol, ncol, m, "  Haplotype-specific Site-level Q matrix:");
//     }
// }



// TODO: rate matrix for haplotype-specific chr-level CNAs using independent Markov chain model
// 0	1	2	3	4	5	6	7	8
// 0	1	1	2	2	2	3	3	4
// 0/0	0/1	1/0	0/2	 1/1	2/0	 1/2	 2/1	 2/2
// -2	-1	-1	0	0	0	+1	+1	+2
// -1/-1	-1/0	0/-1	-1/+1	0/0	+1/-1	0/+1	+1/0	+1/+1
// 9 states when at most one change allowed on one chr haplotype
// patterns not so obvious as the matrix for total changes
// void get_rate_matrix_chr_change_haplotype(double* m, double chr_gain_rate, double chr_loss_rate, int max_chr_change_haplotype){
//     int debug = 0;
//     int ncol = compute_haplotype_change_dim(max_chr_change_haplotype);  
//     assert(ncol == 9);      
//     int lcol = ncol - 1;    // last column

//     for(unsigned i = 0; i < ncol; ++ i){
//         for (unsigned j = 0; j < ncol; ++ j){
//             m[i + j*ncol] = 0;
//         }
//     }

//     // the first two rows, after the first row with all 0s
//     for(unsigned i = 1; i <= 2; i++){
//         m[i] = chr_loss_rate;
//         m[i + (i + 2) * ncol] = chr_gain_rate;
//         m[i + i * ncol] = 0 - chr_loss_rate - chr_gain_rate;
//     }

//     // the next two rows, as first row has all 0s
//     for(unsigned i = 3; i <= 4; i++){
//         m[i + (i - 2) * ncol] = chr_loss_rate;
//         m[i + i * ncol] = 0 - chr_loss_rate;
//     }

//     // for normal 1/1 state
//     unsigned i = 5;
//     m[i + (i - 3) * ncol] = chr_loss_rate;
//     m[i + (i - 4) * ncol] = chr_loss_rate;
//     m[i + (i + 1) * ncol] = chr_gain_rate;
//     m[i + (i + 2) * ncol] = chr_gain_rate;
//     m[i + i * ncol] = 0 - 2 * chr_loss_rate - 2 * chr_gain_rate;

//     // the last two rows, as first row has all 0s
//     unsigned i_norm = 5;
//     for(unsigned i = 6; i < lcol; i++){
//         m[i + (i - 3) * ncol] = chr_loss_rate;
//         m[i + i_norm * ncol] = chr_loss_rate;
//         m[i + lcol * ncol] = chr_gain_rate;
//         m[i + i * ncol] = 0 - 2 * chr_loss_rate - chr_gain_rate;
//     }   

//     m[lcol + (lcol - 2)*ncol] = chr_loss_rate;
//     m[lcol + (lcol - 1)*ncol] = chr_loss_rate;
//     m[lcol + lcol*ncol] = 0 - 2 * chr_loss_rate;

//     if(debug){
//         std::cout << m << std::endl;
//         r8mat_print(ncol, ncol, m, "  Haplotype-specific Chr-level Q matrix:");
//     }
// }


//***  RL START ***//
// Discard the previous functions of hard coding the matrix, and use code to generate automatically.
// Build an universal function to get rate matrix for haplotype-specific CNAs, which can be used for both site-level and chr-level changes, by adjusting the input parameters and the dimension of the matrix. 
// The patterns of non-zero entries are similar for both types of changes, with only difference in the rates and the dimension of the matrix.

/**
 * @brief Build transitions for the CTMC based on the state space and rates.
 * 
 * @param states
 * @param idx 
 * @param max_allowed 
 * @param a: add rate, can be either site-level or chr-level change depending on the input parameters. 
                    * If site-level change, the add_rate is actually dup_rate; 
                    * if chr-level change, the add_rate is actually chr_gain_rate.
 * @param s: subtract rate, can be either site-level or chr-level change depending on the input parameters. 
                    * If site-level change, the subtract_rate is actually del_rate; 
                    * if chr-level change, the subtract_rate is actually chr_loss_rate.
 * @return * vector<Transition> 
 */
vector<Transition> build_transitions(const vector<pair<int,int>>& states,
                                     const map<pair<int,int>, int>& idx,
                                     int max_allowed_cn_change,
                                     double cn_increase_rate,
                                     double cn_decrease_rate) {
    vector<Transition> transitions;
    int max_hap_cn = NORM_PLOIDY / 2 + max_allowed_cn_change;  // maximum copy number on one haplotype, used to determine the state space

    for (size_t i = 0; i < states.size(); ++i) {
        int a = states[i].first;
        int b = states[i].second;
        int from = static_cast<int>(i);

        // subtract on The Pair a: (a,b) -> (a-1,b)
        if (a > 0) {
            pair<int,int> next = {a - 1, b};
            transitions.push_back({from, idx.at(next), cn_decrease_rate, "subtract_a"});
        }

        // add on The Pair a: (a,b) -> (a+1,b)
        // biological rule: a must already exist, so a > 0
        if (a > 0 && a < max_hap_cn && (a + 1 + b) <= max_allowed_cn_change) {
            pair<int,int> next = {a + 1, b};
            transitions.push_back({from, idx.at(next), cn_increase_rate, "add_a"});
        }

        // subtract on The Pair b: (a,b) -> (a,b-1)
        if (b > 0) {
            pair<int,int> next = {a, b - 1};
            transitions.push_back({from, idx.at(next), cn_decrease_rate, "subtract_b"});
        }

        // add on The Pair b: (a,b) -> (a,b+1)
        // biological rule: b must already exist, so b > 0
        if (b > 0 && b < max_hap_cn && (a + b + 1) <= max_allowed_cn_change) {
            pair<int,int> next = {a, b + 1};
            transitions.push_back({from, idx.at(next), cn_increase_rate, "add_b"});
        }
    }

    return transitions;
}

/**
 * @brief Get the rate matrix change haplotype object
 * 
 */
void get_rate_matrix_change_haplotype(double* m, double cn_increase_rate, double cn_decrease_rate,  int max_change_haplotype) {    
    int debug = 0;

    vector<pair<int,int>> states = build_pair_states_cn(max_change_haplotype);
    int n_states = states.size();
    map<pair<int,int>, int> idx = build_pair_index_map(states);
    vector<Transition> transitions = build_transitions(states, idx, max_change_haplotype, cn_increase_rate, cn_decrease_rate);

    // initialize the matrix with 0s
    for (int i = 0; i < n_states; ++i) {
    for (int j = 0; j < n_states; ++j) {
        m[i + j * n_states] = 0.0;
    }
}

    // off-diagonal entries
    for (const auto& tr : transitions) {
        m[tr.from + tr.to * n_states] += tr.rate;
    }

    // diagonal entries
    for (int i = 0; i < n_states; ++i) {
        double row_sum = 0.0;

        for (int j = 0; j < n_states; ++j) {
            if (i != j) {
                row_sum += m[i + j * n_states];
            }
        }

        m[i + i * n_states] = -row_sum;
    }

    if (debug) {
        std::cout << m << std::endl;
        r8mat_print(n_states, n_states, m, "  Haplotype-specific Q matrix:");
    }
}

//***  RL END ***//


// Maximum allowed number of WGD events over a time interval: 2
void get_rate_matrix_wgd(double* m, double wgd_rate, int max_wgd){
    int debug = 0;
    int ncol = max_wgd + 1;

    for(unsigned i = 0; i < ncol; ++ i){
        for(unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }
    for(unsigned i = 0; i < max_wgd; i++){
        m[i+(i + 1)*ncol] = wgd_rate;
        m[i+i*ncol] = 0 - wgd_rate;
    }

    if(debug){
        std::cout << m << std::endl;
        r8mat_print(ncol, ncol, m, "  WGD P matrix:");
    }
}





