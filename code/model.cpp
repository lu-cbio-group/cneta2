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


// rate matrix for total copy number 
void get_rate_matrix_bounded(double* m, const double dup_rate, const double del_rate, const int cn_max){
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
void get_rate_matrix_haplotype_specific(double* m, const double dup_rate, const double del_rate, const int cn_max){
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


// rate matrix for segment-level CNAs using independent Markov chain model
void get_rate_matrix_site_change(double* m, const double dup_rate, const double del_rate, const int site_change_max){
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
void get_rate_matrix_chr_change(double* m, const double chr_gain_rate, const double chr_loss_rate, const int chr_change_max){
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
void get_rate_matrix_site_change_haplotype(double* m, const double dup_rate, const double del_rate, const int site_change_max_haplotype){
    int debug = 0;
    int ncol = compute_haplotype_change_dim(site_change_max_haplotype);  
    assert(ncol == 16);     // hard coded for now, can be extended to higher copy number if needed
    int lcol = ncol - 1;    // last column

    for(unsigned i = 0; i < ncol; ++ i){
        for (unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }

    // the first two rows, after the first row with all 0s
    for(unsigned i = 1; i <= 2; i++){
        m[i] = del_rate;
        m[i + (i + 2) * ncol] = dup_rate;
        m[i + i * ncol] = 0 - del_rate - dup_rate;
    }
    for(unsigned i = 3; i <= 4; i++){
        m[i + (i - 2) * ncol] = del_rate;
        m[i + (i + 5) * ncol] = dup_rate;
        m[i + i * ncol] = 0 - del_rate - dup_rate;
    }

    for(unsigned i = 5; i <= 7; i++){
        m[i + i * ncol] = 0 - 2 * del_rate - 2 * dup_rate;
    }

    // for normal 1/1 state
    unsigned i = 5;
    m[i + (i - 3) * ncol] = del_rate;
    m[i + (i - 4) * ncol] = del_rate;
    m[i + (i + 1) * ncol] = dup_rate;
    m[i + (i + 2) * ncol] = dup_rate;  

    // for 1/2 and 2/1
    for(unsigned i = 6; i <= 7; i++){
        m[i + (i - 3) * ncol] = del_rate;
        m[i + (i + 5) * ncol] = dup_rate;
    }

    i = 6;
    m[i + (i - 3) * ncol] = del_rate;
    m[i + (i + 4) * ncol] = dup_rate;

    i = 7;
    m[i + (i - 2) * ncol] = del_rate;
    m[i + (i + 3) * ncol] = dup_rate;

    for(unsigned i = 8; i <= 9; i++){
        m[i + (i - 5) * ncol] = del_rate;
        m[i + i * ncol] = 0 - del_rate;
    }  

    // 2/2
    i = 10;
    m[i + (i - 3) * ncol] = del_rate;
    m[i + (i - 4) * ncol] = del_rate;
    m[i + (i + 3) * ncol] = dup_rate;
    m[i + (i + 4) * ncol] = dup_rate;  
    m[i + i * ncol] = 0 - 2 * del_rate - 2 * dup_rate;

    for(unsigned i = 11; i <= 12; i++){
        m[i + (i - 3) * ncol] = del_rate;
        m[i + (i - 5) * ncol] = del_rate;
        m[i + (i + 2) * ncol] = dup_rate;
        m[i + i * ncol] = 0 - 2 * del_rate - dup_rate;
    }   

    unsigned i_wgd = 10; // 2/2
    for(unsigned i = 13; i < lcol; i++){
        m[i + (i - 2) * ncol] = del_rate;
        m[i + i_wgd * ncol] = del_rate;
        m[i + lcol * ncol] = dup_rate;
        m[i + i * ncol] = 0 - 2 * del_rate - dup_rate;
    }   

    m[lcol + (lcol - 2)*ncol] = del_rate;
    m[lcol + (lcol - 1)*ncol] = del_rate;
    m[lcol + lcol*ncol] = 0 - 2 * del_rate;

    if(debug){
        std::cout << m << std::endl;
        r8mat_print(ncol, ncol, m, "  Haplotype-specific Site-level Q matrix:");
    }
}



// TODO: rate matrix for haplotype-specific chr-level CNAs using independent Markov chain model
// 0	1	2	3	4	5	6	7	8
// 0	1	1	2	2	2	3	3	4
// 0/0	0/1	1/0	0/2	 1/1	2/0	 1/2	 2/1	 2/2
// -2	-1	-1	0	0	0	+1	+1	+2
// -1/-1	-1/0	0/-1	-1/+1	0/0	+1/-1	0/+1	+1/0	+1/+1
// 9 states when at most one change allowed on one chr haplotype
// patterns not so obvious as the matrix for total changes
void get_rate_matrix_chr_change_haplotype(double* m, const double chr_gain_rate, const double chr_loss_rate, const int chr_change_max_haplotype){
    int debug = 0;
    int ncol = compute_haplotype_change_dim(chr_change_max_haplotype);  
    assert(ncol == 9);      
    int lcol = ncol - 1;    // last column

    for(unsigned i = 0; i < ncol; ++ i){
        for (unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }

    // the first two rows, after the first row with all 0s
    for(unsigned i = 1; i <= 2; i++){
        m[i] = chr_loss_rate;
        m[i + (i + 2) * ncol] = chr_gain_rate;
        m[i + i * ncol] = 0 - chr_loss_rate - chr_gain_rate;
    }

    // the next two rows, as first row has all 0s
    for(unsigned i = 3; i <= 4; i++){
        m[i + (i - 2) * ncol] = chr_loss_rate;
        m[i + i * ncol] = 0 - chr_loss_rate;
    }

    // for normal 1/1 state
    unsigned i = 5;
    m[i + (i - 3) * ncol] = chr_loss_rate;
    m[i + (i - 4) * ncol] = chr_loss_rate;
    m[i + (i + 1) * ncol] = chr_gain_rate;
    m[i + (i + 2) * ncol] = chr_gain_rate;
    m[i + i * ncol] = 0 - 2 * chr_loss_rate - 2 * chr_gain_rate;

    // the last two rows, as first row has all 0s
    unsigned i_norm = 5;
    for(unsigned i = 6; i < lcol; i++){
        m[i + (i - 3) * ncol] = chr_loss_rate;
        m[i + i_norm * ncol] = chr_loss_rate;
        m[i + lcol * ncol] = chr_gain_rate;
        m[i + i * ncol] = 0 - 2 * chr_loss_rate - chr_gain_rate;
    }   

    m[lcol + (lcol - 2)*ncol] = chr_loss_rate;
    m[lcol + (lcol - 1)*ncol] = chr_loss_rate;
    m[lcol + lcol*ncol] = 0 - 2 * chr_loss_rate;

    if(debug){
        std::cout << m << std::endl;
        r8mat_print(ncol, ncol, m, "  Haplotype-specific Chr-level Q matrix:");
    }
}


// Maximum allowed number of WGD events over a time interval: 2
void get_rate_matrix_wgd(double* m, const double wgd_rate, const int wgd_max){
    int debug = 0;
    int ncol = wgd_max + 1;

    for(unsigned i = 0; i < ncol; ++ i){
        for(unsigned j = 0; j < ncol; ++ j){
            m[i + j*ncol] = 0;
        }
    }
    for(unsigned i = 0; i < wgd_max; i++){
        m[i+(i + 1)*ncol] = wgd_rate;
        m[i+i*ncol] = 0 - wgd_rate;
    }

    if(debug){
        std::cout << m << std::endl;
        r8mat_print(ncol, ncol, m, "  WGD P matrix:");
    }
}


// // http://www.guwi17.de/ublas/examples/
// // Use eigendecomposition to compute matrix exponential
// double get_transition_prob_bounded(const boost::numeric::ublas::matrix<double>& q, double t, const int sk, const int sj ){
//     //
//     boost::numeric::ublas::matrix<double> tmp = q * t;
//     boost::numeric::ublas::matrix<double> p = expm_pad(tmp);
//     cout << "t: " << t << endl;
//     cout << "Q matrix: " << q << endl;
//     cout << "tmp matrix: " << tmp << endl;
//     cout << "P matrix: " << p << endl;
//     return p(sk, sj);
// }
// n = cn_max + 1 for model 1 (total copy number)
/**  
 * @brief Compute transition probability matrix from rate matrix          
 * @param q: rate matrix
 * @param p: transition probability matrix
 * @param t: branch length
 * @param n: dimension of the matrix
*/
void get_transition_matrix_bounded(double* q, double* p, const double t, const int n){
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
double get_transition_prob_bounded(double* p, const int sk, const int sj, const int n){
    int debug = 0;

    int i = sk  + sj * n;
    double v = p[i];

    if(debug){
        r8mat_print(n, n, p, "  P matrix before access:");
        cout << "Prob at " << i << "("<< sk << ", " << sj << ") is: " << v << endl;
    }

    return v;
}


/******************** building decomposition table for each possible copy number *********************/

// Insert a tuple of copy number changes with 3 parameters indicating WGD, chromosome gain/loss, and segment duplication/deletion
// not used due to too much simplification
void insert_tuple(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int i, int j, int k){
    int sum = pow(2, i + 1) + j + k;
    // cout << i << "\t" << j << "\t" << k << "\t" << sum << endl;
    if(sum >= 0 && sum <= cn_max){
        vector<int> c{i, j, k};
        decomp_table[sum].insert(c);
        comps.insert(c);
    }
}


// Insert a tuple of copy number changes with 3 parameters indicating changes WGD, chromosome gain/loss, and segment duplication/deletion
// Ignoring event order and only considering the change relative to normal copy number
void insert_tuple_change(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int i, int j, int k){
    int cn_change = pow(2, i + 1) + j + k;     // not as the observed total copy number, but change relative to normal copy number 2
    // cout << i << "\t" << j << "\t" << k << "\t" << sum << endl;
    if(cn_change + NORM_PLOIDY >= 0 && cn_change + NORM_PLOIDY <= cn_max){
        vector<int> c{i, j, k};
        decomp_table[cn_change + NORM_PLOIDY].insert(c);
        comps.insert(c);
    }else{
        cout << "Warning: observed copy number " << cn_change + NORM_PLOIDY << " is out of range [0, " << cn_max << "]" << endl;
        exit(EXIT_FAILURE);
    }
}



// Get possible combinations for a total copy number with ordering of different types of events considered (at most 1 WGD)
// Assuming at most 1 WGD event
void insert_tuple_order_withm(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int m1, int m2, int alpha, int beta, int gamma, int beta_before, int gamma_before){
    int sum = pow(2, alpha + 1) + m1 * beta + gamma + 2 * m2 * beta_before + 2 * gamma_before;
    // cout << alpha << "\t" << beta << "\t" << gamma << "\t" << sum << endl;
    if(sum >= 0 && sum <= cn_max){  // state out of range not considered
        // m1, m2 at the end to facilitate compability with other functions
        vector<int> c{alpha, beta, gamma, beta_before, gamma_before, m1, m2};
        decomp_table[sum].insert(c);
        comps.insert(c);
    }
}


// same to insert_tuple_order_withm when m_max = 1
void insert_tuple_order(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int i, int j, int k, int j0, int k0){
    // for(int m = 0; m <= m_max; m++){
    //     for(int m = 0; m <= m_max; m++){
            int sum = pow(2, i + 1) + j + k + 2 * j0 + 2 * k0;
            // cout << i << "\t" << j << "\t" << k << "\t" << sum << endl;
            if(sum >= 0 && sum <= cn_max){
                vector<int> c{i, j, k, j0, k0};
                decomp_table[sum].insert(c);
                comps.insert(c);
            }
        // }
    // }
    // // Cases when there are multiple copies of a segment on a chromosome before chromosome gain/loss
    // if(j > 0){
    //     for(int m = 2; m <= m_max; m++){
    //         sum = pow(2, i + 1) + m * j + k + 2 * j0 + 2 * k0;
    //         if(sum >= 0 && sum <= cn_max){
    //             vector<int> c{i, j, k, j0, k0};
    //             decomp_table[sum].insert(c);
    //             comps.insert(c);
    //         }
    //         if(j0 > 0){
    //             for(int m2 = 2; m2 <= m_max; m2++){
    //                 sum = pow(2, i + 1) + m * j + k + 2 * m2 * j0 + 2 * k0;
    //                 if(sum >= 0 && sum <= cn_max){
    //                     vector<int> c{i, j, k, j0, k0};
    //                     decomp_table[sum].insert(c);
    //                     comps.insert(c);
    //                 }
    //             }
    //         }
    //     }
    // }
    // if(j0 > 0){
    //     for(int m2 = 2; m2 <= m_max; m2++){
    //         sum = pow(2, i + 1) + j + k + 2 * m2 * j0 + 2 * k0;
    //         if(sum >= 0 && sum <= cn_max){
    //             vector<int> c{i, j, k, j0, k0};
    //             decomp_table[sum].insert(c);
    //             comps.insert(c);
    //         }
    //     }
    // }
}


// Get possible haplotype-specific combinations for a total copy number
void insert_tuple_haplotype_specific(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int i, int j1, int j2, int k1, int k2){
    int sum = pow(2, i + 1) + j1 + j2 + k1 + k2;
    // cout << i << "\t" << j1 << "\t" << j2 << "\t" << k1 << "\t" << k2 << "\t" << sum << endl;
    if(sum >= 0 && sum <= cn_max){
        vector<int> c{i, j1, j2, k1, k2};
        decomp_table[sum].insert(c);
        comps.insert(c);
    }
}


// not used for now for simplicity
void adjust_m_max(const vector<int>& sample_num_wgd, const vector<int>& sample_max_cn, int m_max, int max_chr_change, int max_site_change){
    // For samples with WGD but with larger copy number
    // max_sum = pow(2, 2) + 2 * m_max * max_chr_change + 2 * max_site_change;
    // m_max = m_max;
    // for(int i = 0; i < sample_max_cn.size(); i++){
    //     if(sample_num_wgd[i] < 1) continue;
    //     int cn = sample_max_cn[i];
    //     if(cn > max_sum){
    //         m_max = m_max + (cn - max_sum);
    //     }
    // }

    // For samples without WGD but with larger copy number
    int max_sum = pow(2, 1) + m_max * max_chr_change + max_site_change;
    int orig_m_max = m_max;
    int max_cn = max_sum;
    for(int i = 0; i < sample_max_cn.size(); i++){
        if(sample_num_wgd[i] > 0) continue;
        int cn = sample_max_cn[i];
        if(cn > max_cn){
            max_cn = cn;
        }
    }
    m_max = m_max + (max_cn - max_sum);
    int max_copy = 2 * (max_site_change + 1);    // at most 2 * (max_site_change + 1) copies of a segment before chromosome gain/loss
    m_max = m_max > max_copy ? max_copy : m_max;
    cout << "original m_max " << orig_m_max << endl;
    cout << "Max copy number for samples without WGD under original parameter setting is " << max_sum << endl;
    cout << "new m_max " << m_max << endl;
}


// list all the possible decomposition of a copy number, tuples with 3 elements: (num_WGD, chr_gain/loss, site_dup/del)
void build_decomp_table_3d(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int max_wgd, int max_chr_change, int max_site_change, int is_total){
    for(int i = 0; i <= cn_max; i++){
        set<vector<int>> comp;
        decomp_table[i] = comp;
    }
    if(is_total){
        // For 3D decomposition
        for(int i = 0; i <= max_wgd; i++){
            for(int j = 0; j <= max_chr_change; j++){
                for(int k = 0; k <= max_site_change; k++){
                    insert_tuple(decomp_table, comps, cn_max, 1, i, j, k);
                    insert_tuple(decomp_table, comps, cn_max, 1, i, -j, k);
                    insert_tuple(decomp_table, comps, cn_max, 1, i, j, -k);
                    insert_tuple(decomp_table, comps, cn_max, 1, i, -j, -k);
                }
            }
        }
    }
}


// list all the possible decomposition of a copy number, tuples with 5 elements: (num_WGD, chr_gain/loss, site_dup/del, chr_gain/loss_before_WGD, site_dup/del_before_WGD)
void build_decomp_table(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int max_wgd, int max_chr_change, int max_site_change, int is_total){
    for(int i = 0; i <= cn_max; i++){
        set<vector<int>> comp;
        decomp_table[i] = comp;
    }
    if(is_total){
        // For 3D decomposition
        // for(int i = 0; i <= max_wgd; i++){
        //     for(int j = 0; j <= max_chr_change; j++){
        //         for(int k = 0; k <= max_site_change; k++){
        //             insert_tuple(decomp_table, comps, cn_max, i, j, k);
        //             insert_tuple(decomp_table, comps, cn_max, i, -j, k);
        //             insert_tuple(decomp_table, comps, cn_max, i, j, -k);
        //             insert_tuple(decomp_table, comps, cn_max, i, -j, -k);
        //         }
        //     }
        // }
        // For 5D decomposition
        // No need to consider additional terms when no WGD occurs
        for(int j = 0; j <= max_chr_change; j++){
            for(int k = 0; k <= max_site_change; k++){
                insert_tuple_order(decomp_table, comps, cn_max, m_max, 0, j, k, 0, 0);
                insert_tuple_order(decomp_table, comps, cn_max, m_max, 0, -j, k, 0, 0);
                insert_tuple_order(decomp_table, comps, cn_max, m_max, 0, j, -k, 0, 0);
                insert_tuple_order(decomp_table, comps, cn_max, m_max, 0, -j, -k, 0, 0);
            }
        }

        // With WGD
        for(int i = 1; i <= max_wgd; i++){
            for(int j = 0; j <= max_chr_change; j++){
                for(int k = 0; k <= max_site_change; k++){
                    for(int j0 = 0; j0 <= max_chr_change - j; j0++){
                        for(int k0 = 0; k0 < max_site_change - k; k0++){    // remove = to disallow 0 copy before WGD
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, k, j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, k, -j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, k, j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, k, -j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, k, -j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, k, j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, k, -j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, k, j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, -k, j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, -k, -j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, -k, j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, j, -k, -j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, -k, j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, -k, -j0, k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, -k, j0, -k0);
                            insert_tuple_order(decomp_table, comps, cn_max, m_max, i, -j, -k, -j0, -k0);
                        }
                    }
                }
            }
        }
    }else{  // TO REVISE for haplotype-specific copy number
        for(int i = 0; i <= max_wgd; i++){
            for(int j1 = 0; j1 <= max_chr_change; j1++){
                for(int j2 = 0; j2 <= max_chr_change; j2++){
                    for(int k1 = 0; k1 <= max_site_change; k1++){
                        for(int k2 = 0; k2 <= max_site_change; k2++){
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, -k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, -k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, -k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, -k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, -k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, -k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, -k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, -k1, -k2);
                        }
                    }
                }
            }
        }
    }
}


// list all the possible decomposition of a copy number, tuples with 7 elements: (num_WGD, chr_gain/loss, site_dup/del, chr_gain/loss_before_WGD, site_dup/del_before_WGD, num_chr_copies_before_dup/del)
void build_decomp_table_withm(map<int, set<vector<int>>>& decomp_table, set<vector<int>>& comps, int cn_max, int m_max, int max_wgd, int max_chr_change, int max_site_change, int is_total){
    for(int i = 0; i <= cn_max; i++){
        set<vector<int>> comp;
        decomp_table[i] = comp;
    }
    if(is_total){
        int m1_max, m2_max;  // No need to consider m_max when there is no chromosome change
        // No need to consider additional terms when no WGD occurs
        for(int j = 0; j <= max_chr_change; j++){
            for(int k = 0; k <= max_site_change; k++){
                if(j == 0)  // no chromosome change, no need to consider segment copy number before change
                    m1_max = 0; 
                else 
                    m1_max = m_max;

                for(int m1 = 0; m1 <= m1_max; m1++){
                    // all possible combinations of chromosome change and segment change
                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, 0, 0, j, k, 0, 0);
                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, 0, 0, -j, k, 0, 0);
                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, 0, 0, j, -k, 0, 0);
                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, 0, 0, -j, -k, 0, 0);
                }
            }
        }

        // With WGD
        for(int i = 1; i <= max_wgd; i++){
            for(int j = 0; j <= max_chr_change; j++){
                for(int j0 = 0; j0 <= max_chr_change - j; j0++){    // -j to ensure total chromosome change within limit??
                    for(int k = 0; k <= max_site_change; k++){
                        for(int k0 = 0; k0 < max_site_change - k; k0++){
                            if(j == 0) 
                                m1_max = 0; 
                            else 
                                m1_max = m_max;
                            if(j0 == 0) 
                                m2_max = 0; 
                            else 
                                m2_max = m_max;
                            for(int m1 = 0; m1 <= m1_max; m1++){
                                for(int m2 = 0; m2 <= m2_max; m2++){
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, k, j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, k, -j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, k, j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, k, -j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, k, -j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, k, j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, k, -j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, k, j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, -k, j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, -k, -j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, -k, j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, j, -k, -j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, -k, j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, -k, -j0, k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, -k, j0, -k0);
                                    insert_tuple_order_withm(decomp_table, comps, cn_max, m_max, m1, m2, i, -j, -k, -j0, -k0);
                                }
                            }
                        }
                    }
                }
            }
        }
    }else{
        // TODO: to revise for haplotype-specific copy number
        for(int i = 0; i <= max_wgd; i++){
            for(int j1 = 0; j1 <= max_chr_change; j1++){
                for(int j2 = 0; j2 <= max_chr_change; j2++){
                    for(int k1 = 0; k1 <= max_site_change; k1++){
                        for(int k2 = 0; k2 <= max_site_change; k2++){
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, -k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, -k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, -k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, j2, -k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, -k1, k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, j2, -k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, j1, -j2, -k1, -k2);
                            insert_tuple_haplotype_specific(decomp_table, comps, cn_max, m_max, i, -j1, -j2, -k1, -k2);
                        }
                    }
                }
            }
        }
    }
}


void print_comps(const set<vector<int>>& comps){
    cout << "\nAll possible states for copy number decomposition in order:" << endl;
    int count = 0;
    for (auto c : comps){
        count += 1;
        cout << "\tstate " << count << endl;
        for(int i = 0; i < c.size(); i++){
            cout << "\t" << c[i];
        }
        cout << endl;
    }
}


void print_decomp_table(const map<int, set<vector<int>>>& decomp_table){
    int count = 0;
    for(auto item : decomp_table){
        set<vector<int>> comp = item.second;
        cout << "\nCopy number " << item.first << endl;
        for (auto c : comp){
            count += 1;
            for(int i = 0; i < c.size(); i++){
                cout << "\t" << c[i];
            }
            cout << endl;
        }
    }
    cout << "Number of possible states for copy number decomposition " << count << endl;
}
