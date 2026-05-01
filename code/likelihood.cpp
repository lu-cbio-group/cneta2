#include "likelihood.hpp"

#include "timing.hpp"
using timing::TimePoint;
using timing::now;
using timing::elapsed_seconds;

// global variables to track time spent in likelihood computation
static double time_revised_all   = 0.0;
static double time_revised_valid = 0.0;
static long   cnt_revised_all    = 0;
static long   cnt_revised_valid  = 0;

static double time_decomp_all    = 0.0;
static double time_decomp_valid  = 0.0;
static long   cnt_decomp_all     = 0;
static long   cnt_decomp_valid   = 0;



/**
 * @brief Print the likelihood table at the tips
 * @param rtree Evolutionary tree
 * @param obs Observed copy number states at tips
 * @param L_sk_k Likelihood table
 * @param nstate Number of possible states
 */
template <typename T>
void print_lnl_at_tips(const evo_tree& rtree, const vector<T>& obs, const vector<vector<double>>& L_sk_k, int nstate){
    cout << "\nObserved copy numbers at tips:\n";
    for(int i = 0; i < rtree.nleaf - 1; ++i){
        cout<< "\t" << obs[i];
    }
    cout << endl;

    cout << "Likelihood for tips:\n";
    for(int i = 0; i < rtree.nleaf; ++i){
        for(int j = 0; j < nstate; ++j){
            cout << "\t" << L_sk_k[i][j];
        }
        cout << endl;
    }
}


/**
 * @brief likelihood table initialization for models of site duplication/deletion only 
 * 
 * @param L_sk_k: Ns * nstate likelihood table 
 * @param obs: observed copy number states at tips
 * @param rtree: evolutionary tree 
 * @param model: evolutionary model used 
 * @param nstate: number of possible states       
 * @param is_total: whether or not the input is total copy number      
 */
void initialize_lnl_table(vector<vector<double>>& L_sk_k, const vector<int>& obs, const evo_tree& rtree, int model, int nstate, int is_total){
    int debug = 0;

    int Ns = rtree.nleaf - 1;
    if(model > 1){ // for BOUNDA model
        for(int i = 0; i < Ns; ++i){
            // For total copy number, all the possible combinations have to be considered for missing data.
            // Set related allele specific cases to be 1, with index from obs[i] * (obs[i] + 1)/2 to obs[i] * (obs[i] + 1)/2 + obs[i]. The index is computed based on pre-specified order.
            if(is_total){
                int si = (obs[i] * (obs[i] + 1)) / 2;
                int ei = si + obs[i];
                for(int k = si; k <= ei; k++){
                    L_sk_k[i][k] = 1.0;
                }
            }else{ // With haplotype-specific copy number, only the specific site needs to be filled
                L_sk_k[i][obs[i]] = 1.0;
            }
        }
        // set unaltered 1/1
        L_sk_k[Ns][NORM_ALLElE_STATE] = 1.0;
    }else{ // for MK and BOUNDT model
        for(int i = 0; i < Ns; ++i){
          for(int j = 0; j < nstate; ++j){
    	         if(j == obs[i]) L_sk_k[i][j] = 1.0;
          }
        }
        // set unaltered
        L_sk_k[Ns][NORM_PLOIDY] = 1.0;
    }

    if(debug){
      print_lnl_at_tips<int>(rtree, obs, L_sk_k, nstate);
    }
}


/** 
 *  @brief Get the likelihood on one site of a chromosome (assuming higher level events on nodes)
 *  TODO: likelihood assuming events order not well tested when WGD and/or chr-level changes occur
 *  @param L_sk_k Likelihood table
 *  @param rtree Evolutionary tree
 *  @param knodes List of internal nodes in preorder traversal
 *  @param blens List of unique branch lengths
 *  @param pmat_per_blen List of transition probability matrices corresponding to branch lengths, precomputed
 *  @param has_wgd Indicator for whole genome duplication event
 *  @param z Possible changes in copy number caused by chromosome gain/loss
 *  @param model Evolutionary model used
 *  @param nstate Number of possible states
 */
void get_likelihood_site(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const vector<int>& knodes, const vector<double>& blens, const vector<double*>& pmat_per_blen, int has_wgd, int z, int model, int nstate){
  int debug = 0;
  if(debug){
      cout << "Computing likelihood for one site" << endl;
  }

  for(int kn = 0; kn < knodes.size(); ++kn){
    int k = knodes[kn];
    int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
    double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
    int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
    double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

    // find index of bli and blj in blens
    auto pi = std::equal_range(blens.begin(), blens.end(), bli);
    // assert(distance(pi.first, pi.second) == 1);
    int idx_bli = std::distance(blens.begin(), pi.first);
    auto pj = std::equal_range(blens.begin(), blens.end(), blj);
    // assert(distance(pj.first, pj.second) == 1);
    int idx_blj = std::distance(blens.begin(), pj.first);

    double* pbli = pmat_per_blen[idx_bli];
    double* pblj = pmat_per_blen[idx_blj];

    if(debug){
      cout << "branch lengths so far:";
      for(auto b : blens){
          cout << "\t" << b;
      }
      cout << endl;
      cout << "node:" << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , " << bli << "\t" <<  nj + 1 << " , " << blj << endl;
      cout << "Get P matrix for branch length " << bli << "\t" << blens[idx_bli] << endl;
      r8mat_print(nstate, nstate, pmat_per_blen[idx_bli], "  P matrix:");
      cout << "Get P matrix for branch length " << blj << "\t" << blens[idx_blj] << endl;
      r8mat_print(nstate, nstate, pmat_per_blen[idx_blj], "  P matrix:");
    }

    //loop over possible values of sk
    if(k == rtree.nleaf){    // root node is always normal
        if(debug) cout << "Getting likelihood for root node " << k << endl;
        int nsk = NORM_PLOIDY;     // node sk
        if(model == BOUNDA) nsk = NORM_ALLElE_STATE;
        L_sk_k[k][nsk] = get_prob_children(L_sk_k, rtree, pbli, pblj, nsk, ni, nj, bli, blj, model, nstate);
    }else{
        for(int sk = 0; sk < nstate; ++sk){  // loop over possible states at internal node sk
            int nsk = sk;  // state after changes by other large scale events
            if(has_wgd) nsk = 2 * sk; 
            nsk += z; // incorporate chr gain/loss  
            if(debug) cout << "likelihood for state " << nsk << endl;
            if(nsk < 0 || nsk >= nstate) continue;
            // cout << " getting likelihood of children nodes " << endl;
            L_sk_k[k][nsk] = get_prob_children(L_sk_k, rtree, pbli, pblj, nsk, ni, nj, bli, blj, model, nstate);
        }
    }
  }
  if(debug){
    print_tree_lnl(rtree, L_sk_k, nstate);
  }
}



/** 
 * @brief Get the likelihood on the tree using matrix exponentiation for site duplication/deletion only
 *  used for bounded models BOUNDT and BOUNDA
 * @param vobs Observed copy number data
 * @param rtree Evolutionary tree
 * @param knodes List of internal nodes in preorder traversal
 * @param blens List of unique branch lengths
 * @param pmat_per_blen List of transition probability matrices corresponding to branch lengths, precomputed
 * @param has_wgd Indicator for whole genome duplication event
 * @param cn_type Copy number type indicator
 * @param use_repeat Indicator for reusing computed likelihood for repeated site patterns
 * @param model Evolutionary model used
 * @param nstate Number of possible states
 * @param is_total Indicator for total copy number or haplotype-specific copy number
 * @return double Log likelihood on the tree
 */
double get_likelihood_chr(const map<int, vector<vector<int>>>& vobs, const evo_tree& rtree, const vector<int>& knodes, const vector<double>& blens, const vector<double*>& pmat_per_blen, int has_wgd, int cn_type, int use_repeat, int model, int nstate, int is_total){
    int debug = 0;
    double logL = 0.0;    // for all chromosomes
    double chr_gain = 0.0;
    double chr_loss = 0.0;
    // Use a map to store computed log likelihood to save computation on duplicated site patterns
    map<vector<int>, vector<vector<double>>> sites_lnl_map;

    // for each chromosome
    for(auto vcn : vobs){
      int nchr = vcn.first;
      if(debug) cout << "Computing likelihood on Chr " << nchr << " with  " << vobs.at(nchr).size() << " sites " << endl;
      double chr_logL = 0.0;  // for one chromosome
      double chr_logL_normal = 0.0, chr_logL_gain = 0.0, chr_logL_loss = 0.0;   // when considering chr gain/loss
      double site_logL = 0.0;   // log likelihood for all sites on a chromosome
      int z = 0;    // no chr gain/loss
      // cout << " chromosome number change is " << 0 << endl;

      for(int nc = 0; nc < vobs.at(nchr).size(); nc++){    // for each segment on the chromosome
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs.at(nchr).at(nc);
          vector<vector<double>> L_sk_k(2 * rtree.nleaf - 1, vector<double>(nstate, 0.0));

          if(use_repeat){
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);
                  get_likelihood_site(L_sk_k, rtree, knodes, blens, pmat_per_blen, has_wgd, z, model, nstate);
                  sites_lnl_map[obs] = L_sk_k;
              }else{
                  // cout << "sites repeated" << end1;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);
              get_likelihood_site(L_sk_k, rtree, knodes, blens, pmat_per_blen, has_wgd, z, model, nstate);
          }

          site_logL += extract_tree_lnl(L_sk_k, rtree.nleaf - 1, model);

          if(debug){
              // cout << "\nLikelihood for site " << nc << " is " << lnl << endl;
              print_tree_lnl(rtree, L_sk_k, nstate);
          }
      }

      // likelihood calculation with chr gain/loss
      double chr_normal = 1.0;
      if(cn_type != ONLY_SEG){
          chr_gain = rtree.chr_gain_rate;
          chr_loss = rtree.chr_loss_rate;

          if(debug){
              cout << "chromosome gain rate " << chr_gain << endl;
              cout << "chromosome loss rate " << chr_loss << endl;
              cout << "Number of chr so far " << vobs.size() << endl;
          }

          if(fabs(chr_loss) > SMALL_VAL){
             chr_normal -= chr_loss;
          }
          if(fabs(chr_gain) > SMALL_VAL){
             chr_normal -= chr_gain;
          }
      }

      chr_logL_normal = log(chr_normal) + site_logL;
      chr_logL += chr_logL_normal;

      if(debug){
          cout << "Likelihood without chr gain/loss for " << nchr << " is "  << chr_normal << endl;
          cout << "Site Likelihood for " << nchr << " is "  << site_logL << endl;
          cout << "Likelihood without chr gain/loss: " << chr_logL_normal << endl;
      }

      if(cn_type != ONLY_SEG){
          if(fabs(chr_loss) > SMALL_VAL){
              z = -1;
              double site_logL = 0.0;   // log likelihood for all sites on a chromosome
              // cout << " chromosome number change is " << z << endl;
              for(int nc = 0; nc < vobs.at(nchr).size(); nc++){
                  // cout << "Number of sites for this chr " << vobs.at(nchr).size() << endl;
                  // for each site of the chromosome
                  vector<int> obs = vobs.at(nchr).at(nc);
                  vector<vector<double>> L_sk_k(2 * rtree.nleaf - 1, vector<double>(nstate, 0.0));
                  initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);

                  get_likelihood_site(L_sk_k, rtree, knodes, blens, pmat_per_blen, has_wgd, z, model, nstate);
                  site_logL += extract_tree_lnl(L_sk_k, rtree.nleaf - 1, model);

                  if(debug){
                      print_tree_lnl(rtree, L_sk_k, nstate);
                  }
              } // for all sites on a chromosome

              chr_logL_loss = log(chr_loss) + site_logL;
              chr_logL += log(1 + exp(chr_logL_loss - chr_logL_normal));
              if(debug){
                  cout << "\nLikelihood before chr loss for " << nchr << " is " << site_logL << endl;
                  cout << "\nLikelihood after chr loss: " << chr_logL_loss << endl;
              }
          } // for all chromosome loss

          if(fabs(chr_gain) > SMALL_VAL){
              z = 1;
              double site_logL = 0.0;   // log likelihood for all sites on a chromosome
              // cout << " chromosome number change is " << z << endl;
              for(int nc = 0; nc < vobs.at(nchr).size(); nc++){
                  // cout << "Number of sites for this chr " << vobs.at(nchr).size() << endl;
                  // for each site of the chromosome
                  vector<int> obs = vobs.at(nchr).at(nc);
                  vector<vector<double>> L_sk_k(2 * rtree.nleaf - 1, vector<double>(nstate, 0.0));
                  initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);
                  get_likelihood_site(L_sk_k, rtree, knodes, blens, pmat_per_blen, has_wgd, z, model, nstate);
                  site_logL += extract_tree_lnl(L_sk_k, rtree.nleaf - 1, model);

                  if(debug){
                      print_tree_lnl(rtree, L_sk_k, nstate);
                  }
              } // for all sites on a chromosome

              chr_logL_gain = log(chr_gain) + site_logL;
              if(chr_logL_loss > 0){
                  chr_logL += log(1 + 1 / (exp(chr_logL_normal - chr_logL_gain) + exp(chr_logL_loss - chr_logL_gain)));
              }else{
                  chr_logL += log(1 + exp(chr_logL_gain - chr_logL_normal));
              }

              if(debug){
                  cout << "\nLikelihood before chr gain for " << nchr << " is " << site_logL << endl;
                  cout << "\nLikelihood after chr gain: " << chr_logL_gain << endl;
              }
          } // for all chromosome loss
          // chr_logL = chr_logL_normal + log(1 + exp(chr_logL_loss-chr_logL_normal)) + log(1 + 1 / (exp(chr_logL_normal-chr_logL_gain) + exp(chr_logL_loss-chr_logL_gain)));
      }

      logL += chr_logL;

      if(debug){
          cout << "\nLikelihood after considering chr gain/loss for  " << nchr << " is " << logL << endl;
      }
    } // for each chromosome

    if(debug){
        cout << "\nLikelihood with chr gain/loss for all chromosmes: " << logL << endl;
    }

    return logL;
}




/** 
 *  @brief Get the likelihood on the tree using matrix exponentiation for site duplication/deletion only
 *  used for bounded models BOUNDT and BOUNDA
 *  @param rtree Evolutionary tree
 *  @param vobs Observed copy number data
 *  @param lnl_type Likelihood computation parameters
 */
double get_likelihood_revised(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, LNL_TYPE& lnl_type){
    // int debug = 0;
    // if(debug) cout << "\tget_likelihood by matrix exponential" << endl;

    TimePoint start_time = now();
    ++cnt_revised_all;

    if(!is_tree_valid(rtree, lnl_type.max_tobs, lnl_type.patient_age, lnl_type.cons)){// invalid tree only count toward "all" calls
        TimePoint end_all_invalid = now();
        time_revised_all += elapsed_seconds(start_time, end_all_invalid);
        return SMALL_LNL;
    }

    int model = lnl_type.model;
    int cn_max = lnl_type.cn_max;
    int is_total = lnl_type.is_total;

    // For copy number instantaneous changes
    int nstate = cn_max + 1;
    if(model == BOUNDA) nstate = (cn_max + 1) * (cn_max + 2) / 2;

    // may be different for each tree due to change of mutation rate
    int dim_mat = nstate * nstate;
    double *qmat = new double[dim_mat];
    memset(qmat, 0.0, dim_mat * sizeof(double));

    // model == BOUNDT || model == BOUNDA
    assert(model > 0);
    if(model == BOUNDA){
        get_rate_matrix_haplotype_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
    }else{
        get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
    }

    // if(debug){
    //   r8mat_print(nstate, nstate, qmat, "  Q matrix:" );
    //   check_matrix_row_sum(qmat, nstate);
    // }

    // Find the transition probability matrix for each branch at once to save computation for branches with same lengths
    vector<int> knodes = lnl_type.knodes; // internal nodes in preorder traversal, including the root
    vector<double> blens;
    vector<double*> pmat_per_blen;    // transition probability matrices for each unique branch length, each element is a nstate * nstate matrix

    for(int kn = 0; kn < knodes.size(); ++kn){
        int k = knodes[kn];
        // two outgoing branches of the current internal node k
        double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
        double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

        if(find(blens.begin(), blens.end(), bli) == blens.end()){
            double *pmati = new double[dim_mat];
            memset(pmati, 0.0, dim_mat * sizeof(double));
            get_transition_matrix_bounded(qmat, pmati, bli, nstate);
            pmat_per_blen.push_back(pmati);
            blens.push_back(bli);
            // if(debug){
            //     cout << "Get Pmatrix for branch length " << bli << endl;
            //     r8mat_print(nstate, nstate, pmati, "  P matrix:" );
            // }
        }
        if(find(blens.begin(), blens.end(), blj) == blens.end()){
            double *pmatj = new double[dim_mat];
            memset(pmatj, 0.0, dim_mat * sizeof(double));
            get_transition_matrix_bounded(qmat, pmatj, blj, nstate);
            pmat_per_blen.push_back(pmatj);
            blens.push_back(blj);
            // if(debug){
            //     cout << "Get Pmatrix for branch length " << blj << endl;
            //     r8mat_print(nstate, nstate, pmatj, "  P matrix:" );
            // }
        }
    }

    // sort pmats according to branch lengths
    auto p = sort_permutation(blens, [&](const double a, const double b){ return a < b; });
    blens = apply_permutation(blens, p);
    pmat_per_blen = apply_permutation(pmat_per_blen, p);

    // if(debug){
    //     for(int i = 0; i < pmat_per_blen.size(); ++i){
    //         double blen = blens[i];
    //         cout << "Get P matrix for branch length " << blen << endl;
    //         r8mat_print(nstate, nstate, pmat_per_blen[i], "  P matrix:");
    //     }
    // }

    double logL = 0.0;

    if(lnl_type.cn_type){         // only consider segment-level copy number changes
        // if(debug) cout << "Computing the likelihood without consideration of WGD" << endl;
        logL += get_likelihood_chr(vobs, rtree, knodes, blens, pmat_per_blen, 0, lnl_type.cn_type, lnl_type.use_repeat, model, nstate, is_total);
    }else{
        // if(debug) cout << "Computing the likelihood with consideration of WGD" << endl;
        logL += (1 - rtree.wgd_rate) * get_likelihood_chr(vobs, rtree, knodes, blens, pmat_per_blen, 0, lnl_type.cn_type, lnl_type.use_repeat, model, nstate, is_total);
        logL += rtree.wgd_rate * get_likelihood_chr(vobs, rtree, knodes, blens, pmat_per_blen, 1, lnl_type.cn_type, lnl_type.use_repeat, model, nstate, is_total);
    }

    // if(debug) cout << "Final likelihood before correcting acquisition bias: " << logL << endl;

    if(lnl_type.correct_bias){
        // if(debug) cout << "Correcting for the skip of invariant sites" << endl;

        // Compute the likelihood of dummy sites consisting entirely of 2s for the tree
        double lnl_invar = 0.0;
        // Suppose the value is 2 for all samples
        int normal_cn = 2;
        if(!is_total){
            normal_cn = 4;    // state ID for haplotype-specific CN 1/1
        }

        vector<int> obs(rtree.nleaf - 1, normal_cn);
        vector<vector<double>> L_sk_k(2 * rtree.nleaf - 1, vector<double>(nstate, 0.0));
        initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);

        for(int kn = 0; kn < knodes.size(); ++kn){
            int k = knodes[kn];
            int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
            double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
            int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
            double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

            DBIterPair pi = equal_range(blens.begin(), blens.end(), bli);
            // assert(distance(pi.first, pi.second) == 1);
            int idx_bli = distance(blens.begin(), pi.first);
            DBIterPair pj = equal_range(blens.begin(), blens.end(), blj);
            // assert(distance(pj.first, pj.second) == 1);
            int idx_blj = distance(blens.begin(), pj.first);

            //loop over possible values of sk
            for(int sk = 0; sk < nstate; ++sk){
                double Li = 0.0;
                for(int si = 0; si < nstate; ++si){
                    if(model == MK){
                    Li += get_transition_prob(rtree.mu, bli, sk, si) * L_sk_k[ni][si];
                    }else{
                    Li += pmat_per_blen[idx_bli][sk + si * nstate] * L_sk_k[ni][si];
                    }
                }
                double Lj = 0.0;
                for(int sj = 0; sj < nstate; ++sj){
                    if(model == MK){
                        Lj += get_transition_prob(rtree.mu, blj, sk, sj) * L_sk_k[nj][sj];
                    }else{
                        Lj += pmat_per_blen[idx_blj][sk + sj * nstate] * L_sk_k[nj][sj];
                    }
                }

                L_sk_k[k][sk] = Li * Lj;
            }
        }

        lnl_invar = extract_tree_lnl(L_sk_k, rtree.nleaf - 1, model);
        cout << "Likelihood of an invariant bin: " << lnl_invar << endl;
        double bias = lnl_type.num_invar_bins * lnl_invar;
        cout << "Number of invariant bins " << lnl_type.num_invar_bins << endl;
        logL = logL + bias;
        cout << "Bias to correct " << bias << endl;
        cout << "Final likelihood after correcting acquisition bias: " << logL << endl;


        // if(debug){
        //     cout << "Likelihood of an invariant bin: " << lnl_invar << endl;
        //     cout << "Number of invariant bins " << lnl_type.num_invar_bins << endl;
        //     cout << "Bias to correct " << bias << endl;
        //     cout << "Final likelihood after correcting acquisition bias: " << logL << endl;
        // }
    }

    if(std::isnan(logL) || logL < SMALL_LNL) logL = SMALL_LNL;
    // if(debug){
    //     cout << "Final likelihood: " << logL << endl;
    // }

    delete [] qmat;
    for_each(pmat_per_blen.begin(), pmat_per_blen.end(), DeleteObject());

    // end time measurement, only count toward "valid" calls
    TimePoint end_all_valid = now();
    double dt_valid = elapsed_seconds(start_time, end_all_valid);
    
    time_revised_all   += dt_valid; // valid also count toward all, so time for all is always updated
    time_revised_valid += dt_valid; // only update time for valid calls
    ++cnt_revised_valid; // only update count for valid calls

    return logL;
}



/** 
 *  @brief Compute the likelihood without grouping sites by chromosome, only considering site duplication/deletion
 *  The initial version with more simplifications, used for testing
 *  @param vobs the observed data matrix
 *  @param rtree the given tree
 *  @param knodes list of internal nodes in preorder traversal
 *  @param model model of evolution
 *  @param cn_max maximum copy number
 *  @param is_total indicator of total copy number or haplotype-specific copy number
 */
double get_likelihood(const vector<vector<int>>& vobs, evo_tree& rtree, const vector<int>& knodes, int model, int cn_max, int is_total){
  int debug = 0;
  if(debug) cout << "\tget_likelihood" << endl;

  assert(vobs.size() > 1);
  int Nchar = vobs[0].size();  // Nchar: number of characters for each sample

  int nstate = cn_max + 1;
  if(model == BOUNDA) nstate = (cn_max + 1) * (cn_max + 2) / 2;

  int dim_mat = nstate * nstate;
  double *qmat = new double[dim_mat];
  memset(qmat, 0.0, dim_mat * sizeof(double));

  if(model > 0){
      if(debug){
          cout << "Getting rate matrix" << endl;
      }
      if(model == BOUNDA){
          get_rate_matrix_haplotype_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
      }else{
          get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
      }
  }

  double *pmati = new double[dim_mat];
  double *pmatj = new double[dim_mat];
  memset(pmati, 0.0, dim_mat * sizeof(double));
  memset(pmatj, 0.0, dim_mat * sizeof(double));

  double logL = 0;
  for(int nc = 0; nc < Nchar; ++nc){
    vector<int> obs = vobs[nc];
    if(debug){
      cout << "char: " << nc << endl;
      for(int i = 0; i < rtree.nleaf - 1; ++i) cout << "\t" << obs[i];
      cout << endl;
    }

    vector<vector<double>> L_sk_k(2 * rtree.nleaf - 1, vector<double>(nstate, 0.0));
    initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);
    if(debug){
        print_tree_lnl(rtree, L_sk_k, nstate);
    }

    for(int kn = 0; kn < knodes.size(); ++kn){
      int k = knodes[kn];
      int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
      double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
      int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
      double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;
      if(model > 0){
           get_transition_matrix_bounded(qmat, pmati, bli, nstate);
           get_transition_matrix_bounded(qmat, pmatj, blj, nstate);
           if(debug){
               cout << "Get Pmatrix for branch length " << bli << endl;
               r8mat_print(nstate, nstate, pmati, "  P matrix after change:" );
               cout << "Get Pmatrix for branch length " << blj << endl;
               r8mat_print(nstate, nstate, pmatj, "  P matrix after change:" );
           }
      }

      if(debug) cout << "node:" << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , " << bli << "\t" <<  nj + 1 << " , " << blj << endl;

      //loop over possible values of sk
      for(int sk = 0; sk < nstate; ++sk){
    	  double Li = 0;
    	  // loop over possible si
    	  for(int si = 0; si < nstate; ++si){
            if (model == MK){
                Li += get_transition_prob(rtree.mu, bli, sk, si) * L_sk_k[ni][si];
            }
            else{
                Li += pmati[sk + si * nstate] * L_sk_k[ni][si];
            }
    	       //cout << "\tscoring: Li\t" << li << "\t" << get_transition_prob(mu, bli, sk, si ) << "\t" << L_sk_k[ni][si] << endl;
        }

    	  double Lj = 0;
    	  // loop over possible sj
    	  for(int sj = 0; sj < nstate; ++sj){
            if (model == MK){
    	        Lj += get_transition_prob(rtree.mu, blj, sk, sj) * L_sk_k[nj][sj];
            }
            else{
                Lj += pmatj[sk + sj * nstate] * L_sk_k[nj][sj];
            }
    	  }

	      //cout << "scoring: sk" << sk << "\t" <<  Li << "\t" << Lj << endl;
	      L_sk_k[k][sk] = Li * Lj;
      }

      if(debug){
    	  print_tree_lnl(rtree, L_sk_k, nstate);
      }
    }

    logL += extract_tree_lnl(L_sk_k, rtree.nleaf - 1, model);
  }

  if(debug) cout << "Final likelihood: " << logL << endl;

  delete [] qmat;
  delete [] pmati;
  delete [] pmatj;

  return logL;
}


/** 
 *  @brief Extract the log likelihood from the likelihood table at the root
 *  @param L_sk_k Likelihood table
 *  @param Ns Number of samples
 *  @param model Evolutionary model used, based on observed copy number values
 */
double extract_tree_lnl(vector<vector<double>>& L_sk_k, int Ns, int model){
    int debug = 0;

    if(debug){
        for(int j = 0; j < L_sk_k[Ns + 1].size(); ++j){
          cout << "\t" << L_sk_k[Ns + 1][j];
        }
        cout << endl;
    }

    if(model == BOUNDA){
        // The index is changed from 2 to 4 (1/1)
        if(debug) cout << "Likelihood for root is " << L_sk_k[Ns + 1][NORM_ALLElE_STATE] << endl;

        if(L_sk_k[Ns + 1][NORM_ALLElE_STATE] > 0) return log(L_sk_k[Ns + 1][NORM_ALLElE_STATE]);
        else return LARGE_LNL;
    }else{  // BOUNDT or MK
        if(debug) cout << "Likelihood for root is " << L_sk_k[Ns + 1][NORM_PLOIDY] << endl;

        if(L_sk_k[Ns + 1][NORM_PLOIDY] > 0) return log(L_sk_k[Ns + 1][NORM_PLOIDY]);
        else return LARGE_LNL;
    }
}


/**
 * @brief Print the likelihood table on the tree
 * 
 * @param rtree 
 * @param L_sk_k 
 * @param nstate: total number of possible states on each node 
 */
void print_tree_lnl(const evo_tree& rtree, vector<vector<double>>& L_sk_k, int nstate){
    cout << "\nLikelihood so far:\n";

    int ntotn = 2 * rtree.nleaf - 1;
    for(int i = 0; i < ntotn; ++i){
        for(int j = 0; j < nstate; ++j){
          cout << "\t" << L_sk_k[i][j];
        }
        cout << endl;
    }
}


/************* Funtions when decomposing observed copy numbers onto multiple level changes **********************/
// cn_total: copy number change state plus 2 to make it non-negative
// get the start and end index for the observed total copy number change in the rate matrix
void get_state_index(int cn_total, int peak_sum_haplotype, int& si, int& ei){
    // number of elements for each state
    vector<int> n_states = make_peak_vector(peak_sum_haplotype); // assume 2 in total and 1 for each haplotype at chr level
    int last_index = cn_total;

    if(cn_total < 0){
        cout << "Total copy number is negative, resetting to 0" << endl;
        last_index = 0;
    }    

    int max_tcn = 2 * peak_sum_haplotype - 2;
    if(cn_total > max_tcn){
        cout << "Total copy number " << cn_total << " exceeds the allowed maximum value, resetting to maximum valid value" << endl;
        last_index = max_tcn;
    }
    for (int k = 0; k < last_index; k++){
        si += n_states[k];
    }
    ei = si + n_states[last_index] - 1;
}


// helper function to set likelihood table entries for different levels, considering the indexing of haplotype-specific states
// peak_sum_haplotype: 3 when at most +1 change for each haplotype, 4 when at most +2 change for each haplotype, and so on
void set_lnl_table_change(int cn_total, int row, int max_change_haplotype, vector<vector<double>>& L_sk_k){     
    // only observed total copy number change
    int si = 0; // compute partial sum to get the index of the observed state
    int ei = 0;
    int peak_sum_haplotype = max_change_haplotype + 2;
    get_state_index(cn_total, peak_sum_haplotype, si, ei);

    for (int k = si; k <= ei; k++){
        L_sk_k[row][k] = 1.0;
    }
}


/**
 * @brief likelihood table initialization for models with WGD, chromosome gain/loss, and site duplication/deletion 
 * 
 * @param L_sk_k: Ns * nstate likelihood tables for WGD, chromosome gain/loss, and site duplication/deletion       
 * @param obs_change: observed changes in copy number states at tips
 * @param rtree: evolutionary tree 
 * @param nstates: number of states for WGD, chromosome gain/loss, and site duplication/deletion
 * @param dim_decomp: dimensions for different levels
 * @param debug: debug flag 
 */
void initialize_lnl_table_chr(vector<vector<double>>& lnl_table_chr, const evo_tree& rtree, const vector<int>& change_chr, const DIM_DECOMP& dim_decomp, const LNL_TYPE& lnl_type, int debug){
    if(debug) cout << "\tinitialize likelihood table at tips for chromosome gain/loss" << endl;
    int Ns = rtree.nleaf - 1;

    assert(dim_decomp.dim_chr > 0);

    lnl_table_chr.resize(2 * Ns + 1, vector<double>(dim_decomp.dim_chr, 0.0));
    lnl_table_chr[Ns][NO_CHANGE_HAPLOTYPE] = 1.0;  // no chr gain/loss for haplotype-specific case in normal sample

    for(int i = 0; i < Ns; ++i){
        // For chromosome gain/loss       
        // chr state is likely -2, -1, 0, 1, -1 for total copy number, 
        // and -1/-1	-1/0	0/-1	-1/+1	+1/-1	0/0	0/+1	+1/0	+1/+1 for haplotype-specific copy number
        // TODO: As the matrix is indexed by haplotype-specific changes, all related state need to be set.
        if(lnl_type.is_total){
            if(debug) cout << "Observed total copy number change for sample " << i + 1 << " is " << change_chr[i] << endl;
            int state_chr = 0;             // can lose at most two copies
            // check ambuigous encoding, only minor effect on likelihood
            if(change_chr[i] != 0 && change_chr[i] % CHANGE_CHR == 0){
                if(debug) cout << "Ambiguous encoding for chromosome gain/loss, including cases either before or after WGD." << endl;
                // after WGD, -2 was normalized to -1
                state_chr = change_chr[i] / CHANGE_CHR - MIN_CHANGE;    // CN normalied by ploidy previously
                set_lnl_table_change(state_chr, i, lnl_type.max_chr_change_haplotype, lnl_table_chr);
                // before WGD, -2
                state_chr = (state_chr + MIN_CHANGE) * NORM_PLOIDY - MIN_CHANGE;   // add WGD back
                set_lnl_table_change(state_chr, i, lnl_type.max_chr_change_haplotype, lnl_table_chr);            
            }else{
                state_chr = change_chr[i] - MIN_CHANGE;        
                set_lnl_table_change(state_chr, i, lnl_type.max_chr_change_haplotype, lnl_table_chr);    
            }            
        }else{
            if(debug) cout << "Observed haplotype-specific copy number change state for sample " << i + 1 << " is " << change_chr[i] << endl;
            // TODO after increasing rate matrix
            // after WGD: -2/0 was normalized to -1/0, or -2/-2 to -1/-1
            if(change_chr[i] % CHANGE_CHR == 0){                
                if(debug) cout << "Ambiguous encoding for chromosome gain/loss, including cases either before or after WGD." << endl;
                // before WGD, -2 was normalized to -1, +2 was normalized to +1
                int state_chr = change_chr[i] / CHANGE_CHR;    // CN normalied by ploidy previously
                lnl_table_chr[i][state_chr] = 1.0; 
                // after WGD, +2, 2 concecutive changes, not supported for now
                // need to distinguish changes on haplotypes to get exact haplotype-specific copy numbers
            }else{    
                assert(change_chr[i] >= 0 && change_chr[i] < dim_decomp.dim_chr);   // can gain or lose at most one copy for each haplotype
                // haplotype-specific observed change, the index for rate matrix       
                lnl_table_chr[i][change_chr[i]] = 1.0;  
            }
        }

    } 
    if(debug){
        cout << "Likelihood table at tips for chromosome gain/loss:";
        print_lnl_at_tips<int>(rtree, change_chr, lnl_table_chr, dim_decomp.dim_chr);
    }
}



/** 
 * @brief Compute the likelihood at a child node given the state at the parent node in a decomposed model 
 * 
 * @param node: index of the child node
 * @param cn_change_chr: index for chromosome gain/loss states at the parent node
 * @param lnl_table_chr: likelihood table for each node and each possible combination of states
 * @param pbl_chr: transition probabilities for different levels  
 * @param dim_decomp: dimensions of different levels  
 * @param debug: debug flag
 * @return doulbe: likelihood value
 * When the rate matrix is defined by haplotype-specific changes, sk include the index for haplotype-specific changes
 */
double compute_child_likelihood_chr(int node, int cn_change_chr, vector<vector<double>>& lnl_table_chr, const double* pbl_chr, const DIM_DECOMP& dim_decomp, int debug) {
    double L_chr = 0.0;

    // loop over all possible end states on a branch
    for (int e = 0; e < dim_decomp.dim_chr; ++e) {
        // as changes can be negative, need to transform to positive values first
        double p = pbl_chr[(cn_change_chr) + e * dim_decomp.dim_chr];
        L_chr += p * lnl_table_chr[node][e];
    }

    if (debug > 1) {
        std::cout << "\tLikelihood scoring for chr-level chagnes on node " << node << " given parent state "
                  << cn_change_chr << ", "
                  << L_chr << std::endl;
    }

    return L_chr;
}



/** 
 * @brief Get the likelihood on one site of a chromosome (assuming each observed copy number is composed of three type of events)   
 * 
 * @param lnl_table_chr: likelihood table for each node and each possible state at each copy number change level
 * @param rtree: evolutionary tree
 * @param knodes: list of internal nodes in post-order traversal
 * @param pmat_decomp: transition probability matrices for different levels
 * @param dim_decomp: dimensions of different levels
 * @param debug: debug flag
 */
// Sum over all possible states for internal nodes
// only need to calculate WGD likelihood for one site, as it is the same for all sites in the sample
void get_likelihood_per_chr(vector<vector<double>>& lnl_table_chr, const evo_tree& rtree, const vector<int>& knodes, const PMAT_DECOMP& pmat_decomp, const DIM_DECOMP& dim_decomp, int debug){
  if(debug){
      cout << "Computing likelihood for one chromosome with rate matrix dimensions: " << dim_decomp.dim_chr << endl;
  }

  assert(dim_decomp.dim_chr > 0);

  for(size_t kn = 0; kn < knodes.size(); ++kn){
        int k = knodes[kn];
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

      
        map<double, double*> pmats_chr = pmat_decomp.pmats_chr;
        double* pbli_chr = pmats_chr[bli];
        double* pblj_chr = pmats_chr[blj];
        
        if(debug > 1) cout << "node: " << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , " << bli << "\t" <<  nj + 1 << " , " << blj << endl;

        // loop over possible states of internal nodes
        if(k == rtree.nleaf){    // root node is always normal
            if(debug > 1) cout << "Getting likelihood for root node " << k << endl;
            int sk = NO_CHANGE_HAPLOTYPE;   
            // int sk_idx = change_to_state(sk, cn_max);
            double li_vals = compute_child_likelihood_chr(ni, sk, lnl_table_chr, pbli_chr, dim_decomp, debug);
            double lj_vals = compute_child_likelihood_chr(nj, sk, lnl_table_chr, pblj_chr, dim_decomp, debug);
            double lnl_chr = li_vals * lj_vals;  
            lnl_table_chr[k][sk] = lnl_chr;
        }else{
            if(debug > 1) cout << "Getting likelihood for internal node " << k << endl;

            // fill all possible haplotype-specific states 
            for(int sk = 0; sk < dim_decomp.dim_chr; ++sk){
                // use haplotype-specific changes for placeholder
                if(debug > 1) cout << "Computing likelihood for chromosome change " << sk << endl;   
                double li_vals = compute_child_likelihood_chr(ni, sk, lnl_table_chr, pbli_chr, dim_decomp, debug);
                double lj_vals = compute_child_likelihood_chr(nj, sk, lnl_table_chr, pblj_chr, dim_decomp, debug);
                double lnl_chr = li_vals * lj_vals;  
                lnl_table_chr[k][sk] = lnl_chr;                      
            }     
        }
  }
  
  if(debug > 0){
    cout << "The likelihood tables after get_likelihood_per_chr:" << endl;
    print_tree_lnl(rtree, lnl_table_chr, dim_decomp.dim_chr);
  }
}


/**
 * @brief likelihood table initialization for models with WGD, chromosome gain/loss, and site duplication/deletion 
 * 
 * @param lnl_table_wgd: Ns * nstate likelihood tables for WGD, chromosome gain/loss, and site duplication/deletion       
 * @param rtree: evolutionary tree 
 * @param sample_num_wgd: number of states for WGD across samples
 * @param dim_decomp: dimensions for different levels
 * @param debug: debug flag
 */
void initialize_lnl_table_wgd(vector<vector<double>>& lnl_table_wgd, const evo_tree& rtree, const vector<int>& sample_num_wgd, const DIM_DECOMP& dim_decomp, int debug){
    if(debug) cout << "\tinitialize likelihood table at tips for WGD" << endl;
    int Ns = rtree.nleaf - 1;

    assert(dim_decomp.dim_wgd > 0);         // number of WGD states  

    lnl_table_wgd.resize(2 * Ns + 1, vector<double>(dim_decomp.dim_wgd, 0.0));
    // for the normal cell 
    lnl_table_wgd[Ns][0] = 1.0;   // no WGD 

    for(int i = 0; i < Ns; ++i){
        // WGD state is 0,1,2
        int state_wgd = sample_num_wgd[i];
        lnl_table_wgd[i][state_wgd] = 1.0; 
    }   

    if(debug){
        cout << "Likelihood table at tips for WGD:";
        print_lnl_at_tips<int>(rtree, sample_num_wgd, lnl_table_wgd, dim_decomp.dim_wgd);
    }
}


/** 
 * @brief Compute the likelihood at a child node given the state at the parent node in a decomposed model 
 * 
 * @param node: index of the child node
 * @param num_wgd: index for WGD state at the parent node
 * @param lnl_table_wgd: likelihood table for each node and each possible combination of states
 * @param pbl_wgd: transition probabilities  
 * @param dim_decomp: dimensions of different levels  
 * @param debug: debug flag
 * @return double: likelihood value
 * When the rate matrix is defined by haplotype-specific changes, sk include the index for haplotype-specific changes
 */
double compute_child_likelihood_wgd(int node, int num_wgd, vector<vector<double>>& lnl_table_wgd, const double* pbl_wgd, const DIM_DECOMP& dim_decomp, int debug) {
    double L_wgd = 0.0;

    // loop over all possible end states on a branch
    for (int e = 0; e < dim_decomp.dim_wgd; ++e) {
        double p = pbl_wgd[num_wgd + e * dim_decomp.dim_wgd];
        L_wgd += p * lnl_table_wgd[node][e];
    }

    if (debug > 1) {
        std::cout << "\tLikelihood scoring for WGD on node " << node << " given parent state "
                  << num_wgd << ", "
                  << L_wgd << std::endl;
    }

    return L_wgd;
}


/** 
 * @brief Get the likelihood on one site of a chromosome (assuming each observed copy number is composed of three type of events)   
 * 
 * @param lnl_table_wgd: likelihood table for each node and each possible state at each copy number change level
 * @param rtree: evolutionary tree
 * @param knodes: list of internal nodes in post-order traversal
 * @param pmat_decomp: transition probability matrices for different levels
 * @param dim_decomp: dimensions of different levels
 */
// Sum over all possible states for internal nodes
// only need to calculate WGD likelihood for one site, as it is the same for all sites in the sample
void get_likelihood_wgd(vector<vector<double>>& lnl_table_wgd, const evo_tree& rtree, const vector<int>& knodes, const PMAT_DECOMP& pmat_decomp, const DIM_DECOMP& dim_decomp, int debug){
  if(debug){
      cout << "Computing likelihood for WGD with rate matrix dimensions: " << dim_decomp.dim_wgd << endl;
  }

  assert(dim_decomp.dim_wgd > 0);

  for(size_t kn = 0; kn < knodes.size(); ++kn){
        int k = knodes[kn];
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;
   
        map<double, double*> pmats_wgd = pmat_decomp.pmats_wgd;
        double* pbli_wgd = pmats_wgd[bli];
        double* pblj_wgd = pmats_wgd[blj];

        if(debug > 1) cout << "node: " << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , " << bli << "\t" <<  nj + 1 << " , " << blj << endl;

        // loop over possible states of internal nodes
        if(k == rtree.nleaf){    // root node is always normal
            if(debug > 1) cout << "Getting likelihood for root node " << k << endl;
            int sk = 0;   
            // int sk_idx = change_to_state(sk, cn_max);
            double li_vals = compute_child_likelihood_wgd(ni, sk, lnl_table_wgd, pbli_wgd, dim_decomp, debug);
            double lj_vals = compute_child_likelihood_wgd(nj, sk, lnl_table_wgd, pblj_wgd, dim_decomp, debug);
            double lnl_wgd = li_vals * lj_vals;     
            lnl_table_wgd[k][sk] = lnl_wgd;
        }else{
            if(debug > 1) cout << "Getting likelihood for internal node " << k << endl;

            for(int i = 0; i < dim_decomp.dim_wgd; ++i){   
                int sk = i;                 
                double li_vals = compute_child_likelihood_wgd(ni, sk, lnl_table_wgd, pbli_wgd, dim_decomp, debug);
                double lj_vals = compute_child_likelihood_wgd(nj, sk, lnl_table_wgd, pblj_wgd, dim_decomp, debug);
                double lnl_wgd = li_vals * lj_vals; 
                lnl_table_wgd[k][sk] = lnl_wgd;
                // cout << "WGD " << i << ": " << sk << " likelihood=" << lnl_wgd << endl;
            }      
        }
  }
  
  if(debug > 0){
    cout << "The likelihood tables after get_likelihood_wgd:" << endl;
    print_tree_lnl(rtree, lnl_table_wgd, dim_decomp.dim_wgd);
  }
}


/**
 * @brief likelihood table initialization for models with WGD, chromosome gain/loss, and site duplication/deletion 
 * 
 * @param lnl_table_seg: Ns * nstate likelihood tables for site duplication/deletion       
 * @param obs_change: observed changes in copy number states at tips
 * @param rtree: evolutionary tree 
 * @param dim_decomp: dimensions for different levels
 */
void initialize_lnl_table_site(vector<vector<double>>& lnl_table_seg, const evo_tree& rtree, const vector<CN_CHANGE>& obs_change, const DIM_DECOMP& dim_decomp, const LNL_TYPE& lnl_type, int debug){
    if(debug) cout << "\tinitialize likelihood table at tips for multiple levels" << endl;
    int Ns = rtree.nleaf - 1;

    // number of site duplication/deletion states, 16 for max CN 6 or 25 for max CN 8
    assert(dim_decomp.dim_seg > 0);

    lnl_table_seg.resize(2 * Ns + 1, vector<double>(dim_decomp.dim_seg, 0.0));
    lnl_table_seg[Ns][NO_CHANGE_HAPLOTYPE] = 1.0;  // no site gain/loss for haplotype-specific case

    for(int i = 0; i < Ns; ++i){       
        // For site duplication/deletion
        if(lnl_type.is_total){
            if(debug) cout << "Observed total copy number change for sample " << i + 1 << " is " << obs_change[i].cn_change_site << endl;
            // site state is likely -2, -1, 0, 1, 2
            int state_site = obs_change[i].cn_change_site - MIN_CHANGE;                    // can lost at most two copies
            // TODO: check ambuigous encoding, a bit complicated due to chr-level changes
            set_lnl_table_change(state_site, i, lnl_type.max_site_change_haplotype,lnl_table_seg);             
        }else{
            if(debug) cout << "Observed haplotype-specific copy number change state for sample " << i + 1 << " is " << obs_change[i].cn_change_site << endl;
            assert(obs_change[i].cn_change_site >= 0 && obs_change[i].cn_change_site < dim_decomp.dim_seg);
            // haplotype-specific observed change, index for rate matrix
            lnl_table_seg[i][obs_change[i].cn_change_site] = 1.0; 
        }        
    }   

    if(debug){
        cout << "Likelihood table at tips for site duplication/deletion:";
        print_lnl_at_tips<CN_CHANGE>(rtree, obs_change, lnl_table_seg, dim_decomp.dim_seg);   
    }
}


/**
 * @brief Build rate matrices for different levels of copy number changes.
 * 
 * @param rtree: the evolutionary tree 
 * @param dim_decomp: dimensions of different levels    
 * @return QMAT_DECOMP: rate matrices for different levels 
 */
void build_rate_matrices(QMAT_DECOMP& qmat_decomp, const evo_tree& rtree, const DIM_DECOMP& dim_decomp, const LNL_TYPE& lnl_type, int debug){
    int dim_wgd = dim_decomp.dim_wgd;
    int dim_chr = dim_decomp.dim_chr;
    int dim_seg = dim_decomp.dim_seg;

    int dim_mat_wgd = dim_wgd * dim_wgd;
    int dim_mat_chr = dim_chr * dim_chr;
    int dim_mat_seg = dim_seg * dim_seg;

    
    double *qmat_wgd = nullptr;
    double *qmat_chr = nullptr;
    double *qmat_seg = nullptr;

    if(dim_wgd > 0){
        qmat_wgd = new double[dim_mat_wgd];  // WGD
        memset(qmat_wgd, 0.0, dim_mat_wgd * sizeof(double));
        get_rate_matrix_wgd(qmat_wgd, rtree.wgd_rate, lnl_type.max_wgd);
    }

    if(dim_chr > 0){
        qmat_chr = new double[dim_mat_chr];   // chromosome gain/loss
        memset(qmat_chr, 0.0, dim_mat_chr * sizeof(double));
        // get_rate_matrix_chr_change_haplotype(qmat_chr, rtree.chr_gain_rate, rtree.chr_loss_rate, lnl_type.max_chr_change_haplotype);
        get_rate_matrix_change_haplotype(qmat_chr, rtree.chr_gain_rate, rtree.chr_loss_rate, lnl_type.max_chr_change_haplotype);
    }

    if(dim_seg > 0){
        qmat_seg = new double[dim_mat_seg];  // site duplication/deletion
        memset(qmat_seg, 0.0, dim_mat_seg * sizeof(double));
        // get_rate_matrix_site_change_haplotype(qmat_seg, rtree.dup_rate, rtree.del_rate, lnl_type.max_site_change_haplotype);
        get_rate_matrix_change_haplotype(qmat_seg, rtree.dup_rate, rtree.del_rate, lnl_type.max_site_change_haplotype);
    }

    qmat_decomp.free_all();
    qmat_decomp.qmat_wgd = std::move(qmat_wgd);
    qmat_decomp.qmat_chr = std::move(qmat_chr);
    qmat_decomp.qmat_seg = std::move(qmat_seg);

    if(debug){
        r8mat_print(dim_wgd, dim_wgd, qmat_wgd, "  Q-WGD matrix:");
        r8mat_print(dim_chr, dim_chr, qmat_chr, "  Q-CHR matrix:");
        r8mat_print(dim_seg, dim_seg, qmat_seg, "  Q-SEG matrix:");
    }
}



/**
 * @brief Build transition probability matrices for different levels of copy number changes.
 * @param PMAT_DECOMP: transition probability matrices for different levels 
 * @param rtree: the evolutionary tree 
 * @param knodes: list of internal nodes in post-order traversal
 * @param qmat_decomp: rate matrices for different levels 
 * @param dim_decomp: dimensions of rate matrices at different levels 
 */
void build_transition_matrices(PMAT_DECOMP& pmat_decomp, const evo_tree& rtree, const vector<int>& knodes, const QMAT_DECOMP& qmat_decomp, const DIM_DECOMP& dim_decomp, int debug){
    double *qmat_wgd = qmat_decomp.qmat_wgd; 
    double *qmat_chr = qmat_decomp.qmat_chr; 
    double *qmat_seg = qmat_decomp.qmat_seg;

    int dim_wgd = dim_decomp.dim_wgd;
    int dim_chr = dim_decomp.dim_chr;
    int dim_seg = dim_decomp.dim_seg;

    int dim_mat_wgd = dim_wgd * dim_wgd;
    int dim_mat_chr = dim_chr * dim_chr;
    int dim_mat_seg = dim_seg * dim_seg;    

    
    // Find the transition probability matrix for each branch, indexed by branch length
    map<double, double*> pmats_wgd;
    map<double, double*> pmats_chr;
    map<double, double*> pmats_seg;

    double *pmati_wgd, *pmatj_wgd;
    double *pmati_chr, *pmatj_chr;
    double *pmati_seg, *pmatj_seg;

    for(size_t kn = 0; kn < knodes.size(); ++kn){
         int k = knodes[kn];
         double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
         double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

         // For WGD
         if(dim_wgd > 0){
             if(pmats_wgd.find(bli) == pmats_wgd.end()){
                pmati_wgd = new double[dim_mat_wgd];
                memset(pmati_wgd, 0.0, dim_mat_wgd*sizeof(double));
                get_transition_matrix_bounded(qmat_wgd, pmati_wgd, bli, dim_wgd);
                pmats_wgd[bli] = pmati_wgd;
             }
             if(pmats_wgd.find(blj) == pmats_wgd.end()){
                pmatj_wgd = new double[dim_mat_wgd];
                memset(pmatj_wgd, 0.0, dim_mat_wgd*sizeof(double));
                get_transition_matrix_bounded(qmat_wgd, pmatj_wgd, blj, dim_wgd);
                pmats_wgd[blj] = pmatj_wgd;
             }
         }

         // For chr gain/loss
         if(dim_chr > 0){
             if(pmats_chr.find(bli) == pmats_chr.end()){
                 pmati_chr = new double[dim_mat_chr];
                 memset(pmati_chr, 0.0, dim_mat_chr*sizeof(double));
                 get_transition_matrix_bounded(qmat_chr, pmati_chr, bli, dim_chr);
                 pmats_chr[bli] = pmati_chr;
             }
             if(pmats_chr.find(blj) == pmats_chr.end()){
                 pmatj_chr = new double[dim_mat_chr];
                 memset(pmatj_chr, 0.0, dim_mat_chr*sizeof(double));
                 get_transition_matrix_bounded(qmat_chr, pmatj_chr, blj, dim_chr);
                 pmats_chr[blj] = pmatj_chr;
             }
         }

         // For site duplication/deletion
         if(dim_seg > 0){
             if(pmats_seg.find(bli) == pmats_seg.end()){
                 pmati_seg = new double[dim_mat_seg];
                 memset(pmati_seg, 0.0, dim_mat_seg*sizeof(double));
                 get_transition_matrix_bounded(qmat_seg, pmati_seg, bli, dim_seg);
                 pmats_seg[bli] = pmati_seg;
             }
             if(pmats_seg.find(blj) == pmats_seg.end()){
                 pmatj_seg = new double[dim_mat_seg];
                 memset(pmatj_seg, 0.0, dim_mat_seg*sizeof(double));
                 get_transition_matrix_bounded(qmat_seg, pmatj_seg, blj, dim_seg);
                 pmats_seg[blj] = pmatj_seg;
            }
        }
   }

    if(debug){
      for(auto it = pmats_wgd.begin(); it != pmats_wgd.end(); ++it){
          double key = it->first;
          cout << "Get P matrix for branch length " << key << endl;
          r8mat_print(dim_wgd, dim_wgd, it->second, "  P-WGD matrix:");
      }
      for(auto it = pmats_chr.begin(); it != pmats_chr.end(); ++it){
          double key = it->first;
          cout << "Get P matrix for branch length " << key << endl;
          r8mat_print(dim_chr, dim_chr, it->second, "  P-CHR matrix:");
      }
      for(auto it = pmats_seg.begin(); it != pmats_seg.end(); ++it){
          double key = it->first;
          cout << "Get P matrix for branch length " << key << endl;
          r8mat_print(dim_seg, dim_seg, it->second, "  P-SEG matrix:");
      }
    }

    pmat_decomp.free_all();  // or equivalent manual loop
    pmat_decomp.pmats_wgd = std::move(pmats_wgd);
    pmat_decomp.pmats_chr = std::move(pmats_chr);
    pmat_decomp.pmats_seg = std::move(pmats_seg);
}



/** 
 * @brief Compute the likelihood at a child node given the state at the parent node in a decomposed model 
 * 
 * @param node: index of the child node
 * @param sk: starting index for copy number change state at the parent node
 * @param L_sk_k: likelihood table for each node and each possible combination of states
 * @param prob_decomp: transition probabilities for different levels  
 * @param dim_decomp: dimensions of different levels  
 * @param debug: debug flag
 * @return double: likelihood value
 * When the rate matrix is defined by haplotype-specific changes, sk include the index for haplotype-specific changes
 */
double compute_child_likelihood_change(int node, int cn_change_site, const vector<vector<double>>& lnl_table_seg, const PROB_DECOMP1& prob_decomp, const DIM_DECOMP& dim_decomp, int debug) {
    double L_seg = 0.0;

    for (int e = 0; e < dim_decomp.dim_seg; ++e) {
        double p = prob_decomp.pbl_seg[cn_change_site + e * dim_decomp.dim_seg];
        L_seg += p * lnl_table_seg[node][e];
    }

    if (debug > 1) {
        std::cout << "\tLikelihood scoring for node " << node << " given parent state ("
                  << cn_change_site << "): " 
                  << L_seg << std::endl;
    }

    return L_seg;
}



/** 
 * @brief Compute the probability of children nodes given the state at the parent node in a decomposed model 
 * @param L_sk_k: likelihood table for each node and each possible combination of states
 * @param rtree: evolutionary tree  
 * @param sk: state at the parent node         
 * @param prob_decompi, prob_decompj: transition probabilities for different branches to the two child nodes  
 * @param dim_decomp: dimensions of different levels  
 * @param ni: index of the first child node
 * @param nj: index of the second child node
 * @param debug: debug flag
 */
double get_prob_children_change(vector<vector<double>>& lnl_table_seg, const evo_tree& rtree, int cn_change_site, PROB_DECOMP1& prob_decompi, PROB_DECOMP1& prob_decompj, const DIM_DECOMP& dim_decomp, int ni, int nj, int debug){
    if(debug > 1) cout << "\tget_prob_children_change" << endl;
    
    double li_vals = compute_child_likelihood_change(ni, cn_change_site, lnl_table_seg, prob_decompi, dim_decomp, debug);
    
    double lj_vals = compute_child_likelihood_change(nj, cn_change_site, lnl_table_seg, prob_decompj, dim_decomp, debug);

    double result = li_vals * lj_vals;

    return result;
}



/** 
 * @brief Get the likelihood on one site of a chromosome (assuming each observed copy number is composed of three type of events)   
 * 
 * @param L_sk_k: likelihood table for each node and each possible state at each copy number change level
 * @param rtree: evolutionary tree
 * @param comps: set of possible decomposed component combinations
 * @param knodes: list of internal nodes in post-order traversal
 * @param pmat_decomp: transition probability matrices for different levels
 * @param dim_decomp: dimensions of different levels
 * @param debug: debug flag
 */
// Sum over all possible states for internal nodes
// only need to calculate WGD likelihood for one site, as it is the same for all sites in the sample
void get_likelihood_site_change(vector<vector<double>>& lnl_table_seg, const evo_tree& rtree, const vector<int>& knodes, const PMAT_DECOMP& pmat_decomp, const DIM_DECOMP& dim_decomp, const LNL_TYPE& lnl_type, int debug){
  if(debug){
      cout << "Computing likelihood for one site with rate matrix dimensions: " << dim_decomp.dim_wgd << "\t" << dim_decomp.dim_chr << "\t"  << dim_decomp.dim_seg << endl;
  }

  for(size_t kn = 0; kn < knodes.size(); ++kn){
        int k = knodes[kn];
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

        double *pbli_wgd, *pblj_wgd;
        double *pbli_chr, *pblj_chr;
        double *pbli_seg, *pblj_seg;

        if(dim_decomp.dim_wgd > 0){
            map<double, double*> pmats_wgd = pmat_decomp.pmats_wgd;
            pbli_wgd = pmats_wgd[bli];
            pblj_wgd = pmats_wgd[blj];
        }
        if(dim_decomp.dim_chr > 0){
            map<double, double*> pmats_chr = pmat_decomp.pmats_chr;
            pbli_chr = pmats_chr[bli];
            pblj_chr = pmats_chr[blj];
        }
        if(dim_decomp.dim_seg > 0){
            map<double, double*> pmats_seg = pmat_decomp.pmats_seg;
            pbli_seg = pmats_seg[bli];
            pblj_seg = pmats_seg[blj];
        }

        PROB_DECOMP1 prob_decompi;
        // prob_decompi.pbl_wgd = pbli_wgd;
        prob_decompi.pbl_chr = pbli_chr;
        prob_decompi.pbl_seg = pbli_seg;

        PROB_DECOMP1 prob_decompj;
        // prob_decompj.pbl_wgd = pblj_wgd;       
        prob_decompj.pbl_chr = pblj_chr;        
        prob_decompj.pbl_seg = pblj_seg;

        if(debug  > 1) cout << "node: " << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , " << bli << "\t" <<  nj + 1 << " , " << blj << endl;

        // loop over possible states of internal nodes
        if(k == rtree.nleaf){    // root node is always normal
            if(debug > 1) cout << "Getting likelihood for root node " << k << endl;
            int sk = NO_CHANGE_HAPLOTYPE;   // normal state in terms of decomposed copy number changes 
            // int sk_idx = change_to_state(sk, cn_max);
            double lnl_seg = get_prob_children_change(lnl_table_seg, rtree, sk, prob_decompi, prob_decompj, dim_decomp, ni, nj, debug);
            if(dim_decomp.dim_seg > 0){
                lnl_table_seg[k][NO_CHANGE_HAPLOTYPE] = lnl_seg;   
            }       
        }else{
            if(debug > 1) cout << "Getting likelihood for internal node " << k << endl;
            if(dim_decomp.dim_seg > 0){
                // fill all possible haplotype-specific states related to the observed total copy number change
                for(int sk = 0; sk < dim_decomp.dim_seg; ++sk){
                    // l is the index in the rate, so it may not be the real change, just for placeholder
                    if(debug > 1) cout << "state index for site change: " << sk << endl;
                    double lnl_seg = get_prob_children_change(lnl_table_seg, rtree, sk, prob_decompi, prob_decompj, dim_decomp, ni, nj, debug);
                    lnl_table_seg[k][sk] = lnl_seg;                
                }  
            }        
        }
  }
  
  if(debug > 0){
    cout << "The likelihood tables after get_likelihood_site_change:" << endl;
    print_tree_lnl(rtree, lnl_table_seg, dim_decomp.dim_seg);
  }
}


/** 
 *  @brief Extract the log likelihood from the likelihood table at the root for decomposed model
 *
 *  @param L_sk_k Likelihood tables for WGD, chromosome gain/loss, and site duplication/deletion
 *  @param Ns Number of samples
 */
// Get the likelihood of the tree from likelihood table of state combinations
double extract_tree_lnl_change(const vector<vector<double>>& lnl_table_seg, int Ns, int debug){
    if(debug) cout << "Extracting likelihood for the root" << endl;

    double log_lnl = 0.0;

    // Segment, no duplication/deletion
    if(!lnl_table_seg.empty()){
        double prob_seg = lnl_table_seg[Ns + 1][NO_CHANGE_HAPLOTYPE]; 
        if (prob_seg > 0.0) {
            log_lnl += log(prob_seg);
        } else {
            if (debug) cout << "Zero likelihood at root: segment change" << endl;
        } 
    } else {
        if (debug) cout << "No segment duplication/deletion events in the tree" << endl;
    }

    if(log_lnl == 0.0){
        log_lnl = LARGE_LNL;
    }

    if(debug) cout << "The final likelihood for the root: " << log_lnl << endl;

    return log_lnl;
}


/**
 * @brief Compute the likelihood by chromosome, assuming each observed copy number is decomposed into three types of events
 * @param rtree: the evolutionary tree 
 * @param vobs_change: observed copy number changes at different chromosomes, real change for total CN and state index for haplotype-specific CN
 * @param knodes: list of internal nodes in post-order traversal
 * @param pmat_decomp: transition probability matrices for different levels
 * @param dim_decomp: dimensions of different levels
 * @param lnl_type: tags when computing likelihood, such as indicator of whether to use repeated sites to speed up computation
 * @param is_total: whether to consider total copy number only
 * @return double: log likelihood value
 */
double get_likelihood_chr_change(const evo_tree& rtree, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const vector<int>& knodes, const PMAT_DECOMP& pmat_decomp, const DIM_DECOMP& dim_decomp, const LNL_TYPE& lnl_type, int debug){
    double logL = 0.0;    // for all chromosomes
    // recompute to save parameter change, also do sanity check   
    vector<int> sample_num_wgd(rtree.nleaf - 1, -1);         // use -1 to indicate uninitialized
    vector<vector<int>> chr_sample_change;   // indexed by chromosome first and then sample, to facilitate likelihood computation
    map<vector<CN_CHANGE>, vector<vector<double>>> sites_lnl_map;    // Use a map to store computed log likelihood
    map<vector<int>, vector<vector<double>>> chr_lnl_map;    // Use a map to store computed log likelihood

    for(auto vcn : vobs_change){
        int nchr = vcn.first;
        const vector<vector<CN_CHANGE>> obs_chr = vcn.second;   // all segments on this chromosome
        vector<int> change_chr(rtree.nleaf - 1, -100);   // total copy number change for each sample on this chromosome, use -100 to indicate uninitialized
        vector<vector<double>> lnl_table_chr;
        double site_logL = 0.0;   // log likelihood for all sites on a chromosome

        if(debug){
            cout << "Computing likelihood on Chr " << nchr <<  " with " << obs_chr.size() << " sites" << endl;
            // print_nested_vector<CN_CHANGE>(obs_chr);
        }
        
        for(size_t nc = 0; nc < obs_chr.size(); nc++){    // for each site on the chromosome (may be repeated)
            const vector<CN_CHANGE> obs = obs_chr.at(nc);
            vector<vector<double>> lnl_table_seg;     // one table for one site

            assert(obs.size() == rtree.nleaf - 1);  // sanity check, the number of samples should be the same as the number of tips in the tree
            for(size_t i = 0; i < obs.size(); i++){
                if(sample_num_wgd[i] == -1){
                    sample_num_wgd[i] = obs.at(i).num_wgd;
                }else{
                    assert(sample_num_wgd[i] == obs.at(i).num_wgd);  // sanity check, all sites should have the same WGD state for each sample
                }
                if(change_chr[i] == -100){
                    change_chr[i] = obs.at(i).cn_change_chr;
                }else{
                    assert(change_chr[i] == obs.at(i).cn_change_chr);  // sanity check, all sites should have the same WGD state for each sample
                }                
            }
            chr_sample_change.push_back(change_chr);

            if(debug > 1){
                cout << "Chr " << nchr << " Site " << nc << " observed changes: ";
                for(size_t i = 0; i < obs.size(); i++){
                    cout << "\t" << obs.at(i);
                }
                cout << endl;
            }

            if(dim_decomp.dim_seg > 0){
                if(debug){
                    cout << "Site changes observed, computing site_level likelihood for ";
                    print_vector<CN_CHANGE>(obs);
                }
                if(lnl_type.use_repeat){ 
                    if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                    if(debug) cout << "sites new" << endl;
                        initialize_lnl_table_site(lnl_table_seg, rtree, obs, dim_decomp, lnl_type, debug);
                        get_likelihood_site_change(lnl_table_seg, rtree, knodes, pmat_decomp, dim_decomp, lnl_type, debug);
                        sites_lnl_map[obs] = lnl_table_seg;
                    }else{
                        if(debug) cout << "sites repeated" << endl;
                        lnl_table_seg = sites_lnl_map[obs];
                    }
                }else{
                    if(debug) cout << "sites no repeat consideration" << endl;
                    initialize_lnl_table_site(lnl_table_seg, rtree, obs, dim_decomp, lnl_type, debug);
                    get_likelihood_site_change(lnl_table_seg, rtree, knodes, pmat_decomp, dim_decomp, lnl_type, debug);
                }

                site_logL += extract_tree_lnl_change(lnl_table_seg, rtree.nleaf - 1, debug);

                if(debug > 1){
                    cout << "\nLikelihood computed for site: ";
                    for(size_t i = 0; i < obs.size(); i++){
                        cout << "\t" << obs.at(i);
                    }
                    cout << endl;
                    print_tree_lnl(rtree, lnl_table_seg, dim_decomp.dim_seg);  
                }
            }
        }

        logL += site_logL;

        // compute likelihood for chr-level changes, which is the same for all sites on the chromosome, so only compute once
        if(dim_decomp.dim_chr > 0){
            if(debug){
                cout << "\nChr-level changes observed, computing chromosome-level likelihood for Chr " << nchr << endl;
            }
            if(lnl_type.use_repeat){ 
                if(chr_lnl_map.find(change_chr) == chr_lnl_map.end()){
                    if(debug){
                        cout << "chr sites new on chr " << nchr << endl;
                        print_vector<int>(change_chr);
                    }
                    initialize_lnl_table_chr(lnl_table_chr, rtree, change_chr, dim_decomp, lnl_type, debug);
                    get_likelihood_per_chr(lnl_table_chr, rtree, knodes, pmat_decomp, dim_decomp, debug);
                    chr_lnl_map[change_chr] = lnl_table_chr;
                }else{
                    if(debug){
                        cout << "chr sites repeated on chr " << nchr << endl;
                        print_vector<int>(change_chr);
                    }
                    lnl_table_chr = chr_lnl_map[change_chr];
                }
            }else{
                if(debug){
                    cout << "chr sites no repeat consideration" << endl;
                    print_vector<int>(change_chr);
                }
                initialize_lnl_table_chr(lnl_table_chr, rtree, change_chr, dim_decomp, lnl_type, debug);
                get_likelihood_per_chr(lnl_table_chr, rtree, knodes, pmat_decomp, dim_decomp, debug);
            }

            double chr_logL = log(lnl_table_chr[rtree.nleaf][NO_CHANGE_HAPLOTYPE]);   // likelihood for no change on chromosome, which is the same for all sites on the chromosome
            logL += chr_logL;   

            if(debug){
                cout << "\nLikelihood for chromosome-level changes on " << nchr << " is " << chr_logL << endl;
                cout << "\nSite-level likelihood for chromosome " << nchr << " is " << site_logL << endl;
            }                     
         }
    } // for each chromosome

    // only need to compute once for WGD as it is the same for all sites in the sample
    if(dim_decomp.dim_wgd > 0){
        vector<vector<double>> lnl_table_wgd;
        initialize_lnl_table_wgd(lnl_table_wgd, rtree, sample_num_wgd, dim_decomp, debug);    
        get_likelihood_wgd(lnl_table_wgd, rtree, knodes, pmat_decomp, dim_decomp, debug);  
        double wgd_logL = log(lnl_table_wgd[rtree.nleaf][0]);
        logL += wgd_logL;   // likelihood for WGD is the same for all sites, so only compute once

        if(debug){
            cout << "\nWGD observed, computed WGD likelihood for all samples  " << wgd_logL << endl;
        }       
    }

    if(debug){
        cout << "\nTotal likelihood is " << logL << endl;
    }

    return logL;
}



/**
 * @brief Computing likelihood when WGD and chr gain/loss are encoded by copy number changes.
 * 
 * @param rtree: the evolutionary tree 
 * @param obs: observed copy number changes at the tips 
 * @param obs_decomp: decomposition of observed copy numbers
 * @param lnl_type: tags when computing likelihood, such as indicator of whether to use repeated sites to speed up computation
 * @param debug: debug flag
 */
double get_likelihood_change(evo_tree& rtree, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const LNL_TYPE& lnl_type, int debug){
    if(debug) cout << "\tget likelihood using multiple chains based on copy number changes" << endl;

    // start calculating time
    TimePoint start_time = now();
    ++cnt_decomp_all;

    if(!is_tree_valid(rtree, lnl_type.max_tobs, lnl_type.patient_age, lnl_type.cons)){
        // invalid tree only count toward "all" calls
        TimePoint end_all_invalid = now();
        time_decomp_all += elapsed_seconds(start_time, end_all_invalid);

        return SMALL_LNL;
    }

    int model = lnl_type.model;
    int is_total = lnl_type.is_total;

    int dim_wgd = 0;
    int dim_chr = 0;
    int dim_seg = 0;

    // only consider states when changes are observed 
    if(obs_decomp.max_wgd > 0) dim_wgd = lnl_type.max_wgd + 1;
    // TODO:may consider gain/loss separately and arm-level changes in the future
    // based on maximum haplotype-specific change: +1: 9; +2: 16; +3: 25: +4: 36
    if(obs_decomp.max_chr_change > 0){      
        dim_chr = get_matrix_dim(lnl_type.max_chr_change_haplotype);
    }
    if(obs_decomp.max_site_change > 0){
        dim_seg = get_matrix_dim(lnl_type.max_site_change_haplotype);
    }

    int dim_mat_wgd = dim_wgd * dim_wgd;
    int dim_mat_chr = dim_chr * dim_chr;
    int dim_mat_seg = dim_seg * dim_seg;  

    DIM_DECOMP dim_decomp;
    dim_decomp.dim_wgd = dim_wgd;
    dim_decomp.dim_chr = dim_chr;
    dim_decomp.dim_seg = dim_seg;

    if(debug){
        cout << "Dimensions of rate/transition matrices for WGD, chr gain/loss, site gain/loss: " << dim_wgd << "\t" << dim_chr << "\t" << dim_seg << endl;
        cout << "\tBuilding Q rate matrices for multiple levels" << endl;
    }

    QMAT_DECOMP qmat_decomp;
    build_rate_matrices(qmat_decomp, rtree, dim_decomp, lnl_type, debug);

    if(debug) cout << "\tBuilding P transition matrices for multiple levels" << endl;
    PMAT_DECOMP pmat_decomp;
    build_transition_matrices(pmat_decomp, rtree, lnl_type.knodes, qmat_decomp, dim_decomp, debug);

    double logL = get_likelihood_chr_change(rtree, vobs_change, lnl_type.knodes, pmat_decomp, dim_decomp, lnl_type, debug);

    if(debug) cout << "Final likelihood before correcting acquisition bias: " << logL << endl;
    if(lnl_type.correct_bias){
        if(debug) cout << "Correcting for the skip of invariant sites" << endl;          
        vector<vector<double>> lnl_table_seg;
        CN_CHANGE invariant{2, 0, 0, 0};
        if (!is_total) {
            invariant = {2, 0, NO_CHANGE_HAPLOTYPE, NO_CHANGE_HAPLOTYPE};
        }
        vector<CN_CHANGE> obs(rtree.nleaf - 1, invariant);
        initialize_lnl_table_site(lnl_table_seg, rtree, obs, dim_decomp, lnl_type, debug);
        get_likelihood_site_change(lnl_table_seg, rtree, lnl_type.knodes, pmat_decomp, dim_decomp, lnl_type, debug);

        double lnl_invar = extract_tree_lnl_change(lnl_table_seg, rtree.nleaf - 1, debug);

        double bias = lnl_type.num_invar_bins * lnl_invar;
        logL = logL + bias;

        if(debug){
            cout << "Likelihood of an invariant bin " << lnl_invar << endl;
            cout << "Number of invariant bins " << lnl_type.num_invar_bins << endl;
            cout << "Bias to correct " << bias << endl;
            cout << "Final likelihood after correcting acquisition bias: " << logL << endl;
        }
    }

    if(std::isnan(logL) || logL < SMALL_LNL) logL = SMALL_LNL;
    if(debug){
        cout << "Final likelihood: " << logL << endl;
        cout << "Free memory" << endl;
    }

    TimePoint end_all_valid = now();
    double dt_valid = elapsed_seconds(start_time, end_all_valid);
    
    time_decomp_all += dt_valid; // valid also count toward all, so time for all is always updated
    time_decomp_valid += dt_valid; // only update time for valid calls
    ++cnt_decomp_valid; // only update count for valid calls

    if(debug > 1)  cout << "Likelihood calculation time for this tree: " << dt_valid << " seconds." << endl;

    return logL;
}



/**
 * @brief TODO: Computing likelihood when WGD and chr gain/loss are encoded by copy number changes assuming variable rates across branches
 * 
 * @param rtree: the evolutionary tree 
 * @param obs: observed copy number changes at the tips 
 * @param obs_decomp: decomposition of observed copy numbers
 * @param lnl_type: tags when computing likelihood, such as indicator of whether to use repeated sites to speed up computation
 * @param debug: debug flag
 */
double get_likelihood_change_variable_rate(evo_tree& rtree, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const OBS_DECOMP& obs_decomp, const LNL_TYPE& lnl_type, int debug){
    if(debug) cout << "\tget likelihood using multiple chains based on copy number changes" << endl;

    // start calculating time
    TimePoint start_time = now();
    ++cnt_decomp_all;

    if(!is_tree_valid(rtree, lnl_type.max_tobs, lnl_type.patient_age, lnl_type.cons)){
        // invalid tree only count toward "all" calls
        TimePoint end_all_invalid = now();
        time_decomp_all += elapsed_seconds(start_time, end_all_invalid);

        return SMALL_LNL;
    }

    int model = lnl_type.model;
    int cn_max = lnl_type.cn_max;
    int is_total = lnl_type.is_total;

    int dim_wgd = 0;
    int dim_chr = 0;
    int dim_seg = 0;

    // only consider states that can reach observed changes
    if(obs_decomp.max_wgd > 0) dim_wgd = lnl_type.max_wgd + 1;
    // TODO:may consider gain/loss separately and arm-level changes in the future
    // based on maximum haplotype-specific change: +1: 9; +2: 16; +3: 25: +4: 36
    if(obs_decomp.max_chr_change > 0){      
        // dim_chr = 9; 
        dim_chr = get_matrix_dim(lnl_type.max_chr_change_haplotype);
    }
    if(obs_decomp.max_site_change > 0){
        // dim_seg = 16;
        dim_seg = get_matrix_dim(lnl_type.max_site_change_haplotype);
    }

    int dim_mat_wgd = dim_wgd * dim_wgd;
    int dim_mat_chr = dim_chr * dim_chr;
    int dim_mat_seg = dim_seg * dim_seg;  

    DIM_DECOMP dim_decomp;
    dim_decomp.dim_wgd = dim_wgd;
    dim_decomp.dim_chr = dim_chr;
    dim_decomp.dim_seg = dim_seg;

    if(debug){
        cout << "Dimensions of rate/transition matrices for WGD, chr gain/loss, site gain/loss: " << dim_wgd << "\t" << dim_chr << "\t" << dim_seg << endl;
        cout << "\tBuilding Q rate matrices for multiple levels" << endl;
    }

    // TODO: 
    QMAT_DECOMP qmat_decomp;
    build_rate_matrices(qmat_decomp, rtree, dim_decomp, lnl_type, debug);

    if(debug) cout << "\tBuilding P transition matrices for multiple levels" << endl;
    PMAT_DECOMP pmat_decomp;
    build_transition_matrices(pmat_decomp, rtree, lnl_type.knodes, qmat_decomp, dim_decomp, debug);

    double logL = get_likelihood_chr_change(rtree, vobs_change, lnl_type.knodes, pmat_decomp, dim_decomp, lnl_type, debug);

    if(debug) cout << "Final likelihood before correcting acquisition bias: " << logL << endl;
    if(lnl_type.correct_bias){
        if(debug) cout << "Correcting for the skip of invariant sites" << endl;   

        CN_CHANGE invariant{2, 0, 0, 0};
        if (!lnl_type.is_total) {
            invariant = {2, 0, NO_CHANGE_HAPLOTYPE, NO_CHANGE_HAPLOTYPE};
        }
        vector<CN_CHANGE> obs(rtree.nleaf - 1, invariant);

        vector<vector<double>> lnl_table_seg;
        initialize_lnl_table_site(lnl_table_seg, rtree, obs, dim_decomp, lnl_type, debug);   
        get_likelihood_site_change(lnl_table_seg, rtree, lnl_type.knodes, pmat_decomp, dim_decomp, lnl_type, debug);
        double lnl_invar = extract_tree_lnl_change(lnl_table_seg, rtree.nleaf - 1, debug);

        double bias = lnl_type.num_invar_bins * lnl_invar;
        logL = logL + bias;

        if(debug){
            cout << "Likelihood of an invariant bin " << lnl_invar << endl;
            cout << "Number of invariant bins " << lnl_type.num_invar_bins << endl;
            cout << "Bias to correct " << bias << endl;
            cout << "Final likelihood after correcting acquisition bias: " << logL << endl;
        }
    }

    if(std::isnan(logL) || logL < SMALL_LNL) logL = SMALL_LNL;
    if(debug){
        cout << "Final likelihood: " << logL << endl;
        cout << "Free memory" << endl;
    }

    TimePoint end_all_valid = now();
    double dt_valid = elapsed_seconds(start_time, end_all_valid);
    
    time_decomp_all += dt_valid; // valid also count toward all, so time for all is always updated
    time_decomp_valid += dt_valid; // only update time for valid calls
    ++cnt_decomp_valid; // only update count for valid calls

    if(debug > 1)  cout << "Likelihood calculation time for this tree: " << dt_valid << " seconds." << endl;

    return logL;
}




