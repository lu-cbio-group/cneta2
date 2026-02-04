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
void get_likelihood_site(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const vector<int>& knodes, const vector<double>& blens, const vector<double*>& pmat_per_blen, const int& has_wgd, const int& z, const int& model, const int& nstate){
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
double get_likelihood_chr(const map<int, vector<vector<int>>>& vobs, const evo_tree& rtree, const vector<int>& knodes, const vector<double>& blens, const vector<double*>& pmat_per_blen, const int& has_wgd, const int& cn_type, const int& use_repeat, const int& model, const int& nstate, const int& is_total){
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
      get_rate_matrix_allele_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
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
  auto p = sort_permutation(blens, [&](const double& a, const double& b){ return a < b; });
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
 * @brief Compute the likelihood without grouping sites by chromosome, only considering site duplication/deletion
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
          get_rate_matrix_allele_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
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

/**
 * @brief likelihood table initialization for models with WGD, chromosome gain/loss, and site duplication/deletion 
 * 
 * @param L_sk_k: Ns * nstate likelihood tables for WGD, chromosome gain/loss, and site duplication/deletion       
 * @param obs_change: observed changes in copy number states at tips
 * @param rtree: evolutionary tree 
 * @param nstates: number of states for WGD, chromosome gain/loss, and site duplication/deletion
 * @param dim_decomp: dimensions for different levels
 * @param obs_decomp: observed changes for different levels, including maximum changes 
 */
void initialize_lnl_table_change(LNL_TABLE& L_sk_k, const evo_tree& rtree, const vector<CN_CHANGE>& obs_change, const DIM_DECOMP& dim_decomp, const OBS_DECOMP& obs_decomp, int debug){
    if(debug) cout << "\tinitialize likelihood table at tips for multiple levels" << endl;
    int Ns = rtree.nleaf - 1;

    int nstate_wgd = dim_decomp.dim_wgd;         // number of WGD states  
    int nstate_chr = dim_decomp.dim_chr;         // number of chromosome gain/loss states
    int nstate_site = dim_decomp.dim_seg;        // number of site gain/loss

    // index shifts not used when using separate tables
    // int idx_chr_shift = 1 + nstate_wgd;
    // int idx_site_shift = 1 + nstate_wgd + nstate_chr;

    if(nstate_wgd > 0){
        L_sk_k.lnl_table_wgd.resize(2 * Ns + 1, vector<double>(nstate_wgd, 0.0));
        // for the normal cell 
        L_sk_k.lnl_table_wgd[Ns][0] = 1.0;   // no WGD 
    }

    if(nstate_chr > 0){
        L_sk_k.lnl_table_chr.resize(2 * Ns + 1, vector<double>(nstate_chr, 0.0));
        L_sk_k.lnl_table_chr[Ns][abs(obs_decomp.max_chr_change)] = 1.0;  // no chr gain/loss
    }
    if(nstate_site > 0){
        L_sk_k.lnl_table_seg.resize(2 * Ns + 1, vector<double>(nstate_site, 0.0));
        L_sk_k.lnl_table_seg[Ns][abs(obs_decomp.max_site_change)] = 1.0;  // no site gain/loss  
    }

    for(int i = 0; i < Ns; ++i){
        // For WGD
        // WGD state is 0,1,2
        if(nstate_wgd > 0){
            int state_wgd = obs_change[i].num_wgd;
            L_sk_k.lnl_table_wgd[i][state_wgd] = 1.0; 
        }
       
        // For chromosome gain/loss       
        // chr state is likely -1, 0, 1
        if(nstate_chr > 0){
            int state_chr = obs_change[i].cn_change_chr + abs(obs_decomp.max_chr_change); // change value to positive index
            // compute index of the observed state        
            L_sk_k.lnl_table_chr[i][state_chr] = 1.0; 
        }
        
        // For site gain/loss
        // site state is likely -2, -1, 0, 1, 2
        if(nstate_site > 0){
            int state_site = obs_change[i].cn_change_site + abs(obs_decomp.max_site_change); // change value to positive index
            L_sk_k.lnl_table_seg[i][state_site] = 1.0; 
        }
    }   

    if(debug){
        cout << "Likelihood table at tips for WGD:";
        print_lnl_at_tips<CN_CHANGE>(rtree, obs_change, L_sk_k.lnl_table_wgd, nstate_wgd);
        cout << "Likelihood table at tips for chromosome gain/loss:";
        print_lnl_at_tips<CN_CHANGE>(rtree, obs_change, L_sk_k.lnl_table_chr, nstate_chr);
        cout << "Likelihood table at tips for site duplication/deletion:";
        print_lnl_at_tips<CN_CHANGE>(rtree, obs_change, L_sk_k.lnl_table_seg, nstate_site);   
    }
}


/**
 * @brief Build rate matrices for different levels of copy number changes.
 * 
 * @param rtree: the evolutionary tree 
 * @param obs_decomp: observed changes for different levels, including maximum changes 
 * @param dim_decomp: dimensions of different levels    
 * @return QMAT_DECOMP: rate matrices for different levels 
 */
void build_rate_matrices(QMAT_DECOMP& qmat_decomp, const evo_tree& rtree, const OBS_DECOMP& obs_decomp, const DIM_DECOMP& dim_decomp, int debug){
    int max_wgd = obs_decomp.max_wgd;
    int max_chr_change = obs_decomp.max_chr_change;
    int max_site_change = obs_decomp.max_site_change;

    int dim_wgd = dim_decomp.dim_wgd;
    int dim_chr = dim_decomp.dim_chr;
    int dim_seg = dim_decomp.dim_seg;

    int dim_mat_wgd = dim_wgd * dim_wgd;
    int dim_mat_chr = dim_chr * dim_chr;
    int dim_mat_seg = dim_seg * dim_seg;

    
    double *qmat_wgd = nullptr;
    double *qmat_chr = nullptr;
    double *qmat_seg = nullptr;

    if(max_wgd > 0){
        qmat_wgd = new double[dim_mat_wgd];  // WGD
        memset(qmat_wgd, 0.0, dim_mat_wgd * sizeof(double));
        get_rate_matrix_wgd(qmat_wgd, rtree.wgd_rate, max_wgd);
    }

    if(max_chr_change > 0){
        qmat_chr = new double[dim_mat_chr];   // chromosome gain/loss
        memset(qmat_chr, 0.0, dim_mat_chr * sizeof(double));
        get_rate_matrix_chr_change(qmat_chr, rtree.chr_gain_rate, rtree.chr_loss_rate, max_chr_change);
    }

    if(max_site_change > 0){
        qmat_seg = new double[dim_mat_seg];  // site duplication/deletion
        memset(qmat_seg, 0.0, dim_mat_seg * sizeof(double));
        get_rate_matrix_site_change(qmat_seg, rtree.dup_rate, rtree.del_rate, max_site_change);
    }

    qmat_decomp.free_all();
    qmat_decomp.qmat_wgd = qmat_wgd;
    qmat_decomp.qmat_chr = qmat_chr;
    qmat_decomp.qmat_seg = qmat_seg;

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
 * @param qmat: rate matrices for different levels 
 * @param dim: dimensions of different levels 
 * @param max_wgd: maximum number of WGD events 
 * @param max_chr_change: maximum number of chromosome gain/loss events 
 * @param max_site_change: maximum number of site duplication/deletion events 

 */
void build_transition_matrices(PMAT_DECOMP& pmat, const evo_tree& rtree, const vector<int>& knodes, const QMAT_DECOMP& qmat_decomp, const DIM_DECOMP& dim_decomp, const OBS_DECOMP& obs_decomp, int debug){
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
         if(obs_decomp.max_wgd > 0){
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
         if(obs_decomp.max_chr_change > 0){
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
         if(obs_decomp.max_site_change > 0){
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

    pmat.free_all();  // or equivalent manual loop
    pmat.pmats_wgd = std::move(pmats_wgd);
    pmat.pmats_chr = std::move(pmats_chr);
    pmat.pmats_seg = std::move(pmats_seg);

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
}



/** 
 * @brief Compute the likelihood at a child node given the state at the parent node in a decomposed model 
 * 
 * @param node: index of the child node
 * @param s_wgd: starting index for WGD states at the parent node
 * @param s_chr: starting index for chromosome gain/loss states at the parent node
 * @param s_seg: starting index for site duplication/deletion states at the parent node
 * @param L_sk_k: likelihood table for each node and each possible combination of states
 * @param prob_decomp: transition probabilities for different levels  
 * @param dim_decomp: dimensions of different levels  
 * @param debug: debug flag
 * @return LNL_VAL: likelihood values for different levels 
 */
LNL_VAL compute_child_likelihood(int node, const CN_CHANGE& sk, const LNL_TABLE& L_sk_k, const PROB_DECOMP1& prob_decomp, const DIM_DECOMP& dim_decomp, int debug) {
    double L_wgd = 0.0;
    double L_chr = 0.0;
    double L_seg = 0.0;

    for (int e = 0; e < dim_decomp.dim_wgd; ++e) {
        double p = prob_decomp.pbl_wgd[sk.num_wgd + e * dim_decomp.dim_wgd];
        L_wgd += p * L_sk_k.lnl_table_wgd[node][e];
    }

    for (int e = 0; e < dim_decomp.dim_chr; ++e) {
        double p = prob_decomp.pbl_chr[sk.cn_change_chr + e * dim_decomp.dim_chr];
        L_chr += p * L_sk_k.lnl_table_chr[node][e];
    }

    for (int e = 0; e < dim_decomp.dim_seg; ++e) {
        double p = prob_decomp.pbl_seg[sk.cn_change_site + e * dim_decomp.dim_seg];
        L_seg += p * L_sk_k.lnl_table_seg[node][e];
    }

    if (debug) {
        std::cout << "\tLikelihood scoring for node " << node << " given parent state ("
                  << sk.num_wgd << ", "
                  << sk.cn_change_chr << ", "
                  << sk.cn_change_site << "): " 
                  << L_wgd << "\t"
                  << L_chr << "\t"
                  << L_seg << std::endl;
    }

    LNL_VAL result = {L_wgd, L_chr, L_seg};

    return result;
}



/** 
 * @brief Compute the probability of children nodes given the state at the parent node in a decomposed model 
 * @param L_sk_k: likelihood table for each node and each possible combination of states
 * @param rtree: evolutionary tree  
 * @param decomp_table: decomposition table for different copy numbers
 * @param sk: state at the parent node
 * @param cn_max: maximum copy number allowed
 * @param nstate: number of possible states             
 * @param prob_decomp: transition probabilities for different levels  
 * @param dim_decomp: dimensions of different levels  
 * @param ni: index of the first child node
 * @param nj: index of the second child node
 * @param bli: branch length from node sk to the first child node ni
 * @param blj: branch length from node sk to the second child node nj
 * @param is_total: whether the observed data is total copy number      
*/
LNL_VAL get_prob_children_change(LNL_TABLE& L_sk_k, const evo_tree& rtree, const CN_CHANGE& sk, PROB_DECOMP1& prob_decompi, PROB_DECOMP1& prob_decompj, const DIM_DECOMP& dim_decomp, int ni, int nj, int is_total, int debug){
    if(debug) cout << "\tget_prob_children_change" << endl;
    
    LNL_VAL li_vals = compute_child_likelihood(ni, sk, L_sk_k, prob_decompi, dim_decomp, debug);
    double Li = li_vals.lnl_wgd * li_vals.lnl_chr * li_vals.lnl_seg;
    if(debug) cout << "\tscoring: Li\t" << li_vals.lnl_wgd << "\t" << li_vals.lnl_chr << "\t" << li_vals.lnl_seg << endl;
    
    LNL_VAL lj_vals = compute_child_likelihood(nj, sk, L_sk_k, prob_decompj, dim_decomp, debug);
    double Lj = lj_vals.lnl_wgd * lj_vals.lnl_chr * lj_vals.lnl_seg;
    if(debug) cout << "\tscoring: Lj\t" << lj_vals.lnl_wgd << "\t" << lj_vals.lnl_chr << "\t" << lj_vals.lnl_seg << endl;

    LNL_VAL result = {li_vals.lnl_wgd * lj_vals.lnl_wgd, li_vals.lnl_chr * lj_vals.lnl_chr, li_vals.lnl_seg * lj_vals.lnl_seg};

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
 * @param is_total: whether to consider total copy number only, to use later
 */
// Sum over all possible states for initial and final nodes
void get_likelihood_site_change(LNL_TABLE& L_sk_k, const evo_tree& rtree, const vector<int>& knodes, const PMAT_DECOMP& pmat_decomp, const DIM_DECOMP& dim_decomp, const OBS_DECOMP& obs_decomp, int is_total, int debug){
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
        prob_decompi.pbl_wgd = pbli_wgd;
        prob_decompi.pbl_chr = pbli_chr;
        prob_decompi.pbl_seg = pbli_seg;

        PROB_DECOMP1 prob_decompj;
        prob_decompj.pbl_wgd = pblj_wgd;       
        prob_decompj.pbl_chr = pblj_chr;        
        prob_decompj.pbl_seg = pblj_seg;

        if(debug) cout << "node: " << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , " << bli << "\t" <<  nj + 1 << " , " << blj << endl;

        // loop over possible states of internal nodes
        if(k == rtree.nleaf){    // root node is always normal
            CN_CHANGE sk = {0, 0, 0};   // normal state in terms of decomposed copy number changes e.g. (0,1,2, -1,0,1, -1,0,1)
            // int sk_idx = change_to_state(sk, cn_max);
            LNL_VAL lnl_val = get_prob_children_change(L_sk_k, rtree, sk, prob_decompi, prob_decompj, dim_decomp, ni, nj, is_total, debug);
            if(dim_decomp.dim_wgd > 0) L_sk_k.lnl_table_wgd[k][0] = lnl_val.lnl_wgd;
            if(dim_decomp.dim_chr > 0) L_sk_k.lnl_table_chr[k][0] = lnl_val.lnl_chr;
            if(dim_decomp.dim_seg > 0) L_sk_k.lnl_table_seg[k][0] = lnl_val.lnl_seg;
            if(debug) cout << "Getting likelihood for root node " << k << endl;
        }else{
            for(int i = 0; i < obs_decomp.max_wgd; ++i){               
                 CN_CHANGE sk = {i, 0, 0}; 
                 LNL_VAL lnl_val = get_prob_children_change(L_sk_k, rtree, sk, prob_decompi, prob_decompj, dim_decomp, ni, nj, is_total, debug);
                 L_sk_k.lnl_table_wgd[k][i] = lnl_val.lnl_wgd;
                 if(dim_decomp.dim_chr > 0) L_sk_k.lnl_table_chr[k][0] = lnl_val.lnl_chr;
                 if(dim_decomp.dim_seg > 0) L_sk_k.lnl_table_seg[k][0] = lnl_val.lnl_seg;
            }
            for(int j = -obs_decomp.max_chr_change; j < obs_decomp.max_chr_change; ++j){
                 CN_CHANGE sk = {0, j, 0}; 
                 LNL_VAL lnl_val = get_prob_children_change(L_sk_k, rtree, sk, prob_decompi, prob_decompj, dim_decomp, ni, nj, is_total, debug);
                 if(dim_decomp.dim_wgd > 0) L_sk_k.lnl_table_wgd[k][0] = lnl_val.lnl_wgd;
                 L_sk_k.lnl_table_chr[k][j] = lnl_val.lnl_chr;
                 if(dim_decomp.dim_seg > 0) L_sk_k.lnl_table_seg[k][0] = lnl_val.lnl_seg;
            }
            for(int l = -obs_decomp.max_site_change; l < obs_decomp.max_site_change; ++l){
                 CN_CHANGE sk = {0, 0, l}; 
                 LNL_VAL lnl_val = get_prob_children_change(L_sk_k, rtree, sk, prob_decompi, prob_decompj, dim_decomp, ni, nj, is_total, debug);
                 if(dim_decomp.dim_wgd > 0) L_sk_k.lnl_table_wgd[k][0] = lnl_val.lnl_wgd;
                 if(dim_decomp.dim_chr > 0) L_sk_k.lnl_table_chr[k][0] = lnl_val.lnl_chr;
                 L_sk_k.lnl_table_seg[k][l] = lnl_val.lnl_seg;  
            }
            if(debug) cout << "Getting likelihood for internal node " << k << endl;
        }
  }
  
  if(debug > 0){
    print_tree_lnl(rtree, L_sk_k.lnl_table_wgd, dim_decomp.dim_wgd);
    print_tree_lnl(rtree, L_sk_k.lnl_table_chr, dim_decomp.dim_chr);
    print_tree_lnl(rtree, L_sk_k.lnl_table_seg, dim_decomp.dim_seg);
  }
}


/** 
 *  @brief Extract the log likelihood from the likelihood table at the root for decomposed model
 *
 *  @param L_sk_k Likelihood tables for WGD, chromosome gain/loss, and site duplication/deletion
 *  @param Ns Number of samples
 */
// Get the likelihood of the tree from likelihood table of state combinations
double extract_tree_lnl_change(const LNL_TABLE& L_sk_k, int Ns, int debug){
    if(debug) cout << "Extracting likelihood for the root" << endl;

    double log_lnl = 0.0;

    // WGD
    if (!L_sk_k.lnl_table_wgd.empty()) {
        double prob_wgd = L_sk_k.lnl_table_wgd[Ns + 1][0];   
        if(prob_wgd > 0.0) {
            log_lnl += log(prob_wgd);
        } else {
            if (debug) cout << "Zero likelihood at root: WGD" << endl;
        }
    } else {
        if (debug) cout << "No WGD events in the tree" << endl;
    }

    // Chromosome
    if(!L_sk_k.lnl_table_chr.empty()){
        double prob_chr = L_sk_k.lnl_table_chr[Ns + 1][0];  
        if (prob_chr > 0.0) {
            log_lnl += log(prob_chr);
        } else {
            if (debug) cout << "Zero likelihood at root: chromosome change" << endl;
        }
    } else {
        if (debug) cout << "No chromosome gain/loss events in the tree" << endl;
    }

    // Segment
    if(!L_sk_k.lnl_table_seg.empty()){
        double prob_seg = L_sk_k.lnl_table_seg[Ns + 1][0]; 
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
 * 
 * @param vobs: observed copy number changes at different chromosomes
 * @param obs_decomp: decomposition of observed copy numbers
 * @param rtree: the evolutionary tree 
 * @param comps: set of possible decomposed component combinations
 * @param knodes: list of internal nodes in post-order traversal
 * @param pmat_decomp: transition probability matrices for different levels
 * @param dim_decomp: dimensions of different levels
 * @param lnl_type: tags when computing likelihood, such as indicator of whether to use repeated sites to speed up computation
 * @param cn_max: maximum copy number allowed
 * @param is_total: whether to consider total copy number only
 * @return double: log likelihood value
 */
double get_likelihood_chr_change(const evo_tree& rtree, const map<int, vector<vector<CN_CHANGE>>>& vobs, const vector<int>& knodes, const PMAT_DECOMP& pmat_decomp, const DIM_DECOMP& dim_decomp, const OBS_DECOMP& obs_decomp, const LNL_TYPE& lnl_type, int is_total, int debug){
    double logL = 0.0;    // for all chromosomes

    // Use a map to store computed log likelihood
    map<vector<CN_CHANGE>, LNL_TABLE> sites_lnl_map;

    // for each chromosome
    for(auto vcn : vobs){
      int nchr = vcn.first;
      vector<vector<CN_CHANGE>> obs_chr = vcn.second;
      double site_logL = 0.0;   // log likelihood for all sites on a chromosome

      if(debug){
        cout << "Computing likelihood on Chr " << nchr <<  " with " << obs_chr.size() << " sites" << endl;
      }
     
      for(size_t nc = 0; nc < obs_chr.size(); nc++){    // for each site on the chromosome (may be repeated)
          vector<CN_CHANGE> obs = obs_chr[nc];
          LNL_TABLE L_sk_k;

          if(debug){
              cout << "Site " << nc + 1 << " observed changes: ";
              for(int i = 0; i < obs.size(); i++){
                  cout << "\t" << obs[i];
              }
              cout << endl;
          }
         
          if(lnl_type.use_repeat){ 
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                if(debug) cout << "\tsites new" << endl;
                  initialize_lnl_table_change(L_sk_k, rtree, obs, dim_decomp, obs_decomp, debug);
                  get_likelihood_site_change(L_sk_k, rtree, knodes, pmat_decomp, dim_decomp, obs_decomp, is_total, debug);
                  sites_lnl_map[obs] = L_sk_k;
              }else{
                  if(debug) cout << "\tsites repeated" << endl;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
                if(debug) cout << "\tsites no repeat consideration" << endl;
                initialize_lnl_table_change(L_sk_k, rtree, obs, dim_decomp, obs_decomp, debug);
                get_likelihood_site_change(L_sk_k, rtree, knodes, pmat_decomp, dim_decomp, obs_decomp, is_total, debug);
          }

          site_logL += extract_tree_lnl_change(L_sk_k, rtree.nleaf - 1, debug);
          if(debug){
              // cout << "Crtree.nleaf - 1 at this site: ";
              for(int i = 0; i < obs.size(); i++){
                  cout << "\t" << obs[i];
              }
              cout << endl;
              print_tree_lnl(rtree, L_sk_k.lnl_table_wgd, dim_decomp.dim_wgd);
              print_tree_lnl(rtree, L_sk_k.lnl_table_chr, dim_decomp.dim_chr);
              print_tree_lnl(rtree, L_sk_k.lnl_table_seg, dim_decomp.dim_seg);  
          }
      }

      logL += site_logL;

      if(debug){
          cout << "\nLikelihood for chromosome " << nchr << " is " << site_logL << endl;
      }
    } // for each chromosome

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
    int cn_max = lnl_type.cn_max;
    int is_total = lnl_type.is_total;

    int max_wgd = obs_decomp.max_wgd;
    // TODO: may consider gain/loss separately and arm-level changes in the future
    int max_chr_change = obs_decomp.max_chr_change;
    int max_site_change = obs_decomp.max_site_change;

    //   // use obs_decomp to avoid new parameters
    //   MAX_CHANGE max_change = {max_wgd, max_chr_change, max_site_change};

    int dim_wgd = 0;
    int dim_chr = 0;
    int dim_seg = 0;

    // only consider states that can reach observed changes
    if(max_wgd > 0) dim_wgd = max_wgd + 1;
    if(max_chr_change > 0) dim_chr = 2 * max_chr_change + 1;
    if(max_site_change > 0) dim_seg = 2 * max_site_change + 1;

    int dim_mat_wgd = dim_wgd * dim_wgd;
    int dim_mat_chr = dim_chr * dim_chr;
    int dim_mat_seg = dim_seg * dim_seg;  

    DIM_DECOMP dim_decomp;
    dim_decomp.dim_wgd = dim_wgd;
    dim_decomp.dim_chr = dim_chr;
    dim_decomp.dim_seg = dim_seg;

    cout << "Dimensions of rate/transition matrices for WGD, chr gain/loss, site gain/loss: " << dim_wgd << "\t" << dim_chr << "\t" << dim_seg << endl;

    if(debug) cout << "\tBuilding Q rate matrices for multiple levels" << endl;
    QMAT_DECOMP qmat_decomp;
    build_rate_matrices(qmat_decomp, rtree, obs_decomp, dim_decomp, debug);


    if(debug) cout << "\tBuilding P transition matrices for multiple levels" << endl;
    PMAT_DECOMP pmat_decomp;
    build_transition_matrices(pmat_decomp, rtree, lnl_type.knodes, qmat_decomp, dim_decomp, obs_decomp, debug);

    double logL = get_likelihood_chr_change(rtree, vobs_change, lnl_type.knodes, pmat_decomp, dim_decomp, obs_decomp, lnl_type, is_total, debug);

    if(debug) cout << "Final likelihood before correcting acquisition bias: " << logL << endl;
    if(lnl_type.correct_bias){
        if(debug) cout << "Correcting for the skip of invariant sites" << endl;
        vector<CN_CHANGE> obs(rtree.nleaf - 1, {0, 0, 0});   // invariant site
        LNL_TABLE L_sk_k;
        initialize_lnl_table_change(L_sk_k, rtree, obs, dim_decomp, obs_decomp, debug);   
        get_likelihood_site_change(L_sk_k, rtree, lnl_type.knodes, pmat_decomp, dim_decomp, obs_decomp, is_total, debug);
        double lnl_invar = extract_tree_lnl_change(L_sk_k, rtree.nleaf - 1, debug);

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
    double dt_valid= elapsed_seconds(start_time, end_all_valid);
    
    time_decomp_all   += dt_valid; // valid also count toward all, so time for all is always updated
    time_decomp_valid += dt_valid; // only update time for valid calls
    ++cnt_decomp_valid; // only update count for valid calls

    if(debug > 1)  cout << "Likelihood calculation time for this tree: " << dt_valid << " seconds." << endl;

    return logL;
}



/************* Funtions when decomposing all possible copy number states **********************/
/**
 * @brief Initialize the likelihood table for decomposed copy number states.
 * 
 * @param obs: observed copy number changes at the tips 
 * @param obs_decomp: decomposition information of observed copy numbers
 * @param chr: chromosome number starting from 1
 * @param rtree: the evolutionary tree 
 * @param comps: set of possible decomposed component combinations
 * @param infer_wgd: whether to consider WGD inference
 * @param infer_chr: whether to consider chromosome gain/loss inference
 * @param cn_max: maximum copy number allowed
 * @param is_total: whether to consider total copy number only
 * @return vector<vector<double>> 
 */
// L_sk_k has one row for each tree node and one column for each possible state; chr starting from 1
// This function is critical in obtaining correct likelihood. 
// If one tip is not initialized, the final likelihood will be 0.
vector<vector<double>> initialize_lnl_table_decomp(vector<int>& obs, const OBS_DECOMP& obs_decomp, int chr, const evo_tree& rtree, const set<vector<int>>& comps, int infer_wgd, int infer_chr, int cn_max, int is_total){
    int debug = 0;
    // construct a table for each state of each node
    int ntotn = 2 * rtree.nleaf - 1;
    // each state corresponds to one combination of components for a copy number
    // a copy number can be decomposed into multiple combinations of components when no constraint is applied
    int nstate = comps.size();      // number of possible states, depending on decomposition
    vector<vector<double>> L_sk_k(ntotn, vector<double>(nstate, 0.0));

    for(int i = 0; i < rtree.nleaf - 1; ++i){
        // For total copy number, all the possible combinations have to be considered.
        // Set related allele specific cases to be 1, with index from obs[i] * (obs[i] + 1)/2 to obs[i] * (obs[i] + 1)/2 + obs[i]. The index is computed based on pre-specified order.
        int cn = obs[i];
        int num_change = 0;  // record chromosome CN changes
        if(infer_chr && chr > 0){
            num_change = obs_decomp.sample_change_chr[i][chr-1];
        }
        if(!is_total){    // changing input haplotype-specific copy number to total
            cn = state_to_total_cn(obs[i], cn_max);
        }
        if(debug) cout << "\nStart filling likelihood table for sample "  << i + 1 << " chromosome " << chr  << " copy number " << cn << endl;
        // Fill all the possible state combinations
        int k = 0;
        for (auto c : comps){
            if(debug){
                cout << "\tcn vector:";
                for(int k = 0; k < c.size(); k++){
                    cout << "\t" << c[k];
                }
                cout << endl;
            }

            int alpha = c[0];  // wgd component
            // If a sample has one WGD event, the corresponding component must be the same
            if(infer_wgd && alpha != obs_decomp.sample_num_wgd[i]) continue;  // skip incompatible WGD cases to save time

            // WGD may occur before, at least one chromosome gain before or after WGD
            if(num_change >= 1 && (c[1] <= 0 && c[3] <= 0)){
                if(debug){
                    cout << "\t\tpotential chromosome gain in sample "  << i + 1 << " chromosome " << chr << " is " << num_change << endl;
                }
                continue;
            }
            if(num_change <= -1 && (c[1] >= 0 && c[3] >= 0)){
                if(debug){
                    cout << "\t\tpotential chromosome loss in sample "  << i + 1 << " chromosome " << chr << " is " << num_change << endl;
                }
                continue;
            }

            int sum = pow(2, alpha + 1) + c[5] * c[1] + c[2] + 2 * c[6] * c[3] + 2 * c[4];
            if(sum == cn){
                if(debug) cout << "\t\tfilling 1 when sum == " << cn << endl;
                L_sk_k[i][k] = 1.0;
            }

            // // assuming m_max >= 1. It is likely that all copies of a segment is lost before chromosome gain/loss
            // for(int m1 = 0; m1 <= obs_decomp.m_max; m1++){
            //     for(int m2 = 0; m2 <= obs_decomp.m_max; m2++){
            //         int sum = pow(2, alpha + 1) + m1 * c[1] + c[2] + 2 * m2 * c[3] + 2 * c[4];
            //         if(sum == cn){
            //             if(debug) cout << "\t\tfilling 1 when sum == " << cn << endl;
            //             L_sk_k[i][k] = 1.0;
            //         }
            //     }
            // }

            k++;
        }
        // each row should have one entry being 1
        int sum = 0;
        for(int j = 0; j < L_sk_k[i].size(); j++){
            sum += L_sk_k[i][j];
        }
        if(sum < 1){
            cout << "Error in filling table for copy number " << cn << " in sample " << i + 1 << " chromosome " << chr << endl << endl;
        }
    }
    // set likelihood for normal sample
    // Fill all the possible state combinations
    for(int j = 0; j < comps.size(); j++){
        int k = 0;
        for(auto v : comps){
            bool zeros = all_of(v.begin(), v.end(), [](int i) { return i == 0; });
            if(zeros){
                L_sk_k[rtree.nleaf - 1][k] = 1.0;
                break;
            }
            k++;
        }
    }

    if(debug){
        print_lnl_at_tips<int>(rtree, obs, L_sk_k, nstate);
    }

    return L_sk_k;
}


// Assume the likelihood table is for each combination of states
double get_prob_children_decomp(vector<vector<double>>& L_sk_k, const evo_tree& rtree, map<int, set<vector<int>>>& decomp_table, int sk, int cn_max, int nstate, PROB_DECOMP& prob_decomp, DIM_DECOMP& dim_decomp, int ni, int nj, int bli, int blj, int is_total){
    int debug = 0;
    int s_wgd, s_chr, s_seg;
    int e_wgd, e_chr, e_seg;
    // Each copy number is decomposed into a set of tuples
    set<vector<int>> comp_start = decomp_table[sk];


    double Li = 0.0;
    for(auto s : comp_start){
        s_wgd = s[0];
        s_chr = s[1] + s[3];    // add before and after WGD changes
        s_seg = s[2] + s[4];
        // loop over possible si
        for(int si = 0; si < nstate; ++si){
            // get the start and end state for each type
            int cn = state_to_total_cn(si, cn_max);
            set<vector<int>> comp_end = decomp_table[cn];

            for(auto e : comp_end){
                e_wgd = e[0];
                e_chr = e[1] + e[3];
                e_seg = e[2] + e[4];
                double prob_wgd = prob_decomp.pbli_wgd[s_wgd + e_wgd * dim_decomp.dim_wgd];
                double prob_chr = prob_decomp.pbli_chr[s_chr + e_chr * dim_decomp.dim_chr];
                double prob_seg = prob_decomp.pbli_seg[s_seg + e_seg * dim_decomp.dim_seg];
                double prob = prob_wgd * prob_chr * prob_seg * L_sk_k[ni][si];
                Li += prob;
                if(debug) cout << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob << "\n";
            }
          if(debug) cout << "\tscoring: Li\t" << Li << endl;
        }
    }

    double Lj = 0.0;
    for(auto s : comp_start){
        s_wgd = s[0];
        s_chr = s[1] + s[3];
        s_seg = s[2] + s[4];
        // loop over possible sj
        for(int sj = 0; sj < nstate; ++sj){
            int cn = state_to_total_cn(sj, cn_max);
            set<vector<int>> comp_end = decomp_table[cn];

            for(auto e : comp_end){
                e_wgd = e[0];
                e_chr = e[1] + e[3];
                e_seg = e[2] + e[4];
                double prob_wgd = prob_decomp.pblj_wgd[s_wgd + e_wgd * dim_decomp.dim_wgd];
                double prob_chr = prob_decomp.pblj_chr[s_chr + e_chr * dim_decomp.dim_chr];
                double prob_seg = prob_decomp.pblj_seg[s_seg + e_seg * dim_decomp.dim_seg];
                double prob = prob_wgd * prob_chr * prob_seg * L_sk_k[nj][sj];
                Lj += prob;
                if(debug) cout << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob << "\n";
            }
        }
        if(debug) cout << "\tscoring: Lj\t" << Lj << endl;
    }

    return Li * Lj;
}


// Assume the likelihood table is for each combination of states
// distinguish changes before and after WGD
double get_prob_children_decomp2(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const set<vector<int>>& comps, int sk, int cn_max, int nstate, PROB_DECOMP& prob_decomp, DIM_DECOMP& dim_decomp, int ni, int nj, int bli, int blj, int is_total){
    int debug = 0;
    int s_wgd, s_chr, s_seg, s_chr2, s_seg2;
    int e_wgd, e_chr, e_seg, e_chr2, e_seg2;
    double prob_wgd, prob_chr, prob_seg, prob_chr2, prob_seg2, prob_chr_all, prob_seg_all;
    // Each copy number is decomposed into a set of tuples relative to WGD
    set<vector<int>>::iterator iter = comps.begin();
    // It will move forward the passed iterator by passed value
    advance(iter, sk);
    vector<int> s = *iter;

    s_wgd = s[0];
    s_chr = s[1];
    s_seg = s[2];
    s_chr2 = s[3];
    s_seg2 = s[4];

    int dim_wgd = dim_decomp.dim_wgd;
    int dim_chr = dim_decomp.dim_chr;
    int dim_seg = dim_decomp.dim_seg;

    // The indices for chromosome and segment matrix have to be adjusted
    int delta_chr = (dim_chr - 1) / 2;
    int delta_seg = (dim_seg - 1) / 2;

    if(debug){
        cout << "Starting state " << sk << ":\t" << s_wgd << "\t" << s_chr << "\t" << s_seg << "\t" << s_chr2 << "\t" << s_seg2 << "\n";
        cout << "   Offset for chr and segment matrices " << delta_chr << "\t" << delta_seg << "\n";
    }


    double Li = 0.0;
    int si = 0;
    for(auto e : comps){
        double prob = 0;
        prob_wgd = 1;
        prob_chr = 1;
        prob_seg = 1;
        prob_chr2 = 1;
        prob_seg2 = 1;
        prob_chr_all = 1;
        prob_seg_all = 1;
        if(L_sk_k[ni][si] > 0){
            e_wgd = e[0];
            e_chr = e[1];
            e_seg = e[2];
            e_chr2 = e[3];
            e_seg2 = e[4];
            if(dim_wgd > 1) prob_wgd = prob_decomp.pbli_wgd[s_wgd + e_wgd * dim_wgd];
            if(dim_chr > 1){
                prob_chr = prob_decomp.pbli_chr[(s_chr + delta_chr) + (e_chr + delta_chr) * dim_chr];
                prob_chr2 = prob_decomp.pbli_chr[(s_chr2 + delta_chr) + (e_chr2 + delta_chr) * dim_chr];
                prob_chr_all = prob_chr * prob_chr2;
            }
            if(dim_seg > 1){
                prob_seg = prob_decomp.pbli_seg[(s_seg + delta_seg) + (e_seg + delta_seg) * dim_seg];
                prob_seg2 = prob_decomp.pbli_seg[(s_seg2 + delta_seg) + (e_seg2 + delta_seg) * dim_seg];
                prob_seg_all = prob_seg * prob_seg2;
            }
            prob = prob_wgd * prob_chr_all * prob_seg_all * L_sk_k[ni][si];
            if(debug){
                cout << "End state " << si << ":\t" << e_wgd << "\t" << e_chr << "\t" << e_seg << "\t" << e_chr2 << "\t" << e_seg2 << "\n";
                cout << "Probability for each type of event " << "\t" << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob_chr2 << "\t" << prob_seg2 << "\t" << prob << "\n";
            }
        }
        Li += prob;
        si++;
    }
    if(debug) cout << "\tscoring: Li\t" << Li << endl;

    double Lj = 0.0;
    int sj = 0;
    for(auto e : comps){
        double prob = 0;
        prob_wgd = 1;
        prob_chr = 1;
        prob_seg = 1;
        prob_chr2 = 1;
        prob_seg2 = 1;
        prob_chr_all = 1;
        prob_seg_all = 1;
        if(L_sk_k[nj][sj] > 0){
            e_wgd = e[0];
            e_chr = e[1];
            e_seg = e[2];
            e_chr2 = e[3];
            e_seg2 = e[4];
            if(dim_wgd > 1) prob_wgd = prob_decomp.pblj_wgd[s_wgd + e_wgd * dim_wgd];
            if(dim_chr > 1){
                prob_chr = prob_decomp.pblj_chr[(s_chr + delta_chr) + (e_chr + delta_chr) * dim_chr];
                prob_chr2 = prob_decomp.pblj_chr[(s_chr2 + delta_chr) + (e_chr2 + delta_chr) * dim_chr];
                prob_chr_all = prob_chr * prob_chr2;
            }
            if(dim_seg > 1){
                prob_seg = prob_decomp.pblj_seg[(s_seg + delta_seg) + (e_seg + delta_seg) * dim_seg];
                prob_seg2 = prob_decomp.pblj_seg[(s_seg2 + delta_seg) + (e_seg2 + delta_seg) * dim_seg];
                prob_seg_all = prob_seg * prob_seg2;
            }
            prob = prob_wgd * prob_chr_all * prob_seg_all * L_sk_k[nj][sj];
            if(debug){
                cout << "End state " << sj << ": \t" << e_wgd << "\t" << e_chr << "\t" << e_seg << "\t" << e_chr2 << "\t" << e_seg2 << "\n";
                cout << "Probability for each type of event " << "\t" << prob_wgd << "\t" << prob_chr << "\t" << prob_seg << "\t" << prob_chr2 << "\t" << prob_seg2 << "\t" << prob << "\n";
            }
        }
        Lj += prob;
        sj++;
    }
    if(debug) cout << "\tscoring: Lj\t" << Lj << endl;

    return Li * Lj;
}



/** 
 * @brief Get the likelihood on one site of a chromosome (assuming each observed copy number is composed of three type of events)   
 * 
 * @param L_sk_k: likelihood table for each node and each possible state
 * @param rtree: evolutionary tree
 * @param comps: set of possible decomposed component combinations
 * @param knodes: list of internal nodes in post-order traversal
 * @param pmat_decomp: transition probability matrices for different levels
 * @param dim_decomp: dimensions of different levels
 * @param cn_max: maximum copy number allowed
 * @param is_total: whether to consider total copy number only
 */
// Sum over all possible states for initial and final nodes
void get_likelihood_site_decomp(vector<vector<double>>& L_sk_k, const evo_tree& rtree, const set<vector<int>>& comps, const vector<int>& knodes, PMAT_DECOMP& pmat_decomp, DIM_DECOMP& dim_decomp, int cn_max, int is_total){
  int debug = 0;

  int dim_wgd = dim_decomp.dim_wgd;
  int dim_chr = dim_decomp.dim_chr;
  int dim_seg = dim_decomp.dim_seg;
  int nstate = comps.size();

  if(debug){
      cout << "Computing likelihood for one site" << endl;
      cout << dim_wgd << "\t" << dim_chr << "\t"  << dim_seg << "\t"  << nstate << endl;
  }

  for(int kn = 0; kn < knodes.size(); ++kn){
        int k = knodes[kn];
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

        double *pbli_wgd, *pblj_wgd;
        double *pbli_chr, *pblj_chr;
        double *pbli_seg, *pblj_seg;

        if(dim_wgd > 1){
            pbli_wgd = pmat_decomp.pmats_wgd[bli];
            pblj_wgd = pmat_decomp.pmats_wgd[blj];
        }
        if(dim_chr > 1){
            pbli_chr = pmat_decomp.pmats_chr[bli];
            pblj_chr = pmat_decomp.pmats_chr[blj];
        }
        if(dim_seg > 1){
            pbli_seg = pmat_decomp.pmats_seg[bli];
            pblj_seg = pmat_decomp.pmats_seg[blj];
        }

        PROB_DECOMP prob_decomp;
        prob_decomp.pbli_wgd = pbli_wgd;
        prob_decomp.pblj_wgd = pblj_wgd;
        prob_decomp.pbli_chr = pbli_chr;
        prob_decomp.pblj_chr = pblj_chr;
        prob_decomp.pbli_seg = pbli_seg;
        prob_decomp.pblj_seg = pblj_seg;

        if(debug) cout << "node: " << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , " << bli << "\t" <<  nj + 1 << " , " << blj << endl;

        // loop over possible observed states of start nodes
        if(k == rtree.nleaf){    // root node is always normal
            // int sk = 4;
            // L_sk_k[k][sk] = get_prob_children_decomp(L_sk_k, rtree, decomp_table, sk, cn_max, nstate, prob_decomp, dim_decomp, ni, nj, bli, blj, is_total);
            if(debug) cout << "Getting likelihood for root node " << k << endl;
            int sk = 0;
            for(auto v : comps){
                bool zeros = all_of(v.begin(), v.end(), [](int i) { return i == 0; });
                if(zeros){    // normal state
                    L_sk_k[k][sk] = get_prob_children_decomp2(L_sk_k, rtree, comps, sk, cn_max, nstate, prob_decomp, dim_decomp, ni, nj, bli, blj, is_total);
                    break;
                }
                sk++;
            }
        }else{
            for(int sk = 0; sk < nstate; ++sk){
                if(debug) cout << " Getting likelihood of children nodes " << endl;
                // L_sk_k[k][sk] = get_prob_children_decomp(L_sk_k, rtree, decomp_table, sk, cn_max, nstate, prob_decomp, dim_decomp, ni, nj, bli, blj, is_total);
                L_sk_k[k][sk] = get_prob_children_decomp2(L_sk_k, rtree, comps, sk, cn_max, nstate, prob_decomp, dim_decomp, ni, nj, bli, blj, is_total);
            }
        }
  }
  if(debug > 0){
    print_tree_lnl(rtree, L_sk_k, nstate);
  }
}



// Compute likelihood by chromosome, assuming each observed copy number is decomposed into three types of events
double get_likelihood_chr_decomp(const map<int, vector<vector<int>>>& vobs, const OBS_DECOMP& obs_decomp, const evo_tree& rtree, const set<vector<int>>& comps, const vector<int>& knodes, PMAT_DECOMP& pmat_decomp, DIM_DECOMP& dim_decomp, int infer_wgd, int infer_chr, int use_repeat, int cn_max, int is_total){
    int debug = 0;
    double logL = 0.0;    // for all chromosmes
    int nstate = comps.size();
    // Use a map to store computed log likelihood
    map<vector<int>, vector<vector<double>>> sites_lnl_map;

    // for each chromosome
    for(auto vcn : vobs){
      int nchr = vcn.first;
      if(debug){
        cout << "Computing likelihood on Chr " << nchr <<  " with " << vobs.at(nchr).size() << "sites" << endl;
      }
      double site_logL = 0.0;   // log likelihood for all sites on a chromosome

      // cout << " chromosome number change is " << 0 << endl;
      for(int nc = 0; nc < vobs.at(nchr).size(); nc++){    // for each segment on the chromosome
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs.at(nchr).at(nc);
          vector<vector<double>> L_sk_k;
          if(use_repeat){ 
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  L_sk_k = initialize_lnl_table_decomp(obs, obs_decomp, nchr, rtree, comps, infer_wgd, infer_chr, cn_max, is_total);
                  get_likelihood_site_decomp(L_sk_k, rtree, comps, knodes, pmat_decomp, dim_decomp, cn_max, is_total);
                  sites_lnl_map[obs] = L_sk_k;
              }else{
                  if(debug) cout << "\tsites repeated" << endl;
                  L_sk_k = sites_lnl_map[obs];
              }
          }else{
              L_sk_k = initialize_lnl_table_decomp(obs, obs_decomp, nchr, rtree, comps, infer_wgd, infer_chr, cn_max, is_total);
              get_likelihood_site_decomp(L_sk_k, rtree, comps, knodes, pmat_decomp, dim_decomp, cn_max, is_total);
          }

          site_logL += extract_tree_lnl_decomp(L_sk_k, comps, rtree.nleaf - 1);

          if(debug){
              // cout << "Crtree.nleaf - 1 at this site: ";
              for(int i = 0; i < obs.size(); i++){
                  cout << "\t" << obs[i];
              }
              cout << endl;
              print_tree_lnl(rtree, L_sk_k, nstate);
          }
      }

      logL += site_logL;

      if(debug){
          cout << "\nLikelihood for chromosome " << nchr << " is " << site_logL << endl;
      }
    } // for each chromosome

    return logL;
}


// Computing likelihood when WGD and chr gain/loss are incorporated
// Assume likelihood is for haplotype-specific information
// Additional parameters are included for decomposition-based likelihood calculation: 
// obs_decomp: observed decomposition information
// comps: set of all possible component vectors
double get_likelihood_decomp(evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const OBS_DECOMP& obs_decomp, const set<vector<int>>& comps, LNL_TYPE& lnl_type, int debug){
  if(debug) cout << "\tget likelihood using multiple chains based on copy number decomposition" << endl;

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

  int max_wgd = obs_decomp.max_wgd;
  int max_chr_change = obs_decomp.max_chr_change;
  int max_site_change = obs_decomp.max_site_change;

  // For WGD model
  int dim_wgd = max_wgd + 1;
  int dim_chr = 2 * max_chr_change + 1;
  int dim_seg = 2 * max_site_change + 1;

  int dim_mat_wgd = dim_wgd * dim_wgd;
  int dim_mat_chr = dim_chr * dim_chr;
  int dim_mat_seg = dim_seg * dim_seg;

  double *qmat_wgd, *qmat_chr, *qmat_seg;
  double *pmat_wgd, *pmat_chr, *pmat_seg;
  // Find the transition probability matrix for each branch
  map<double, double*> pmats_wgd;
  map<double, double*> pmats_chr;
  map<double, double*> pmats_seg;
  double *pmati_wgd, *pmatj_wgd;
  double *pmati_chr, *pmatj_chr;
  double *pmati_seg, *pmatj_seg;

  if(max_wgd > 0){
      qmat_wgd = new double[dim_mat_wgd];  // WGD
      memset(qmat_wgd, 0.0, dim_mat_wgd * sizeof(double));
      get_rate_matrix_wgd(qmat_wgd, rtree.wgd_rate, max_wgd);
  }

  if(max_chr_change > 0){
      qmat_chr = new double[dim_mat_chr];   // chromosome gain/loss
      memset(qmat_chr, 0.0, dim_mat_chr * sizeof(double));
      get_rate_matrix_chr_change(qmat_chr, rtree.chr_gain_rate, rtree.chr_loss_rate, max_chr_change);
  }

  if(max_site_change > 0){
      qmat_seg = new double[dim_mat_seg];  // site duplication/deletion
      memset(qmat_seg, 0.0, dim_mat_seg * sizeof(double));
      get_rate_matrix_site_change(qmat_seg, rtree.dup_rate, rtree.del_rate, max_site_change);
  }

  if(debug){
    cout << "Dimension of Q matrix at genome, chromosome, and segment level: " << dim_wgd << "\t" << dim_chr << "\t" << dim_seg << "\n";
  }

  vector<int> knodes = lnl_type.knodes;
  for(int kn = 0; kn < knodes.size(); ++kn){
         int k = knodes[kn];
         double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
         double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

         // For WGD
         if(max_wgd > 0){
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
         if(max_chr_change > 0){
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
         if(max_site_change > 0){
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

  PMAT_DECOMP pmat_decomp;
  pmat_decomp.pmats_wgd = pmats_wgd;
  pmat_decomp.pmats_chr = pmats_chr;
  pmat_decomp.pmats_seg = pmats_seg;

  DIM_DECOMP dim_decomp;
  dim_decomp.dim_wgd = dim_wgd;
  dim_decomp.dim_chr = dim_chr;
  dim_decomp.dim_seg = dim_seg;

  // cout << "Number of states is " << nstate << endl;
  double logL = get_likelihood_chr_decomp(vobs, obs_decomp, rtree, comps, knodes, pmat_decomp, dim_decomp, lnl_type.infer_wgd, lnl_type.infer_chr, lnl_type.use_repeat, cn_max, is_total);

  if(debug) cout << "Final likelihood before correcting acquisition bias: " << logL << endl;
  if(lnl_type.correct_bias){
    if(debug) cout << "Correcting for the skip of invariant sites" << endl;
    vector<int> obs(rtree.nleaf - 1, NORM_PLOIDY);
    vector<vector<double>> L_sk_k = initialize_lnl_table_decomp(obs, obs_decomp, 0, rtree, comps, lnl_type.infer_wgd, lnl_type.infer_chr, cn_max, is_total);
    get_likelihood_site_decomp(L_sk_k, rtree, comps, knodes, pmat_decomp, dim_decomp, cn_max, is_total);
    double lnl_invar = extract_tree_lnl_decomp(L_sk_k, comps, rtree.nleaf - 1);

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

  if(max_wgd > 0){
      delete [] qmat_wgd;
      for(auto m : pmats_wgd){
          delete [] m.second;
      }
  }
  if(max_chr_change > 0){
      delete [] qmat_chr;
      for(auto m : pmats_chr){
          delete [] m.second;
      }
  }
  if(max_site_change > 0){
      delete [] qmat_seg;
      for(auto m : pmats_seg){
          delete [] m.second;
      }
  }

  TimePoint end_all_valid = now();
  double dt_valid= elapsed_seconds(start_time, end_all_valid);
  
  time_decomp_all   += dt_valid; // valid also count toward all, so time for all is always updated
  time_decomp_valid += dt_valid; // only update time for valid calls
  ++cnt_decomp_valid; // only update count for valid calls

  if(debug > 1)  cout << "Likelihood calculation time for this tree: " << dt_valid << " seconds." << endl;

  return logL;
}



// Get the likelihood of the tree from likelihood table of state combinations
double extract_tree_lnl_decomp(vector<vector<double>>& L_sk_k, const set<vector<int>>& comps, int Ns){
    int debug = 0;
    if(debug) cout << "Extracting likelihood for the root" << endl;

    double likelihood = 0;
    int k = 0;

    for (auto v : comps){
        // if(debug) cout << k << "\t" << v[0] << "\t" << v[1] << "\t" << v[2] << "\t" << v[3] << "\t" << v[4] << endl;
        bool zeros = all_of(v.begin(), v.end(), [](int i){ return i == 0; });
        if(zeros){
            likelihood = L_sk_k[Ns + 1][k];
            break;
        }
        k++;
    }

    if(debug){
        for(int j = 0; j < comps.size(); ++j){
          cout << "\t" << L_sk_k[Ns + 1][j];
        }
        cout << endl;
    }

    if(likelihood > 0) return log(likelihood);
    else return LARGE_LNL;
}

