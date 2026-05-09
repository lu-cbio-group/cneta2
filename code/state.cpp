#include "state.hpp"


void print_tree_state(const evo_tree& rtree, const vector<vector<int>>& S_sk_k, int nstate){
    cout << "\nStates so far:\n";

    int ntotn = 2 * rtree.nleaf - 1;
    for(int i = 0; i < ntotn; ++i){
        for(int j = 0; j < nstate; ++j){
          cout << "\t" << S_sk_k[i][j];
        }
        cout << endl;
    }
}


/*************** functions for bounded model on site duplication/deletion *****************/


/**
 * @brief Extract the most likely joint estimation of ancestral nodes based on the likelihood and state tables obtained from get_ancestral_states_site 
 * @param rtree     The evolutionary tree     
 * @param knodes    The list of internal node indices in the tree, ordered from tips to root
 * @param comps         The set of possible copy number compositions for the decomposition model, where each composition is represented as a vector of 5 integers (wgd, chr_gain, seg_gain, chr_loss, seg_loss)
 * @param L_sk_k    The likelihood table for each node and each possible state of its parent, obtained from get_ancestral_states_site                
 * @param S_sk_k        The state table for each node and each possible state of its parent, recording the most likely state of the node given the parent state, obtained from get_ancestral_states_site 
 * @param model     The model used for ancestral state reconstruction, MK, BOUNDA, or DECOMP
 * @param cn_max            The maximum copy number considered in the model, only used for BOUNDA model
 * @param is_total      Whether the observed copy number is total copy number (true) or haplotype-specific copy number (false), only used for BOUNDA model
 * @param asr_states        The output map to store the most likely state for each node, where the key is the node index and the value is the state index
 */
void extract_tree_ancestral_state(const evo_tree& rtree, const vector<int>& knodes, const vector<vector<int>>& S_sk_k, int model, map<int, int>& asr_states){
    int debug = 0;

    if(debug){
        cout << "Get most likely joint estimation of ancestral nodes" << endl;
    }

    int parent_state;
    if(model == BOUNDA){
        parent_state = NORM_ALLElE_STATE;
    }else{
        parent_state = NORM_PLOIDY;
    }
    asr_states[rtree.nleaf] = parent_state;   // for root, adaption of step 4 of Pupko algorithm

    // Step 5 of Pupko algorithm
    // Traverse the tree from root to tips, need to know the parent of each node
    for(int i = knodes.size() - 2; i >= 0; i--){  // starting from node ID for MRCA
        int nid = knodes[i];
        // Find the parent node of current node
        int parent = rtree.nodes[nid].parent;
        if(asr_states.find(parent) == asr_states.end()){
            cout << "Cannot find state for the parent of node " << nid + 1 << ", " << parent + 1 << endl;
            exit(EXIT_FAILURE);
        }
        parent_state = asr_states[parent];

        int state = S_sk_k[nid][parent_state];
        asr_states[nid] = state;

        if(debug){
            cout << "\t\tnode " << nid + 1 << " with state " << state << " and parent " << parent + 1 << " whose state is " << parent_state << endl;
        }
    }

}


void print_tip_states(const evo_tree& rtree, int nstate, const vector<vector<int>>& S_sk_k){
    cout << "\nState vector for tips (observed data):\n";
    for (int i = 0; i < rtree.nleaf; ++i){
        for (int j = 0; j < nstate; ++j)
        {
            cout << "\t" << S_sk_k[i][j];
        }
        cout << endl;
    }
}


/** 
 * @brief Initialize the tables of tip nodes for reconstructing joint ancestral states (Step 1 of Pupko algorithm)
 * @param obs observed copy number at tips
 * @param rtree the evolutionary tree
 * @param blens vector of branch lengths for which transition probability matrices are computed
 * @param pmat_per_blen vector of transition probability matrices for each branch length in blens
 * @param L_sk_k likelihood table for each node and each possible state of its parent
 * @param S_sk_k state table for each node and each possible state of its parent, recording the most likely state of the node given the parent state
 * @param model the model used for ancestral state reconstruction, MK or BOUNDA
 * @param nstate number of states in the model
 * @param is_total whether the observed copy number is total copy number (true) or haplotype-specific copy number (false), only used for BOUNDA model
 */
void initialize_asr_table(const vector<int>& obs, const evo_tree& rtree, const vector<double>& blens, const vector<double*>& pmat_per_blen, vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, int model, int nstate, int is_total){
    int debug = 0;
    if(debug){
        cout << "Initializing tables for reconstructing joint ancestral state" << endl;
        print_branch_lengths(blens);
    }

    for(int i = 0; i < rtree.nleaf - 1; ++i){
        // cout << "node " << i + 1  << ", observed CN " << obs[i] << endl;
        // Find the parent node
        // int parent = rtree.edges[rtree.nodes[i].e_in].start;
        double blen = rtree.edges[rtree.nodes[i].e_in].length;
        // cout << "parent " << parent + 1 << ", blen " << blen << endl;

        auto pi = equal_range(blens.begin(), blens.end(), blen);
        // assert(distance(pi.first, pi.second) == 1);
        int idx_blen = distance(blens.begin(), pi.first);

        if(debug){
            cout << "Pmatrix for branch length " << blen << " " << blens[idx_blen] << endl;
            r8mat_print(nstate, nstate, pmat_per_blen[idx_blen], "  P matrix:");
        }

        // Find the state(s) of current node
        vector<int> tip_states;
        set_tip_states(model, is_total, obs[i], tip_states);

        if(debug){
            cout << "There are " << tip_states.size() << " states for copy number " << obs[i] << endl;
        }

        for(int j = 0; j < nstate; ++j){  // For each possible parent state, find the most likely tip states
            vector<double> vec_li(nstate, SMALL_LNL);
            // another loop as there maybe multiple states for a specific total CN
            for(int m = 0; m < tip_states.size(); ++m){
                int k = tip_states[m];
                if(debug){
                    cout << "parent state " << j << ", child state " << k << endl;
                }
                double li = 0.0;
                if(model == MK){
                    li = get_transition_prob(rtree.mu, blen, j, k);
                }else{
                    li = pmat_per_blen[idx_blen][j  + k * nstate];  // assume parent has state j
                }
                // if(li > 0) li = log(li);
                // else li = SMALL_LNL;
                vec_li[k] = li;
            }

            int max_i = distance(vec_li.begin(), max_element(vec_li.begin(), vec_li.end()));
            double max_li = *max_element(vec_li.begin(), vec_li.end());
            assert(max_li == vec_li[max_i]);
            // cout << "for node: i " << i + 1 << ", parent state " << j << ", max state is " << max_i << " with probability " << exp(max_li) << endl;
            S_sk_k[i][j] = max_i;
            L_sk_k[i][j] = max_li;
        }
    }

    if(debug){
      print_lnl_at_tips<int>(rtree, obs, L_sk_k, nstate);
      print_tip_states(rtree, nstate, S_sk_k);
    }
}


void set_tip_states(int model, int is_total, int obs_val, vector<int>& tip_states){
    if (model == BOUNDA){
        if (is_total){ // With total copy number, there may be multiple states corresponding to the same total CN, need to consider all of them
            int si = (obs_val * (obs_val + 1)) / 2;
            int ei = si + obs_val;
            for (int k = si; k <= ei; k++)
            {
                tip_states.push_back(k);
            }
        }else
        { // With haplotype-specific copy number, only the specific site needs to be filled, size 1
            tip_states.push_back(obs_val);
        }
    }else{
        tip_states.push_back(obs_val);
    }
}


// Find the state of a node k that gives the maximum likelihood under the given parent state sp, and record the state in S_sk_k
// return the maximum likelihood for node k given parent state sp
double get_max_prob_children(const vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, double* pblen, int k, int nstate, int sp, int ni, int nj, double blen, int model, int debug){
    vector<double> vec_li;
    double li = 0.0;
    // loop over possible si for a fixed state of parent (sp)
    for(int si = 0; si < nstate; ++si){
        if(model == MK){
          li = get_transition_prob(rtree.mu, blen, sp, si);
        }else{
          li = pblen[sp + si * nstate];
        }

        // if(debug){
        //     cout << "\tfor state " << si << endl;
        //     cout << "\t\tseparate likelihood "  << li << "\t" << L_sk_k[ni][si] << "\t" << L_sk_k[nj][si] << endl;
        // }

        li = li * L_sk_k[ni][si] * L_sk_k[nj][si];

        // if(debug) cout << "\t\tthe likelihood of the best reconstruction of subtree at " << sp + 1 << " is: " << li << endl;

        vec_li.push_back(li);
    }

    // step 2a of Pupko algorithm: find the state that gives the maximum likelihood and record in S_sk_k
    int max_i = distance(vec_li.begin(), max_element(vec_li.begin(), vec_li.end()));
    double max_li = *max_element(vec_li.begin(), vec_li.end());
    assert(max_li == vec_li[max_i]);

    S_sk_k[k][sp] = max_i;

    if(debug){
        cout << "All likelihoods:";
        for(int i = 0; i < vec_li.size(); i++){
            cout << "\t" << vec_li[i];
        }
        cout << endl;
        cout << "For node: k " << k + 1 << ", parent state " << sp << ", max state is " << max_i << " with probability " << max_li << endl;
    }

    return max_li;
}




/**
 * @brief Get the ancestral states site for nonroot internal nodes (Step 2 of Pupko algorithm)
 * 
 * @param L_sk_k 
 * @param S_sk_k 
 * @param rtree 
 * @param knodes 
 * @param blens 
 * @param pmat_per_blen 
 * @param nstate 
 * @param model 
 */
void get_ancestral_states_site(vector<vector<double>>& L_sk_k, vector<vector<int>>& S_sk_k, const evo_tree& rtree, const vector<int>& knodes, const vector<double>& blens, const vector<double*>& pmat_per_blen, int nstate, int model){
  int debug = 0;
  if(debug){
      cout << "Getting ancestral state for one site" << endl;
  }

  int Ns = rtree.nleaf - 1;
  for(int kn = 0; kn < knodes.size() - 1; ++kn){   // not include root which is always normal
    int k = knodes[kn];
    int np = rtree.edges[rtree.nodes[k].e_in].start;
    double blen = rtree.edges[rtree.nodes[k].e_in].length;
    int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
    int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;

    auto pi = equal_range(blens.begin(), blens.end(), blen);
    // assert(distance(pi.first, pi.second) == 1);
    int idx_blen = distance(blens.begin(), pi.first);
    double* pblen = pmat_per_blen[idx_blen];

    if(debug) cout << "node:" << np + 1 << " -> " << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , "  <<  nj + 1 << " , " << blen << endl;

    // loop over possible values of sk
    if(k == 2 * Ns){    // root node is always normal, for edge (Ns + 1, 2Ns)
        if(debug) cout << "Getting states for node MRCA " << k + 1 << endl;
        int sp = NORM_PLOIDY;
        if(model == BOUNDA) sp = NORM_ALLElE_STATE;
        if(debug) cout << "likelihood for state " << sp << endl;
        L_sk_k[k][sp] = get_max_prob_children(L_sk_k, S_sk_k, rtree, pblen, k, nstate, sp, ni, nj, blen, model, debug);
    }else{
        for(int sp = 0; sp < nstate; ++sp){  // looping over all possible states of its parent
            if(debug) cout << "likelihood for state " << sp << endl;
            L_sk_k[k][sp] = get_max_prob_children(L_sk_k, S_sk_k, rtree, pblen, k, nstate, sp, ni, nj, blen, model, debug);
        }
    }
  }

  if(debug){
    print_tree_lnl(rtree, L_sk_k, nstate);
    print_tree_state(rtree, S_sk_k, nstate);
  }
}



/** 
 * @brief Set the transition probability matrices for each branch length
 * @param rtree The evolutionary tree
 * @param Ns Number of internal nodes
 * @param nstate Number of states
 * @param model The model type
 * @param cn_max Maximum copy number
 * @param knodes The list of internal node indices
 * @param blens The list of unique branch lengths
 * @param pmat_per_blen The list of transition probability matrices corresponding to each branch length
 * @param fout The output file stream for logging
 */
void set_pmat(const evo_tree& rtree, int Ns, int nstate, int model, int cn_max, const vector<int>& knodes, vector<double>& blens, vector<double*>& pmat_per_blen, int debug){
  int dim_mat = nstate * nstate;
  double *qmat = new double[dim_mat];
  memset(qmat, 0.0, dim_mat * sizeof(double));

  assert(model == BOUNDT || model == BOUNDA);
  if(debug){
      cout << "Getting rate matrix" << endl;
  }
  if(model == BOUNDA){
      get_rate_matrix_haplotype_specific(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
  }else{
      get_rate_matrix_bounded(qmat, rtree.dup_rate, rtree.del_rate, cn_max);
  }

  for(int kn = 0; kn < knodes.size(); ++kn){
    int k = knodes[kn];
    double bli = rtree.edges[rtree.nodes[k].e_ot[0]].length;
    double blj = rtree.edges[rtree.nodes[k].e_ot[1]].length;

    if(find(blens.begin(), blens.end(), bli) == blens.end()){
      double *pmati = new double[dim_mat];
      memset(pmati, 0.0, dim_mat * sizeof(double));
      get_transition_matrix_bounded(qmat, pmati, bli, nstate);
      pmat_per_blen.push_back(pmati);
      blens.push_back(bli);
    }
    if(find(blens.begin(), blens.end(), blj) == blens.end()){
      double *pmatj = new double[dim_mat];
      memset(pmatj, 0.0, dim_mat * sizeof(double));
      get_transition_matrix_bounded(qmat, pmatj, blj, nstate);
      pmat_per_blen.push_back(pmatj);
      blens.push_back(blj);
    }
  }

  auto p = sort_permutation(blens, [&](const double a, const double b){ return a < b; });
  blens = apply_permutation(blens, p);
  pmat_per_blen = apply_permutation(pmat_per_blen, p);

  if(debug){
      print_branch_lengths(blens);

      for(int i = 0; i < pmat_per_blen.size(); ++i){
          double blen = blens[i];
          cout << "Pmatrix for branch length " << blen << endl;
          r8mat_print(nstate, nstate, pmat_per_blen[i], "  P matrix:");
      }
  }

  delete [] qmat;
}


string get_prob_line(const vector<vector<double>>& L_sk_k, int nid, int nchr, int nc, int model, int is_total, int cn_max){
    string line = to_string(nid + 1);
    if(nchr > 0){
        line += "\t" + to_string(nchr);
    }
    if(nc >= 0){
        line += "_" + to_string(nc + 1);
    }

    if(is_total && (model == BOUNDA || model == DECOMP)){   // need to convert state probability to cn probability under BOUNDA and DECOMP
        vector<double> Lsk_cn(cn_max + 1, 0.0);  // aggregate probabilities over the same total CN
        for(int i = 0; i < L_sk_k[nid].size(); i++){
            // cout << "\t" << to_string(L_sk_k[nid][i]);
            int cni = 0;
            if(model == BOUNDA){
                cni = state_to_total_cn(i, cn_max);
            }else{ // model == DECOMP
                int max_haplotype_change = cn_max / 2 - NORM_PLOIDY / 2;  // maximum haplotype-specific CN change given the maximum total CN change
                cni = change_state_to_total_cn(i, max_haplotype_change);  // for decomposition model, the state encodes the change information, need to convert to total CN change first and then get total CN by adding normal CN
            }
            Lsk_cn[cni] += L_sk_k[nid][i];
        }
        // change to relative probability
        double sum_prob = accumulate(Lsk_cn.begin(), Lsk_cn.end(), 0.0);
        for(int i = 0; i < Lsk_cn.size(); i++){
            double rprob = Lsk_cn[i] / sum_prob;
            line += "\t" + to_string(rprob);
        }
    }else{
        get_prob_line_orig(L_sk_k, nid, line);
    }

    return line;
}


// no need to sum over states for the same total CN
void get_prob_line_orig(const vector<vector<double>>& L_sk_k, int nid, string& line){
    double sum_prob = accumulate(L_sk_k[nid].begin(), L_sk_k[nid].end(), 0.0);
    for (int i = 0; i < L_sk_k[nid].size(); i++){
        double rprob = L_sk_k[nid][i] / sum_prob;
        line += "\t" + to_string(rprob);
    }
}

// get the most likely CN for a site nchr][nc] based on the state probabilities in L_sk_k, and store in cnp
void get_site_cnp(const vector<vector<double>>& L_sk_k, int nid, int nchr, int nc, int is_total, int cn_max, copy_number& cnp){
    // need to convert state probability to cn probability
    if(is_total){
        vector<double> Lsk_cn(cn_max + 1, 0.0);
        for(int i = 0; i < L_sk_k[nid].size(); i++){
            int cni = state_to_total_cn(i, cn_max);
            Lsk_cn[cni] += L_sk_k[nid][i];
        }
        // get the most likely total CN
        int tcn = distance(Lsk_cn.begin(), max_element(Lsk_cn.begin(), Lsk_cn.end()));
        cnp[nchr][nc] = tcn;
    }else{
        // get the index of the most likely haplotype-specific CN
        int state = distance(L_sk_k[nid].begin(), max_element(L_sk_k[nid].begin(), L_sk_k[nid].end()));
        cnp[nchr][nc] = state;
    }
}


void print_node_cnp(ofstream& fout, const copy_number& cnp, int nid, int cn_max, int is_total){
    for(auto site_cn : cnp){
        int nchr = site_cn.first;
        for(auto seg: site_cn.second){
            int nc = seg.first;
            // output CN for all sites: nid, chr, seg, cn
            string line_cn = to_string(nid + 1) + "\t" + to_string(nchr) + "\t" + to_string(nc + 1) + "\t";

            if(is_total){
                int tcn = seg.second;
                line_cn += to_string(tcn);
            }else{
                int cnA;
                int cnB;
                int state = seg.second;
                state_to_allele_cn(state, cn_max, cnA, cnB);
                line_cn += to_string(cnA) + "\t" + to_string(cnB);
            }
            fout << line_cn << endl;
        }
    }
}


// Infer the most likely CN change state for a site [nchr][nc] or a chr [nchr] based on the state probabilities in L_sk_k, and return the state index (for haplotype-specific model) or total CN (for total CN model)   
int get_node_change_state(const vector<vector<double>>& L_sk_k, int nid, int is_total, int max_haplotype_change){
    int state;
    if(is_total){
        int cn_max = (max_haplotype_change + 1) * 2;  // maximum total CN change given the maximum haplotype-specific change
        vector<double> Lsk_cn(cn_max + 1, 0.0);
        for(int i = 0; i < L_sk_k[nid].size(); i++){
            int cni = change_state_to_total_cn(i, max_haplotype_change);
            Lsk_cn[cni] += L_sk_k[nid][i];
        }
        // get the most likely total CN
        state = distance(Lsk_cn.begin(), max_element(Lsk_cn.begin(), Lsk_cn.end()));                    
    }else{                     
        state = distance(L_sk_k[nid].begin(), max_element(L_sk_k[nid].begin(), L_sk_k[nid].end()));;  // for haplotype-specific model, the state itself encodes the CN change information
    }
    return state;
}


string get_prob_line_header(int nstate){
    string header = "";
    for (int i = 0; i < nstate; i++)
    {
        header += "\tprobablity_" + to_string(i);
    }
    return(header);
}



double reconstruct_marginal_ancestral_state(const evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, const vector<int>& knodes, int model, int cn_max, int use_repeat, int is_total, string ofile){
    int debug = 0;
    if(debug) cout << "\treconstruct marginal ancestral state" << endl;

    string ofile_mrca = ofile + ".mrca.state";
    ofstream fout(ofile_mrca);

    string ofile_mrca_cn = ofile + ".mrca.cn";
    ofstream fout_cn(ofile_mrca_cn);

    int Ns = rtree.nleaf - 1;
    int nid = 2 * Ns;    // node ID for MRCA
    // For copy number instantaneous changes
    int nstate = cn_max + 1;
    if(model == BOUNDA) nstate = (cn_max + 1) * (cn_max + 2) / 2;

    string header="node\tsite";
    string prob_header = "";
    if(is_total && model == BOUNDA){
        prob_header = get_prob_line_header(cn_max + 1);
    }else{
        prob_header = get_prob_line_header(nstate);
    }
    fout << header << prob_header << endl;

    // Find the transition probability matrix for each branch
    vector<double> blens;
    vector<double*> pmat_per_blen;

    set_pmat(rtree, Ns, nstate, model, cn_max, knodes, blens, pmat_per_blen, debug);

    double logL = 0.0;    // for all chromosomes
    map<vector<int>, vector<vector<double>>> sites_lnl_map;
    copy_number cn_mrca;  // all CNs for MRCA

    // for each chromosome
    for(auto vcn : vobs){
      int nchr = vcn.first;
      if(debug) cout << "Computing likelihood for chr " << nchr << " with  " << vobs.at(nchr).size() << " sites" << endl;
      double site_logL = 0.0;   // log likelihood for all sites on a chromosome

      for(int nc = 0; nc < vobs.at(nchr).size(); nc++){    // for each segment on the chromosome
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs.at(nchr).at(nc);
          vector<vector<double>> L_sk_k(2 * rtree.nleaf - 1, vector<double>(nstate, 0.0));
          bool is_repeated = false;

          if(use_repeat){
              if(debug) cout << " Use repeated site patterns on site " << nc << endl;
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  if(debug) cout << "sites first seen" << endl;
                  initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);
                  get_likelihood_site(L_sk_k, rtree, knodes, blens, pmat_per_blen, 0, 0, model, nstate);
                  sites_lnl_map[obs] = L_sk_k;
              }else{
                  if(debug) cout << "sites repeated on site " << nc << endl;
                  L_sk_k = sites_lnl_map[obs];
                  is_repeated = true;
              }
          }else{
              initialize_lnl_table(L_sk_k, obs, rtree, model, nstate, is_total);
              get_likelihood_site(L_sk_k, rtree, knodes, blens, pmat_per_blen, 0, 0, model, nstate);
          }

          site_logL += extract_tree_lnl(L_sk_k, Ns, model);

          if(!is_repeated){
            string line = get_prob_line(L_sk_k, nid, nchr, nc, model, is_total, cn_max);
            // fout << setprecision(dbl::max_digits10) << line << endl;
            fout << line << endl;
          }
          get_site_cnp(L_sk_k, nid, nchr, nc, is_total, cn_max, cn_mrca);
      }

      logL += site_logL;

      if(debug){
          cout << "Site Likelihood for " << nchr << " is "  << site_logL << endl;
      }
    } // for each chromosome

    print_node_cnp(fout_cn, cn_mrca, nid, cn_max, is_total);

    fout.close();
    fout_cn.close();

    if(debug){
        cout << "\nLikelihood without correcting acquisition bias: " << logL << endl;
        // cout << "CNs at MRCA is: " << endl;
        // for(int i = 0; i < cn_mrca.size(); i++){
        //     cout << cn_mrca[i] << endl;
        // }
    }

    for_each(pmat_per_blen.begin(), pmat_per_blen.end(), DeleteObject());

    return logL;
}

// Infer the copy number of all internal nodes given a tree at a site, assuming only site duplication/deletion
// Using Pupko 2020 dynamic programming algorithm to reconstruct the joint ancestral state for all internal nodes
// Lx(i) is the likelihood of the best reconstruction of the subtree rooted at node x given that the parent of node x is in state i; 
// Sx(i) is the state of node x in the optimal conditional reconstruction 
void reconstruct_joint_ancestral_state(const evo_tree& rtree, const map<int, vector<vector<int>>>& vobs, vector<int>& knodes, int model, int cn_max, int use_repeat, int is_total, int m_max, string ofile){
    int debug = 0;
    if(debug) cout << "\treconstruct joint ancestral state" << endl;

    string ofile_joint = ofile + ".joint.state";
    ofstream fout(ofile_joint);

    string header;
    if(is_total){
        header = "node\tsite\tcn";
    }else{
        header = "node\tsite\tcnA\tcnB";
    }
    fout << header << endl;

    string ofile_joint_cn = ofile + ".joint.cn";
    ofstream fout_cn(ofile_joint_cn);

    int Ns = rtree.nleaf - 1;
    // For copy number instantaneous changes
    int nstate = cn_max + 1;
    if(model == BOUNDA) nstate = (cn_max + 1) * (cn_max + 2) / 2;
    int max_id = 2 * Ns;    // node ID for MRCA

    // Find the transition probability matrix for each branch
    vector<double> blens;
    vector<double*> pmat_per_blen;

    set_pmat(rtree, Ns, nstate, model, cn_max, knodes, blens, pmat_per_blen, debug);

    if(debug){
        print_branch_lengths(blens);
    }

    int ntotn = 2 * rtree.nleaf - 1;
    map<vector<int>, vector<vector<double>>> sites_lnl_map;   // ignore chr ID as the likelihood is the same
    map<vector<int>, vector<vector<int>>> sites_state_map;
    // map<vector<int>, int> sites_duplicated;
    map<int, copy_number> cnps;  // all CNs for all internal nodes

    // for each chromosome
    for(auto vcn : vobs){
      int nchr = vcn.first;
      if(debug) cout << "Computing likelihood on Chr " << nchr << " with " << vobs.at(nchr).size() << " sites" << endl;

      for(int nc = 0; nc < vobs.at(nchr).size(); nc++){    // for each segment on the chromosome
          if(debug) cout << "\tfor site " << nc << " on chromosome " << nchr << endl;
          // for each site of the chromosome (may be repeated)
          vector<int> obs = vobs.at(nchr).at(nc);
          vector<vector<double>> L_sk_k(ntotn, vector<double>(nstate, 0.0));
          vector<vector<int>> S_sk_k(ntotn, vector<int>(nstate, 0));
          bool is_repeated = false;

          if(use_repeat){
              if(debug) cout << " Use repeated site patterns on site " << nc << endl;
              if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
                  if(debug) cout << "sites first seen" << endl;
                  initialize_asr_table(obs, rtree, blens, pmat_per_blen, L_sk_k, S_sk_k, model, nstate, is_total);
                  get_ancestral_states_site(L_sk_k, S_sk_k, rtree, knodes, blens, pmat_per_blen, nstate, model);
                  sites_lnl_map[obs] = L_sk_k;
                  sites_state_map[obs] = S_sk_k;
                //   sites_duplicated[obs] = 1;
              }else{
                  if(debug) cout << "sites repeated on site " << nc << endl;
                //   sites_duplicated[obs]++;
                  L_sk_k = sites_lnl_map[obs];
                  S_sk_k =  sites_state_map[obs];
                  is_repeated = true;
              }
          }else{
              initialize_asr_table(obs, rtree, blens, pmat_per_blen, L_sk_k, S_sk_k, model, nstate, is_total);
              get_ancestral_states_site(L_sk_k, S_sk_k, rtree, knodes, blens, pmat_per_blen, nstate, model);
          }

          if(debug){
            cout << " Get the likelihood table of all internal nodes" << endl;
            print_tree_lnl(rtree, L_sk_k, nstate);
            print_tree_state(rtree, S_sk_k, nstate);
          }

        map<int, int> asr_states;     // The state ID used in rate matrix
        extract_tree_ancestral_state(rtree, knodes, S_sk_k, model, asr_states);

        if(!is_repeated){
            // int best_state = asr_states[Ns + 1];
            // cout << "optimal state at node "  << Ns + 1 << " is " << best_state << endl;
            // double prob = L_sk_k[max_id][best_state];  // posterior probability of optimal reconstruction
            // cout << "The probability vector at node MRCA " << max_id + 1 << ":";
            // for(int i = 0; i < nstate; i++){
            //     cout << "\t" <<  L_sk_k[max_id][i];
            // }
            // cout << endl;
            // double prob_sum = accumulate(L_sk_k[max_id].begin(), L_sk_k[max_id].end(), 0.0);
            // double prob_rel = prob / prob_sum;
            // cout << prob << "\t" << prob_sum << "\t" << prob_rel << endl;

            write_joint_state_line(max_id, Ns, asr_states, nchr, nc, is_total, cn_max, fout);
        }
        get_inode_cnp(max_id, Ns, asr_states, is_total, cn_max, cnps, nchr, nc);
      }
    } // for each chromosome

     for(int nid = max_id; nid > Ns + 1; nid--){
        print_node_cnp(fout_cn, cnps[nid], nid, cn_max, is_total);
     }

    fout.close();
    fout_cn.close();

    if(debug){
        if(use_repeat){
            write_pattern_all(vobs);
            write_pattern_uniq(sites_lnl_map);

            // cout << "write duplicated patterns" << endl;
            // ofstream fout2("./pattern_dup");
            // for(auto sd : sites_duplicated){
            //     fout2 << sd.second;
            //     for(auto c : sd.first){
            //         fout2 << "\t" << c;
            //     }
            //     fout2 << endl;
            // }
            // fout2.close();
        }
    }

    for_each(pmat_per_blen.begin(), pmat_per_blen.end(), DeleteObject());
}


void get_inode_cnp(int max_id, int Ns, map<int, int> &asr_states, int is_total, int cn_max, map<int, copy_number> &cnps, int nchr, int nc){
    // For all sites
    for (int nid = max_id; nid > Ns + 1; nid--){
        int state = asr_states[nid];

        if (is_total){
            int tcn = state_to_total_cn(state, cn_max);
            cnps[nid][nchr][nc] = tcn;
        }else{
            cnps[nid][nchr][nc] = state;
        }
    }
}


void write_joint_state_line(int max_id, int Ns, map<int, int>& asr_states, int nchr, int nc, int is_total, int cn_max, ofstream& fout){
    for (int nid = max_id; nid > Ns + 1; nid--){
        int state = asr_states[nid]; // state assigned to nid when its parent is optimal
        string line = to_string(nid + 1) + "\t" + to_string(nchr) + "_" + to_string(nc + 1);

        if (is_total){
            int tcn = state_to_total_cn(state, cn_max);
            line += "\t" + to_string(tcn);
        }else{
            int cnA, cnB;
            state_to_allele_cn(state, cn_max, cnA, cnB);
            line += "\t" + to_string(cnA) + "\t" + to_string(cnB);
        }

        // line += "\t" + to_string(prob_rel);

        // fout << setprecision(dbl::max_digits10) << line << endl;
        fout << line << endl;
    }
}


void write_pattern_uniq(map<vector<int>, vector<vector<double>>>& sites_lnl_map){
    cout << "write unique patterns" << endl;
    ofstream fout1("./pattern_uniq");
    for (auto sm : sites_lnl_map)
    {
        vector<int> key = sm.first;
        for (auto c : key)
        {
            fout1 << "\t" << c;
        }
        fout1 << endl;
    }
    fout1.close();
}

void write_pattern_all(const map<int, vector<vector<int>>>& vobs)
{
    ofstream fout3("./pattern_all");
    // for each chromosome
    for (auto vcn : vobs)
    {
        int nchr = vcn.first;
        for (int nc = 0; nc < vobs.at(nchr).size(); nc++)
        { // for each segment on the chromosome
            // for each site of the chromosome (may be repeated)
            vector<int> obs = vobs.at(nchr).at(nc);
            for (auto c : obs)
            {
                fout3 << "\t" << c;
            }
            fout3 << endl;
        }
    }
    fout3.close();
}



/*************** functions for independent chain model on multiple levels of CNAs *****************/


// print copy number changes for a node
void print_node_cnp_decomp(ofstream& fout, const copy_number_change& cnp, int nid, const LNL_TYPE& lnl_type, const DIM_DECOMP& dim_decomp){
    for(auto site_cn : cnp){
        int nchr = site_cn.first;
        for(auto seg: site_cn.second){
            int nc = seg.first;
            // output CN for all sites: nid, chr, seg, cn
            string line_cn = to_string(nid + 1) + "\t" + to_string(nchr) + "\t" + to_string(nc + 1) + "\t";

            CN_CHANGE cc = seg.second;

            if(lnl_type.is_total){               
                line_cn += cc.to_string();
            }else{
                vector<pair<int,int>> states_chr = build_pair_states(lnl_type.max_chr_change_haplotype);   
                pair<int, int> cnAB_chr = states_chr[cc.cn_change_chr];

                vector<pair<int,int>> states_site = build_pair_states(lnl_type.max_site_change_haplotype);   
                pair<int, int> cnAB_site = states_site[cc.cn_change_site];

                line_cn += "CN_CHANGE(cn_state=" + to_string(cc.cn_state);
                if(dim_decomp.dim_wgd > 0){
                    line_cn += ", num_wgd=" + to_string(cc.num_wgd); 
                }
                if(dim_decomp.dim_chr > 0){
                    line_cn += ", cn_change_chr=" + to_string(cnAB_chr.first) + "|" + to_string(cnAB_chr.second); 
                }
                if(dim_decomp.dim_seg > 0){
                    line_cn += ", cn_change_site=" + to_string(cnAB_site.first) + "|" + to_string(cnAB_site.second);
                }
                                      
                line_cn += ")";
            }
            fout << line_cn << endl;
        }
    }
}


double compute_site_likelihood(
    const evo_tree& rtree,
    const vector<CN_CHANGE>& obs,
    const DIM_DECOMP& dim_decomp,
    const LNL_TYPE& lnl_type,
    PMAT_DECOMP& pmat_decomp,   
    map<vector<CN_CHANGE>, vector<vector<double>>>& sites_lnl_map,
    vector<vector<double>>& lnl_table_seg,
    int debug)
{
    if(lnl_type.use_repeat){
        if(sites_lnl_map.find(obs) == sites_lnl_map.end()){
            if(debug) cout << "sites new" << endl;
            initialize_lnl_table_site(lnl_table_seg, rtree, obs, dim_decomp, lnl_type, debug);
            get_likelihood_site_change(lnl_table_seg, rtree, pmat_decomp, dim_decomp, lnl_type, debug);
            sites_lnl_map[obs] = lnl_table_seg;
        }else{
            if(debug) cout << "sites repeated" << endl;
            lnl_table_seg = sites_lnl_map[obs];
        }
    }else{
        if(debug) cout << "sites no repeat consideration" << endl;
        initialize_lnl_table_site(lnl_table_seg, rtree, obs, dim_decomp, lnl_type, debug);
        get_likelihood_site_change(lnl_table_seg, rtree, pmat_decomp, dim_decomp, lnl_type, debug);
    }

    double site_logL = extract_tree_lnl_change(lnl_table_seg, rtree.nleaf, NO_CHANGE_HAPLOTYPE, debug);

    return site_logL;
}


double compute_chr_likelihood(
    const evo_tree& rtree,
    int nchr,
    const vector<int>& change_chr,
    const vector<int>& knodes,
    const DIM_DECOMP& dim_decomp,
    const LNL_TYPE& lnl_type,
    PMAT_DECOMP& pmat_decomp,
    map<vector<int>, vector<vector<double>>>& chr_lnl_map,
    vector<vector<double>>& lnl_table_chr,
    int debug)
{
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

    if(debug) cout << "log likelihood for no change on chromosome " << nchr << " is: " << chr_logL << endl;

    return chr_logL;
}

int get_tcn_max_from_haplotype_change(int max_haplotype_change){
    return (max_haplotype_change + 1) * 2;  // maximum total CN change given the maximum haplotype-specific change
}


// Infer the copy number of the MRCA given a tree at a site, assuming independent Markov chains
// Get states for individual event type first, then figure out how to merge them
// knodes here does not equal to knodes in lnl_type
double reconstruct_marginal_ancestral_state_decomp(const evo_tree& rtree, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, const vector<int>& knodes, const OBS_DECOMP& obs_decomp, const LNL_TYPE& lnl_type_orig, string ofile, int debug){
    // int debug = 0;
    if(debug) cout << "\treconstruct marginal ancestral state with independent chain model" << endl;

    LNL_TYPE lnl_type = lnl_type_orig;
    lnl_type.knodes = knodes;   // update knodes in lnl_type for likelihood computation

    string ofile_mrca_cn = ofile + ".mrca.cn";
    ofstream fout_cn(ofile_mrca_cn);

    copy_number_change cn_mrca; // CNs for MRCA node

    int Ns = rtree.nleaf - 1; 
    int nid = 2 * (Ns + 1) - 2;   // the index of MRCA node in the tree, which is the last one in knodes

    double logL = 0.0;    // for all chromosomes

    DIM_DECOMP dim_decomp;
    set_pmat_decomp_dim(obs_decomp, dim_decomp, lnl_type, debug);

    QMAT_DECOMP qmat_decomp;
    if(debug) cout << "\tBuilding Q rate matrices for multiple levels" << endl;
    build_rate_matrices(qmat_decomp, rtree, dim_decomp, lnl_type, debug);

    if(debug) cout << "\tBuilding P transition matrices for multiple levels" << endl;
    PMAT_DECOMP pmat_decomp; 
    build_transition_matrices(pmat_decomp, rtree, knodes, qmat_decomp, dim_decomp, debug);

    string ofile_mrca = ofile + ".mrca.seg.state";
    ofstream fout_seg;
    if(dim_decomp.dim_seg > 0) fout_seg.open(ofile_mrca);

    ofile_mrca = ofile + ".mrca.chr.state";
    ofstream fout_chr;
    if(dim_decomp.dim_chr > 0) fout_chr.open(ofile_mrca);

    ofile_mrca = ofile + ".mrca.wgd.state";
    ofstream fout_wgd;
    if(dim_decomp.dim_wgd > 0) fout_wgd.open(ofile_mrca);

    string header="node\tsite";
    string prob_header = "";
    if(lnl_type.is_total){
        if(dim_decomp.dim_seg > 0){
            prob_header = get_prob_line_header(get_tcn_max_from_haplotype_change(lnl_type.max_site_change_haplotype) + 1);
            fout_seg << header << prob_header << endl;
        }
        if(dim_decomp.dim_chr > 0){
            prob_header = get_prob_line_header(get_tcn_max_from_haplotype_change(lnl_type.max_chr_change_haplotype) + 1);
            fout_chr << header << prob_header << endl; 
        }       
    }else{
        if(dim_decomp.dim_seg > 0){
            prob_header = get_prob_line_header(dim_decomp.dim_seg);
            fout_seg << header << prob_header << endl;
        }
        if(dim_decomp.dim_chr > 0){
            prob_header = get_prob_line_header(dim_decomp.dim_chr);
            fout_chr << header << prob_header << endl;
        }
    }
    
    if(dim_decomp.dim_wgd > 0){
        header="node";
        prob_header = get_prob_line_header(dim_decomp.dim_wgd);
        fout_wgd << header << prob_header << endl;
    }
  
    // Use a map to store computed log likelihood for each level of CNAs
    // recompute to save parameter change, also do sanity check   
    vector<int> sample_num_wgd(rtree.nleaf - 1, -1);         // use -1 to indicate uninitialized
    map<vector<CN_CHANGE>, vector<vector<double>>> sites_lnl_map;    // Use a map to store computed log likelihood
    map<vector<int>, vector<vector<double>>> chr_lnl_map;    // Use a map to store computed log likelihood

    // vector<vector<int>> mrca_cn_site_changes;  // record mrca states for each site on each chromosome
    // vector<int> mrca_cn_chr_changes;  // record mrca states for each site on each chromosome
    // int wgd_state = 0;
      
    // for each chromosome
    for(auto vcn : vobs_change){
        int nchr = vcn.first;
        if(debug) cout << "Computing likelihood on Chr " << nchr << endl;
        const vector<vector<CN_CHANGE>> obs_chr = vcn.second;   // all segments on this chromosome
        vector<int> change_chr(rtree.nleaf - 1, -100);   // total copy number change for each sample on this chromosome, use -100 to indicate uninitialized
        vector<vector<double>> lnl_table_chr;
        double site_logL = 0.0;   // log likelihood for all sites on a chromosome

        for(size_t nc = 0; nc < vobs_change.at(nchr).size(); nc++){    // for each segment on the chromosome
          // cout << "Number of sites for this chr " << vobs.at(nchr).size() << endl;
          // for each site of the chromosome (may be repeated)
            const vector<CN_CHANGE> obs = obs_chr.at(nc);
            vector<vector<double>> lnl_table_seg;     // one table for one site

            set_change_wgd_chr(obs, rtree, sample_num_wgd, change_chr);

            if(dim_decomp.dim_seg > 0){
                site_logL += compute_site_likelihood(rtree, obs, dim_decomp, lnl_type, pmat_decomp, sites_lnl_map, lnl_table_seg, debug);

                assert(lnl_table_seg.size() > nid);  // sanity check, the likelihood table should have the same number of nodes as the tree

                // Get the likelihood table of MRCA node (with largest ID) in the tree from likelihood table
                int site_state = get_node_change_state(lnl_table_seg, nid, lnl_type.is_total, lnl_type.max_site_change_haplotype);
                cn_mrca[nchr][nc].cn_change_site = site_state;

                // Print the state of MRCA at this site: sample, chromosome, segment, CN, probability of each change state (index for haplotype-specific CNs)
                int cn_max = get_tcn_max_from_haplotype_change(lnl_type.max_site_change_haplotype);  // maximum total CN change given the maximum haplotype-specific change
                string line = get_prob_line(lnl_table_seg, nid, nchr, nc, DECOMP, lnl_type.is_total, cn_max);
                // fout << setprecision(dbl::max_digits10) << line << endl;
                fout_seg << line << endl;

                if(debug){
                    cout << "\t\tSite state with maximum likelihood " << site_state << " at chr " << nchr << " site " << nc << endl;
                }                
            }
        }
        logL += site_logL;

        // compute likelihood for chr-level changes, which is the same for all sites on the chromosome, so only compute once
        if(dim_decomp.dim_chr > 0){
            if(debug){
                cout << "\nChr-level changes observed, computing chromosome-level likelihood for Chr " << nchr << endl;
            }
            double chr_logL = compute_chr_likelihood(rtree, nchr, change_chr, knodes, dim_decomp, lnl_type, pmat_decomp, chr_lnl_map, lnl_table_chr, debug);
            logL += chr_logL;   

            assert(lnl_table_chr.size() > nid);  // sanity check, the likelihood table should have the same number of nodes as the tree

            int chr_state = get_node_change_state(lnl_table_chr, nid, lnl_type.is_total, lnl_type.max_chr_change_haplotype);
            // reuse nc for consistency of indexing
            for(size_t nc = 0; nc < vobs_change.at(nchr).size(); nc++){ 
                cn_mrca[nchr][nc].cn_change_chr = chr_state;
            }
 
            int cn_max = get_tcn_max_from_haplotype_change(lnl_type.max_chr_change_haplotype);  // maximum total CN change given the maximum haplotype-specific change
            string line = get_prob_line(lnl_table_chr, nid, nchr, -1, DECOMP, lnl_type.is_total, cn_max);
            // fout << setprecision(dbl::max_digits10) << line << endl;
            fout_chr << line << endl;

            if(debug){
                cout << "\nChr " << nchr
                        << "  chr-level logL: "  << chr_logL
                        << "  site logL: "        << site_logL
                        << "  MRCA chr state: "   << chr_state << "\n";
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

        assert(lnl_table_wgd.size() > nid);  // sanity check, the likelihood table should have the same number of nodes as the tree

        // There may be multiple possible copy numbers given different ordering of events, so we need to compute the CN for each possible combination of events and then aggregate them to get the final CN for MRCA
        // int cn = pow(2, c[0] + 1) +  c[1] + c[2] + 2 * c[3] + 2 * c[4];
        // the index is the same as the number of WGDs
        int wgd_state = (int)distance(lnl_table_wgd[nid].begin(), max_element(lnl_table_wgd[nid].begin(), lnl_table_wgd[nid].end()));
        for(auto vcn : vobs_change){
            // reuse nchr, nc for consistency of indexing
            int nchr = vcn.first;
            for(size_t nc = 0; nc < vobs_change.at(nchr).size(); nc++){ 
                cn_mrca[nchr][nc].num_wgd = wgd_state;
            }           
        }
            string line = to_string(nid + 1);
            get_prob_line_orig(lnl_table_wgd, nid, line);
            // fout << setprecision(dbl::max_digits10) << line << endl;
            fout_wgd << line << endl;

        if(debug){
            cout << "\nWGD observed, computed WGD likelihood for all samples  " << wgd_logL << endl;
        }   
    }

    print_node_cnp_decomp(fout_cn, cn_mrca, nid, lnl_type, dim_decomp);

    fout_seg.close();
    fout_chr.close();
    fout_wgd.close();

    fout_cn.close();

    if(debug){
        cout << "\nLikelihood without correcting acquisition bias: " << logL << endl;
        print_cn_mrca(cn_mrca);
    }

    return logL;
}


void print_cn_mrca(copy_number_change& cn_mrca){
    cout << "CNs at MRCA is: " << endl;
    for (const auto &chr_entry : cn_mrca)
    {
        int nchr = chr_entry.first;
        for (const auto &seg_entry : chr_entry.second)
        {
            int nc = seg_entry.first;
            cout << "chr " << nchr
                 << " site " << nc + 1
                 << ": " << seg_entry.second << endl;
        }
    }
}

void set_change_wgd_chr(const vector<CN_CHANGE>& obs, const evo_tree& rtree, vector<int>& sample_num_wgd, vector<int>& change_chr)
{
    assert(obs.size() == rtree.nleaf - 1); // sanity check, the number of samples should be the same as the number of tips in the tree
    for (size_t i = 0; i < obs.size(); i++)
    {
        if (sample_num_wgd[i] == -1)
        {
            sample_num_wgd[i] = obs.at(i).num_wgd;
        }
        else
        {
            assert(sample_num_wgd[i] == obs.at(i).num_wgd); // sanity check, all sites should have the same WGD state for each sample
        }
        if (change_chr[i] == -100)
        {
            change_chr[i] = obs.at(i).cn_change_chr;
        }
        else
        {
            assert(change_chr[i] == obs.at(i).cn_change_chr); // sanity check, all sites should have the same WGD state for each sample
        }
    }
}


void set_tip_states_tcn_decomp(int state_chr, int max_change_haplotype, int& si, int& ei, vector<int>& tip_states){
    si = 0; // compute partial sum to get the index of the observed state
    ei = 0;
    int peak_sum_haplotype = max_change_haplotype + 2;
    get_tcn_state_index(state_chr, peak_sum_haplotype, si, ei);
    for (int k = si; k <= ei; k++){
        tip_states.push_back(k);
    }
}


// pmat_per_blen is a map with key as branch length and value as the corresponding transition probability matrix, which is a vector of size dim*dim in row-major order
void record_parent_state_decomp(int i, int dim, const vector<int>& tip_states, const map<double, double*>& pmat_per_blen, double blen, vector<vector<double>>& lnl_table, vector<vector<int>>& state_table, int debug){
    if(debug){
        cout << "recording parent states for node " << i + 1 << " with branch length " << blen << endl;
        // print_pmatrix_map(pmat_per_blen, dim, "X");
    }

    double* pmat_blen = pmat_per_blen.at(blen);
    if(debug > 1){
        r8mat_print(dim, dim, pmat_per_blen.at(blen), "  P matrix:");

        auto pmat_it = pmat_per_blen.find(blen);
        assert(pmat_it != pmat_per_blen.end());
        double* pmat = pmat_it->second;
        assert(pmat != nullptr);        
        r8mat_print(dim, dim, pmat, " Found P matrix:");
    }
    for(int j = 0; j < dim; ++j){  // For each possible parent state, find the most likely tip states
        vector<double> vec_lnl(dim, SMALL_LNL);
        // another loop as there maybe multiple states for a specific total CN
        for(int m = 0; m < tip_states.size(); ++m){
            int k = tip_states[m];
            if(debug > 1){
                cout << "parent state " << j << ", child state " << k << ", branch length " << blen << ", matrix dimensions " << dim << "x" << dim << endl;             
            }
            // use [] will add new entry if not exist, but it should exist as we have initialized the pmat for all possible branch lengths and states
            double li = pmat_blen[j + k * dim];  // assume parent has state j
            vec_lnl[k] = li;
        }
        if(debug > 1) print_vector<double>(vec_lnl);

        int max_i_state = distance(vec_lnl.begin(), max_element(vec_lnl.begin(), vec_lnl.end()));
        double max_lnl = *max_element(vec_lnl.begin(), vec_lnl.end());
        assert(max_lnl == vec_lnl[max_i_state]);

        if(debug) cout << "for node: i " << i + 1 << ", parent state " << j << ", max state is " << max_i_state << " with probability " << max_lnl << endl;
        state_table[i][j] = max_i_state;
        lnl_table[i][j] = max_lnl;
    }
}



void set_tip_states_decomp(const CN_CHANGE& obs_val, vector<int>& tip_states_chr, vector<int>& tip_states_site, const LNL_TYPE& lnl_type, int debug){   
    if(debug){
        cout << "Setting tip states for observed CN change: " << obs_val << endl;
    }  

    if(lnl_type.is_total){   // With total copy number, there may be multiple states corresponding to the same total CN, need to consider all of them
        int si = 0; // compute partial sum to get the index of the observed state
        int ei = 0;            
        if(obs_val.cn_change_chr != 0 && obs_val.cn_change_chr % CHANGE_CHR == 0){
            if(debug) cout << "Ambiguous encoding for chromosome gain/loss, including cases either before or after WGD." << endl;
            
            // after WGD, -2 was normalized to -1
            int state_chr = obs_val.cn_change_chr / CHANGE_CHR - MIN_CHANGE;    // CN normalied by ploidy previously           
            set_tip_states_tcn_decomp(state_chr, lnl_type.max_chr_change_haplotype, si, ei, tip_states_chr);   
                        
            // before WGD, -2
            state_chr = (state_chr + MIN_CHANGE) * NORM_PLOIDY - MIN_CHANGE;   // add WGD back
            set_tip_states_tcn_decomp(state_chr, lnl_type.max_chr_change_haplotype, si, ei, tip_states_chr);
        }else{
            int state_chr = obs_val.cn_change_chr - MIN_CHANGE;        
            set_tip_states_tcn_decomp(state_chr, lnl_type.max_chr_change_haplotype, si, ei, tip_states_chr);     
        } 

        int state_site = obs_val.cn_change_site - MIN_CHANGE; 
        set_tip_states_tcn_decomp(state_site, lnl_type.max_site_change_haplotype, si, ei, tip_states_site);
    }else{ // With haplotype-specific copy number, only the specific site needs to be filled, size 1
        if(obs_val.cn_change_chr % CHANGE_CHR == 0){                
            if(debug) cout << "Ambiguous encoding for chromosome gain/loss, including cases either before or after WGD." << endl;
            int state_chr = obs_val.cn_change_chr / CHANGE_CHR;    // CN normalied by ploidy previously
            tip_states_chr.push_back(state_chr);
            // after WGD, +2, 2 concecutive changes, not supported for now
            // need to distinguish changes on haplotypes to get exact haplotype-specific copy numbers
        }else{    
            tip_states_chr.push_back(obs_val.cn_change_chr);
        }            
        tip_states_chr.push_back(obs_val.cn_change_chr);

        tip_states_site.push_back(obs_val.cn_change_site);
    }
}


/** 
 * @brief Initialize the tables of tip nodes for reconstructing joint ancestral states under independent chain model (Step 1 of Pupko algorithm)
 * state_table: vectors to store tip states
 * pmat_decomp: multiple P matrices for each branches for multiple CNA types
*/
void initialize_asr_table_decomp(const vector<CN_CHANGE>& obs, const evo_tree& rtree, const PMAT_DECOMP& pmat_decomp, const DIM_DECOMP& dim_decomp, LNL_TABLE& lnl_table, STATE_TABLE& state_table, const LNL_TYPE& lnl_type, int debug){
    if(debug){
        cout << "Initializing tables for reconstructing joint ancestral state under independent chain model" << endl;
    }

    for(int i = 0; i < rtree.nleaf - 1; ++i){
        double blen = rtree.edges[rtree.nodes[i].e_in].length;

        if(debug > 1) {
            cout << "node " << i + 1  << ", observed CN " << obs[i] << endl;
            int parent = rtree.edges[rtree.nodes[i].e_in].start;
            cout << "parent " << parent + 1 << ", blen " << blen << endl;
        }
  
        // States recorded by cn change for total copy number and state index for haplotype-specific copy number
        int tip_state_wgd = obs[i].num_wgd;
        vector<int> tip_states_chr;
        vector<int> tip_states_site;   
        set_tip_states_decomp(obs[i], tip_states_chr, tip_states_site, lnl_type, debug);

        if(dim_decomp.dim_wgd > 0){
            if(debug) cout << "Setting tip states for WGD change: " << tip_state_wgd << endl;
            record_parent_state_decomp(i, dim_decomp.dim_wgd, vector<int>{tip_state_wgd}, pmat_decomp.pmats_wgd, blen, lnl_table.lnl_table_wgd, state_table.state_table_wgd, debug);
        }
        if(dim_decomp.dim_chr > 0){
            if (debug){
                cout << "Setting tip states for chromosome change: ";
                print_vector<int>(tip_states_chr);
            }           
            record_parent_state_decomp(i, dim_decomp.dim_chr, tip_states_chr, pmat_decomp.pmats_chr, blen, lnl_table.lnl_table_chr, state_table.state_table_chr, debug);
        }
        if(dim_decomp.dim_seg > 0){
            if (debug){
                cout << "Setting tip states for site change: ";
                print_vector<int>(tip_states_site);
            }
            record_parent_state_decomp(i, dim_decomp.dim_seg, tip_states_site, pmat_decomp.pmats_seg, blen, lnl_table.lnl_table_seg, state_table.state_table_seg, debug);
        }
    }

    if(debug){
      print_lnl_at_tips<CN_CHANGE>(rtree, obs, lnl_table.lnl_table_seg, dim_decomp.dim_seg);
      print_tip_states(rtree, dim_decomp.dim_seg, state_table.state_table_seg);

      print_lnl_at_tips<CN_CHANGE>(rtree, obs, lnl_table.lnl_table_chr, dim_decomp.dim_chr);
      print_tip_states(rtree, dim_decomp.dim_chr, state_table.state_table_chr);

      print_lnl_at_tips<CN_CHANGE>(rtree, obs, lnl_table.lnl_table_wgd, dim_decomp.dim_wgd);
      print_tip_states(rtree, dim_decomp.dim_wgd, state_table.state_table_wgd);     
    }
}



// Find the state of a node k that gives the maximum likelihood under the given parent state sp, and record the state in S_sk_k
// return the maximum likelihood for node k given parent state sp
void get_max_prob_children_decomp(LNL_TABLE& lnl_table, STATE_TABLE& state_table, const evo_tree& rtree, const PMAT_DECOMP& pmat_decomp, int k, const DIM_DECOMP& dim_decomp, const CN_CHANGE& sp, int ni, int nj, double blen, int model, int debug){
    if(debug){
        cout << "Getting max probability under independent model for node " << k + 1 << " with parent state : " << sp << endl;
    }

    if(dim_decomp.dim_seg > 0){
        lnl_table.lnl_table_seg[k][sp.cn_change_site] = get_max_prob_children(lnl_table.lnl_table_seg, state_table.state_table_seg, rtree, pmat_decomp.pmats_seg.at(blen), k, dim_decomp.dim_seg, sp.cn_change_site, ni, nj, blen, DECOMP, debug);
    }

    if(dim_decomp.dim_chr > 0){
        lnl_table.lnl_table_chr[k][sp.cn_change_chr] = get_max_prob_children(lnl_table.lnl_table_chr, state_table.state_table_chr, rtree, pmat_decomp.pmats_chr.at(blen), k, dim_decomp.dim_chr, sp.cn_change_chr, ni, nj, blen, DECOMP, debug);
    }
    
    if(dim_decomp.dim_wgd > 0){
        lnl_table.lnl_table_wgd[k][sp.num_wgd] = get_max_prob_children(lnl_table.lnl_table_wgd, state_table.state_table_wgd, rtree, pmat_decomp.pmats_wgd.at(blen), k, dim_decomp.dim_wgd, sp.num_wgd, ni, nj, blen, DECOMP, debug);
    }
}



/**
 * @brief Get the ancestral states site for nonroot internal nodes under independent Markov chain model (Step 2 of Pupko algorithm)
 * For each node, get the state with maximum likelihood given the parent state, and record the state in the state table for each type of CNAs
 */
void get_ancestral_states_site_decomp(LNL_TABLE& lnl_table, STATE_TABLE& state_table, const evo_tree& rtree, const vector<int>& knodes, const PMAT_DECOMP& pmat_decomp, const DIM_DECOMP& dim_decomp, int debug){
  if(debug){
      cout << "Getting ancestral state for one site under independent Markov chain model" << endl;
      cout << dim_decomp.dim_wgd << "\t" << dim_decomp.dim_chr << "\t"  << dim_decomp.dim_seg << endl;
  }

  for(int kn = 0; kn < knodes.size() - 1; ++kn){
        int k = knodes[kn];
        int np = rtree.edges[rtree.nodes[k].e_in].start;
        int ni = rtree.edges[rtree.nodes[k].e_ot[0]].end;
        int nj = rtree.edges[rtree.nodes[k].e_ot[1]].end;
        double blen = rtree.edges[rtree.nodes[k].e_in].length;

        if(debug) cout << "node:" << np + 1 << " -> " << rtree.nodes[k].id + 1 << " -> " << ni + 1 << " , "  <<  nj + 1 << " , " << blen << endl;

        // loop over possible observed states of start nodes
        int Ns = rtree.nleaf - 1;
        if(k == 2 * Ns){    // root node is always normal, for edge (Ns + 1, 2Ns)
            CN_CHANGE sp {NORM_PLOIDY, 0, NO_CHANGE_HAPLOTYPE, NO_CHANGE_HAPLOTYPE};
            if(debug) cout << "Getting states for MRCA node " << k + 1 << " with state " << sp << endl;
            get_max_prob_children_decomp(lnl_table, state_table, rtree, pmat_decomp, k, dim_decomp, sp, ni, nj, blen, DECOMP, debug);
        }else{  
            if(debug) cout << "Getting multiple-level states for node " << k + 1 << endl; 
            for(int sp = 0; sp < dim_decomp.dim_seg; ++sp){               
                lnl_table.lnl_table_seg[k][sp] = get_max_prob_children(lnl_table.lnl_table_seg, state_table.state_table_seg, rtree, pmat_decomp.pmats_seg.at(blen), k, dim_decomp.dim_seg, sp, ni, nj, blen, DECOMP, debug);
            }

            for(int sp = 0; sp < dim_decomp.dim_chr; ++sp){
                lnl_table.lnl_table_chr[k][sp] = get_max_prob_children(lnl_table.lnl_table_chr, state_table.state_table_chr, rtree, pmat_decomp.pmats_chr.at(blen), k, dim_decomp.dim_chr, sp, ni, nj, blen, DECOMP, debug);
            }

            for(int sp = 0; sp < dim_decomp.dim_wgd; ++sp){
                lnl_table.lnl_table_wgd[k][sp] = get_max_prob_children(lnl_table.lnl_table_wgd, state_table.state_table_wgd, rtree, pmat_decomp.pmats_wgd.at(blen), k, dim_decomp.dim_wgd, sp, ni, nj, blen, DECOMP, debug);
            }
        }
  }
  if(debug){
    cout << "\nFinished getting ancestral states for all nonroot internal nodes, printing the likelihood and states for each level of CNAs" << endl;
    
    cout << "SEG level:" << endl;
    print_tree_lnl(rtree, lnl_table.lnl_table_seg, dim_decomp.dim_seg);
    print_tree_state(rtree, state_table.state_table_seg, dim_decomp.dim_seg);

    cout << "CHR level: " << endl;
    print_tree_lnl(rtree, lnl_table.lnl_table_chr, dim_decomp.dim_chr);
    print_tree_state(rtree, state_table.state_table_chr, dim_decomp.dim_chr);

    cout << "WGD level: " << endl;
    print_tree_lnl(rtree, lnl_table.lnl_table_wgd, dim_decomp.dim_wgd);
    print_tree_state(rtree, state_table.state_table_wgd, dim_decomp.dim_wgd);
  }
}

// TODO: to check validity and update
static pair<int, int> decomp_state_to_change_pair(int state, int max_change_haplotype) {
    if (state == NO_CHANGE_HAPLOTYPE) {
        return {0, 0};
    }

    int width = max_change_haplotype + 2;
    int a = state / width + MIN_CHANGE_HAPLOTYPE;
    int b = state % width + MIN_CHANGE_HAPLOTYPE;
    return {a, b};
}



int get_cn_from_state_decomp(int state_wgd,
                             int state_chr,
                             int state_site,
                             const LNL_TYPE& lnl_type) {
    int baseline_hap = 1;
    for (int i = 0; i < state_wgd; ++i) {
        baseline_hap *= NORM_PLOIDY;
    }

    pair<int, int> chr_change =
        decomp_state_to_change_pair(state_chr, lnl_type.max_chr_change_haplotype);
    pair<int, int> site_change =
        decomp_state_to_change_pair(state_site, lnl_type.max_site_change_haplotype);

    int cnA = max(0, baseline_hap + chr_change.first + site_change.first);
    int cnB = max(0, baseline_hap + chr_change.second + site_change.second);

    if (lnl_type.is_total) {
        return cnA + cnB;
    }

    return allele_cn_to_state(cnA, cnB);
}


bool has_seg_state(const vector<int>& state_table_nid) {
    bool has_seg_state = false;
    // print_vector<int>(state_table_nid);

    for (int state : state_table_nid) {
        if (state != NONE_STATE) {
            has_seg_state = true;
            break;
        }
    }
    // cout << "has_seg_state: " << has_seg_state << endl;
    return has_seg_state;
}


void extract_tree_ancestral_state_decomp(const evo_tree& rtree, const vector<int>& knodes, const STATE_TABLE& state_table, map<int, CN_CHANGE>&asr_states, const LNL_TYPE& lnl_type, int debug){
    if(debug){
        cout << "Exatract most likely joint estimation of ancestral nodes under independent chain model" << endl;
    }

    CN_CHANGE parent_state = {NORM_PLOIDY, 0, NO_CHANGE_HAPLOTYPE, NO_CHANGE_HAPLOTYPE};
    asr_states[rtree.nleaf] = parent_state;   // for root, adaption of step 4 of Pupko algorithm

    // Step 5 of Pupko algorithm
    // Traverse the tree from root to tips, need to know the parent of each node
    for(int i = knodes.size() - 2; i >= 0; i--){  // starting from node ID for MRCA
        int nid = knodes[i];
        // Find the parent node of current node
        int parent = rtree.nodes[nid].parent;
        if(asr_states.find(parent) == asr_states.end()){
            cout << "Cannot find state for the parent of node " << nid + 1 << ", " << parent + 1 << endl;
            exit(EXIT_FAILURE);
        }
        parent_state = asr_states[parent];
        if(debug) cout << "parent_state of node " << nid + 1 << ", " << parent + 1 << " is " << parent_state << endl;

        int state_site = NO_CHANGE_HAPLOTYPE;
        int state_chr = NO_CHANGE_HAPLOTYPE;
        int state_wgd = 0;

        if(has_seg_state(state_table.state_table_seg[nid])){
            if(debug) print_vector<int>(state_table.state_table_seg[nid]);
            assert(parent_state.cn_change_site >= 0 && parent_state.cn_change_site < state_table.state_table_seg[nid].size());
            assert(state_table.state_table_seg[nid][parent_state.cn_change_site] != NONE_STATE);
            state_site = state_table.state_table_seg[nid][parent_state.cn_change_site];
        }

        if(has_seg_state(state_table.state_table_chr[nid])){
            if(debug) print_vector<int>(state_table.state_table_chr[nid]);
            assert(parent_state.cn_change_chr >= 0 && parent_state.cn_change_chr < state_table.state_table_chr[nid].size());
            assert(state_table.state_table_chr[nid][parent_state.cn_change_chr] != NONE_STATE);
            state_chr = state_table.state_table_chr[nid][parent_state.cn_change_chr];
        }

        if(has_seg_state(state_table.state_table_wgd[nid])){
            if(debug) print_vector<int>(state_table.state_table_wgd[nid]);
            assert(parent_state.num_wgd >= 0 && parent_state.num_wgd < state_table.state_table_wgd[nid].size());
            assert(state_table.state_table_wgd[nid][parent_state.num_wgd] != NONE_STATE);
            state_wgd = state_table.state_table_wgd[nid][parent_state.num_wgd];
        }

        // total CN: relative CN change 
        // haplotype-specific CN: state index for haplotype-specific CN
        int cn = get_cn_from_state_decomp(state_wgd, state_chr, state_site, lnl_type);

        // has to follow the order of variables
        CN_CHANGE state = {cn, state_wgd, state_chr, state_site};
        asr_states[nid] = state;

        if(debug){
            cout << "\t\tnode " << nid + 1 << " with state " << state << " and parent " << parent + 1 << " whose state is " << parent_state << endl;
        }
    }

}


void get_inode_cnp_decomp(int max_id, int Ns, map<int, CN_CHANGE>& asr_states, const LNL_TYPE& lnl_type, map<int, copy_number_change>& cnps, int nchr, int nc){
    // For all sites
    for (int nid = max_id; nid > Ns + 1; nid--){
        CN_CHANGE state = asr_states[nid];
        // TODO: compute abosolute copy numbers
        cnps[nid][nchr][nc] = state;       
    }
}


void write_joint_state_line_decomp(int max_id, int Ns, map<int, CN_CHANGE>& asr_states, int nchr, int nc, const LNL_TYPE& lnl_type, ofstream& fout){
    for (int nid = max_id; nid > Ns + 1; nid--){
        CN_CHANGE state = asr_states[nid]; // state assigned to nid when its parent is optimal
        string line = to_string(nid + 1) + "\t" + to_string(nchr) + "_" + to_string(nc + 1) + "\t" + state.to_string();
        // fout << setprecision(dbl::max_digits10) << line << endl;
        fout << line << endl;
    }
}


// Infer the copy number of all internal nodes given a tree at a site, assuming independent Markov chains
// Similar to the dynamic programming of likelihood computation, store the most likely state for each node at each site and then extract the ancestral state for all internal nodes by backtracking from the root to the tips
// Replace “sum over child states” with “max over child states and remember the argmax”
void reconstruct_joint_ancestral_state_decomp(const evo_tree& rtree, const map<int, vector<vector<CN_CHANGE>>>& vobs_change, vector<int>& knodes, const OBS_DECOMP& obs_decomp, const LNL_TYPE& lnl_type_orig, string ofile, int debug){
    if(debug) cout << "\treconstruct joint ancestral state under independent chain model" << endl;

    LNL_TYPE lnl_type = lnl_type_orig;
    lnl_type.knodes = knodes;   // update knodes in lnl_type for likelihood computation
    if(debug){
        cout << "\tknodes for ASR: " << knodes.size() << " nodes, from " << knodes.front() + 1 << " to " << knodes.back() + 1 << endl;
        print_vector<int>(knodes);
    }

    string ofile_joint = ofile + ".joint.state";
    ofstream fout(ofile_joint);

    string header = "node\tsite\tcn_change";
    fout << header << endl;

    // no header to be consistent the main input copy number file, which also does not have header
    string ofile_joint_cn = ofile + ".joint.cn";
    ofstream fout_cn(ofile_joint_cn);

    int Ns = rtree.nleaf - 1;
    int max_id = 2 * Ns;   // node ID for MRCA

    int ntotn = 2 * rtree.nleaf - 1;

    DIM_DECOMP dim_decomp;
    set_pmat_decomp_dim(obs_decomp, dim_decomp, lnl_type, debug);

    QMAT_DECOMP qmat_decomp;
    build_rate_matrices(qmat_decomp, rtree, dim_decomp, lnl_type, debug);

    if(debug) cout << "\tBuilding P transition matrices for multiple levels" << endl;
    PMAT_DECOMP pmat_decomp; 
    build_transition_matrices(pmat_decomp, rtree, lnl_type.knodes, qmat_decomp, dim_decomp, debug);

    // vector<int> sample_num_wgd(rtree.nleaf - 1, -1);         // use -1 to indicate uninitialized
    // vector<vector<int>> chr_sample_change;   // indexed by chromosome first and then sample, to facilitate likelihood computation
    map<vector<CN_CHANGE>, LNL_TABLE> sites_lnl_map;    // Use a map to store computed log likelihood
    map<vector<CN_CHANGE>, STATE_TABLE> sites_state_map;
    // map<vector<int>, vector<vector<double>>> chr_lnl_map;    // Use a map to store computed log likelihood
    map<int, copy_number_change> cnps;  // all CNs for all internal nodes

    // similar to the likelihood computation, but store the state for each node at each site in a separate table, and then extract the ancestral state for all internal nodes by backtracking from the root to the tips
    // do calculation for CNAs at different levels separately, and then combine them to get the final CN for each node
    for(auto vcn :  vobs_change){
      int nchr = vcn.first;
      if(debug) cout << "Computing likelihood on Chr " << nchr << " with " << vobs_change.at(nchr).size() << " sites " << endl;

      for(int nc = 0; nc < vobs_change.at(nchr).size(); nc++){    // for each segment on the chromosome
        vector<CN_CHANGE> obs = vobs_change.at(nchr).at(nc);
        if(debug){
            cout << "\tfor site " << nc << " on chromosome " << nchr << endl;
            print_vector<CN_CHANGE>(obs);
        }

        // A state/likelihood table for each level of CNAs at each site
        LNL_TABLE lnl_table;
        STATE_TABLE state_table;
        initialize_lnl_state_table(lnl_table, state_table, ntotn, dim_decomp);

        bool is_repeated = false;

        if(lnl_type.use_repeat){
            if (sites_lnl_map.find(obs) == sites_lnl_map.end()){
                initialize_asr_table_decomp(obs, rtree, pmat_decomp, dim_decomp, lnl_table, state_table, lnl_type, debug);
                get_ancestral_states_site_decomp(lnl_table, state_table, rtree, knodes, pmat_decomp, dim_decomp, debug);
                sites_lnl_map[obs] = lnl_table;
                sites_state_map[obs] = state_table;
            }else{
                if (debug) cout << "\t\tsites repeated" << endl;
                lnl_table = sites_lnl_map[obs];
                state_table = sites_state_map[obs];
                is_repeated = true;
            }
        }else{
            initialize_asr_table_decomp(obs, rtree, pmat_decomp, dim_decomp, lnl_table, state_table, lnl_type, debug);
            get_ancestral_states_site_decomp(lnl_table, state_table, rtree, knodes, pmat_decomp, dim_decomp, debug);
        }

        map<int, CN_CHANGE> asr_states;
        extract_tree_ancestral_state_decomp(rtree, knodes, state_table, asr_states, lnl_type, debug);

        if(!is_repeated){
            if(debug) cout << "\t\twriting joint state for Chr " << nchr << " site " << nc << endl;
            write_joint_state_line_decomp(max_id, Ns, asr_states, nchr, nc, lnl_type, fout);
        }
        if(debug) cout << "\t\tgetting joint CN for Chr " << nchr << " site " << nc << endl;
        get_inode_cnp_decomp(max_id, Ns, asr_states, lnl_type, cnps, nchr, nc);
      }
    } // for each chromosome

     if(debug) cout << "\twriting joint CN for all internal nodes" << endl;
     for(int nid = max_id; nid > Ns + 1; nid--){
        if(debug) cout << "\t\tfor node " << nid + 1 << endl;
        print_node_cnp_decomp(fout_cn, cnps[nid], nid, lnl_type, dim_decomp);
     }

    fout.close();
    fout_cn.close();
}


void initialize_lnl_state_table(LNL_TABLE& lnl_table, STATE_TABLE& state_table, int ntotn, DIM_DECOMP& dim_decomp)
{
    lnl_table.lnl_table_wgd.assign(ntotn, vector<double>(dim_decomp.dim_wgd, 0.0));
    lnl_table.lnl_table_chr.assign(ntotn, vector<double>(dim_decomp.dim_chr, 0.0));
    lnl_table.lnl_table_seg.assign(ntotn, vector<double>(dim_decomp.dim_seg, 0.0));

    // use NONE_STATE to check if values are changed
    state_table.state_table_wgd.assign(ntotn, vector<int>(dim_decomp.dim_wgd, NONE_STATE));
    state_table.state_table_chr.assign(ntotn, vector<int>(dim_decomp.dim_chr, NONE_STATE));
    state_table.state_table_seg.assign(ntotn, vector<int>(dim_decomp.dim_seg, NONE_STATE));
}
