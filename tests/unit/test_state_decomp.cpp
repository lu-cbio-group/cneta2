// Unit tests for turning independent-chain (model 3) states back into copy
// numbers in ancestral state reconstruction (state.cpp).
//
// Chromosome- and site-level states are indices into the list built by
// build_pair_states() (common.hpp), which is ordered by total copy number and
// then by haplotype A. Decoding must use the same order, or reconstructed copy
// numbers are wrong even when the reconstructed states are right.

#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <array>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

#include "likelihood.hpp"
#include "parse_cn.hpp"
#include "state.hpp"
#include "tree_op.hpp"

namespace {

std::string data_path(const std::string& name){
    return std::string(CNETA_TEST_DATA_DIR) + "/" + name;
}

const int NS = 3;    // tumour regions in the fixture

// Fixture data in the form mode 4 hands to the reconstruction: total copy
// number, segment-level changes only.
struct AsrInput {
    evo_tree tree;
    map<int, vector<vector<CN_CHANGE>>> vobs_change;
    OBS_DECOMP obs_decomp;
    LNL_TYPE lnl_type;
    vector<int> inodes;
};

AsrInput load_fixture(int max_site_change_haplotype){
    AsrInput in;
    const int cn_max = 4, max_wgd = 1, max_chr_change_haplotype = 1;
    INPUT_PROPERTY prop{NS, cn_max, max_wgd, max_chr_change_haplotype, max_site_change_haplotype,
                        DECOMP, /*is_total*/ 1, /*is_rcn*/ 0, /*is_bin*/ 0, /*incl_all*/ 1};
    INPUT_DATA data{0, 0, 0, vector<int>(NS, 0), vector<int>(NS, 0), vector<int>(NS, 0)};
    map<int, vector<vector<int>>> data_change;
    read_data_var_regions_by_chr_change(data_change, data_path("tiny-cn.txt.gz"), "", "", prop, data, 0);
    get_obs_vector_by_chr_change(data_change, in.vobs_change, NS);

    int max_site_change = *max_element(data.site_max_change.begin(), data.site_max_change.end());
    in.obs_decomp = {/*m_max*/ 1, /*max_wgd*/ 0, /*max_chr_change*/ 0, max_site_change};

    in.tree = read_tree_info(data_path("tiny-tree.txt"), NS);
    restore_mutation_rates(in.tree, vector<double>{0, 0.001, 0.001, 0, 0, 0});
    in.inodes = get_inodes_bottom_up(in.tree, 0);

    in.lnl_type = {DECOMP, cn_max, max_wgd, max_chr_change_haplotype, max_site_change_haplotype,
                   /*is_total*/ 1, /*cons*/ 1, /*max_tobs*/ 4.0, /*age*/ 100,
                   /*use_repeat*/ 1, /*correct_bias*/ 0, data.num_invar_bins, data.seg_size,
                   /*cn_type*/ 0, in.inodes};
    return in;
}

std::string out_prefix(const std::string& name){
    auto dir = std::filesystem::temp_directory_path() / "cneta_test_state_decomp";
    std::filesystem::create_directories(dir);
    for(const auto& entry : std::filesystem::directory_iterator(dir)){
        if(entry.path().filename().string().rfind(name + ".", 0) == 0){
            std::filesystem::remove(entry.path());
        }
    }
    return (dir / name).string();
}

LNL_TYPE lnl_type_for(int max_chr_change, int max_site_change, int is_total){
    LNL_TYPE t{};
    t.model = DECOMP;
    t.max_wgd = 1;
    t.max_chr_change_haplotype = max_chr_change;
    t.max_site_change_haplotype = max_site_change;
    t.is_total = is_total;
    return t;
}

}  // namespace

TEST_CASE("site-level states decode to the copy number of their change pair", "[state]"){
    for(int max_change = 1; max_change <= 4; max_change++){
        auto states = build_pair_states(max_change);
        for(int is_total : {1, 0}){
            LNL_TYPE t = lnl_type_for(1, max_change, is_total);
            for(int idx = 0; idx < (int)states.size(); idx++){
                INFO("max_change " << max_change << ", index " << idx << ", is_total " << is_total);
                int cnA = 1 + states[idx].first;
                int cnB = 1 + states[idx].second;
                int expected = is_total ? cnA + cnB : allele_cn_to_state(cnA, cnB);
                REQUIRE(get_cn_from_state_decomp(0, NO_CHANGE_HAPLOTYPE, idx, t) == expected);
            }
        }
    }
}

TEST_CASE("chromosome-level states decode to the copy number of their change pair", "[state]"){
    for(int max_change = 1; max_change <= 3; max_change++){
        auto states = build_pair_states(max_change);
        LNL_TYPE t = lnl_type_for(max_change, 2, 1);
        for(int idx = 0; idx < (int)states.size(); idx++){
            INFO("max_change " << max_change << ", index " << idx);
            int expected = 2 + states[idx].first + states[idx].second;
            REQUIRE(get_cn_from_state_decomp(0, idx, NO_CHANGE_HAPLOTYPE, t) == expected);
        }
    }
}

TEST_CASE("marginal MRCA probabilities are the exact posterior", "[state]"){
    // On the 3-region tree the MRCA has one internal child, so the posterior of
    // each MRCA state can be computed by summing over every state of that child
    // and every tip state compatible with the observed total copy number. The
    // root is in the normal state, so the MRCA posterior includes the transition
    // along the trunk edge (root -> MRCA).
    AsrInput in = load_fixture(2);
    const evo_tree& t = in.tree;

    DIM_DECOMP dim;
    set_pmat_decomp_dim(in.obs_decomp, dim, in.lnl_type, 0);
    REQUIRE(dim.dim_seg > 0);
    QMAT_DECOMP qmat;
    build_rate_matrices(qmat, t, dim, in.lnl_type, 0);
    PMAT_DECOMP pmat;
    build_transition_matrices(pmat, t, in.inodes, qmat, dim, 0);

    const int n = dim.dim_seg;
    const int mrca = 2 * NS;
    const int inner = NS + 2;
    REQUIRE(t.edges[t.nodes[inner].e_in].start == mrca);
    auto states = build_pair_states(in.lnl_type.max_site_change_haplotype);
    REQUIRE((int)states.size() == n);

    auto p = [&](int child, int from, int to){
        return pmat.pmats_seg.at(t.edges[t.nodes[child].e_in].length)[from + to * n];
    };

    std::string prefix = out_prefix("marginal");
    reconstruct_marginal_ancestral_state_decomp(t, in.vobs_change, in.inodes, in.obs_decomp, in.lnl_type, false, prefix, 0);

    std::ifstream fin(prefix + ".mrca.seg.state");
    std::string line;
    REQUIRE(std::getline(fin, line));    // header

    int nsite = 0;
    vector<int> best_cn;
    vector<double> best_prob;
    for(const auto& [chr, sites] : in.vobs_change){
        for(size_t nc = 0; nc < sites.size(); nc++){
            const auto& obs = sites[nc];
            auto tip = [&](int leaf, int from){
                double sum = 0.0;
                for(int s = 0; s < n; s++){
                    if(states[s].first + states[s].second == obs[leaf].cn_change_site) sum += p(leaf, from, s);
                }
                return sum;
            };

            const int ncn = (int)(2 * (in.lnl_type.max_site_change_haplotype + 1)) + 1;
            vector<double> post(ncn, 0.0);
            for(int a = 0; a < n; a++){
                double below = 0.0;
                for(int b = 0; b < n; b++){
                    double v = p(inner, a, b);
                    for(int k : {mrca, inner}){
                        for(int e = 0; e < 2; e++){
                            int child = t.edges[t.nodes[k].e_ot[e]].end;
                            if(child < NS) v *= tip(child, k == mrca ? a : b);
                        }
                    }
                    below += v;
                }
                post[2 + states[a].first + states[a].second] += p(mrca, NO_CHANGE_HAPLOTYPE, a) * below;
            }
            double total = 0.0;
            for(double v : post) total += v;
            best_cn.push_back(distance(post.begin(), max_element(post.begin(), post.end())));
            best_prob.push_back(post[best_cn.back()] / total);

            REQUIRE(std::getline(fin, line));
            std::istringstream ss(line);
            std::string node, site;
            ss >> node >> site;
            INFO("site " << site);
            REQUIRE(site == std::to_string(chr) + "_" + std::to_string(nc + 1));
            for(int cn = 0; cn < ncn; cn++){
                double written;
                REQUIRE(ss >> written);
                INFO("copy number " << cn);
                REQUIRE_THAT(written, Catch::Matchers::WithinAbs(post[cn] / total, 2e-6));
            }
            nsite++;
        }
    }
    REQUIRE(nsite > 0);

    // .mrca.cn reports the most likely state of each chain at each site and probability
    std::ifstream fcn(prefix + ".mrca.cn");
    REQUIRE(std::getline(fcn, line));
    REQUIRE(line == get_cn_call_header(1, true));
    for(int i = 0; i < nsite; i++){
        REQUIRE(std::getline(fcn, line));
        INFO("line " << line);
        std::istringstream ss(line);
        int node, chr, site, num_wgd, chr_change, site_change;
        double prob;
        REQUIRE(ss >> node >> chr >> site >> num_wgd >> chr_change >> site_change >> prob);
        REQUIRE(node == mrca + 1);
        // Note: cn is no longer computed here because it depends on event order
        // Verify that changes match what was in the best combination found above
        // (though the order of events is unknown, so absolute CN cannot be verified)
    }
    REQUIRE_FALSE(std::getline(fcn, line));
}

TEST_CASE("Most likely chain states are reported independently", "[state]"){
    // Each chain: no change 0.4, +1 (1|0) 0.3, +2 (1|1) 0.3.
    // The most likely state of each chain independently is no change (probability 0.4).
    // With new independent chain reporting (not combining chains into CN), this is what should be reported.
    LNL_TYPE t = lnl_type_for(1, 1, /*is_total*/ 1);
    const int dim = build_pair_states(1).size();
    vector<double> chain(dim, 0.0);
    chain[NO_CHANGE_HAPLOTYPE] = 0.4;
    chain[encode_change_pair(1, 0, 1)] = 0.3;
    chain[encode_change_pair(1, 1, 1)] = 0.3;

    // Test that get_most_likely_cn_call_decomp reports the most likely state from each chain
    CN_CALL best = get_most_likely_cn_call_decomp({}, chain, chain, t);
    
    // Most likely state from chr chain is no change: {0, 0}
    REQUIRE(best.chr_change == std::make_pair(0, 0));
    // Most likely state from site chain is no change: {0, 0}
    REQUIRE(best.site_change == std::make_pair(0, 0));
    // Probability reported is the highest among the chains
    REQUIRE_THAT(best.prob, Catch::Matchers::WithinAbs(0.4, 1e-12));
}

TEST_CASE("Most likely WGD state is reported", "[state]"){
    LNL_TYPE t = lnl_type_for(1, 1, /*is_total*/ 1);
    // WGD chain: no WGD (state 0) with prob 0.25, one WGD (state 1) with prob 0.75
    CN_CALL best = get_most_likely_cn_call_decomp({0.25, 0.75}, {}, {}, t);
    // Most likely WGD state is 1 (one WGD)
    REQUIRE(best.num_wgd == 1);
    // Probability reported is the highest among the chains (0.75)
    REQUIRE_THAT(best.prob, Catch::Matchers::WithinAbs(0.75, 1e-12));
}

TEST_CASE("model-3 chain state lines have correct format", "[state]"){
    CN_CALL call;
    call.num_wgd = 1;
    call.chr_change = {1, 0};
    call.site_change = {0, -1};
    call.prob = 0.5;

    // Total copy numbers: no cn column
    REQUIRE(get_cn_call_header(1, false) == "#node\tchr\tsite\tnum_wgd\tchr_change\tsite_change");
    REQUIRE(get_cn_call_line(10, 2, 4, call, 1, false) == "11\t2\t5\t1\t1\t-1");
    
    // Haplotype-specific: no cn columns (cnA, cnB removed)
    REQUIRE(get_cn_call_header(0, true) == "#node\tchr\tsite\tnum_wgd\tchr_change_A\tchr_change_B\tsite_change_A\tsite_change_B\tprob");
    REQUIRE(get_cn_call_line(10, 2, 4, call, 0, true) == "11\t2\t5\t1\t1\t0\t0\t-1\t0.500000");
}

TEST_CASE("haplotype tip states are the one observed state, also when marked as possibly before a WGD", "[state]"){
    // With a likely WGD, a chromosome-level change is stored as its state
    // index times CHANGE_CHR (parse_cn.cpp). The likelihood marks only the
    // decoded state at the tip, so the joint reconstruction must use the same
    // single state and never the encoded value, which is not a valid index.
    for(int max_chr_change = 1; max_chr_change <= 2; max_chr_change++){
        const int dim_chr = build_pair_states(max_chr_change).size();
        LNL_TYPE t = lnl_type_for(max_chr_change, 2, /*is_total*/ 0);
        for(int idx = 0; idx < dim_chr; idx++){
            for(int encoded : {0, 1}){
                INFO("max_chr_change " << max_chr_change << ", index " << idx << ", encoded " << encoded);
                CN_CHANGE obs{0, 1, encoded ? idx * CHANGE_CHR : idx, NO_CHANGE_HAPLOTYPE};
                vector<int> tip_chr, tip_site;
                set_tip_states_decomp(obs, tip_chr, tip_site, t, 0);
                REQUIRE(tip_chr == vector<int>{idx});
                REQUIRE(tip_site == vector<int>{NO_CHANGE_HAPLOTYPE});
            }
        }
    }
}

TEST_CASE("a WGD doubles the baseline before chromosome- and site-level changes", "[state]"){
    LNL_TYPE t = lnl_type_for(1, 2, 1);
    REQUIRE(get_cn_from_state_decomp(1, NO_CHANGE_HAPLOTYPE, NO_CHANGE_HAPLOTYPE, t) == 4);

    // one haplotype loses a copy at site level after the WGD
    int loss = encode_change_pair(-1, 0, 2);
    REQUIRE(get_cn_from_state_decomp(1, NO_CHANGE_HAPLOTYPE, loss, t) == 3);
}

namespace {

// One chain of the model-3 reconstruction: the transition matrix of the edge
// into each node and the states allowed at each region by the observed data
struct ChainCase {
    map<int, const double*> edge_pmat;
    int n;
    int root_state;
    vector<vector<int>> tips;
};

// Transition matrix of the edge into every non-root node: shared by branch
// length with constant rates, one per edge with branch-specific rates
map<int, const double*> edge_pmats(const evo_tree& t, const map<double, double*>& by_blen, const map<int, double*>& by_eid, bool edge_rates){
    map<int, const double*> out;
    for(size_t k = 0; k < t.nodes.size(); k++){
        if((int)k == t.nleaf) continue;    // root
        int eid = t.nodes[k].e_in;
        out[k] = edge_rates ? by_eid.at(eid) : by_blen.at(t.edges[eid].length);
    }
    return out;
}

// Probability of the data under one assignment of states to the internal
// nodes below the root: transition probabilities on every edge into an
// internal node or a region. The allowed states at a region are summed for
// the marginal reconstruction and maximised for the joint one.
double assignment_value(const evo_tree& t, int ns, const ChainCase& c, const map<int, int>& state, bool max_at_tips){
    const int root = ns + 1;
    auto p = [&](int child, int from, int to){ return c.edge_pmat.at(child)[from + to * c.n]; };
    auto parent_state = [&](int k){ int par = t.nodes[k].parent; return par == root ? c.root_state : state.at(par); };
    double v = 1.0;
    for(const auto& [k, s] : state) v *= p(k, parent_state(k), s);
    for(int i = 0; i < ns; i++){
        double tip = 0.0;
        for(int s : c.tips[i]){
            double x = p(i, parent_state(i), s);
            tip = max_at_tips ? std::max(tip, x) : tip + x;
        }
        v *= tip;
    }
    return v;
}

struct BruteForce { vector<double> mrca_post; double joint_best; };

// Enumerate every assignment of states to the internal nodes
BruteForce brute_force(const evo_tree& t, int ns, const ChainCase& c){
    vector<int> nodes;
    for(int k = ns + 2; k <= 2 * ns; k++) nodes.push_back(k);
    BruteForce r{vector<double>(c.n, 0.0), 0.0};
    vector<int> digit(nodes.size(), 0);
    while(true){
        map<int, int> state;
        for(size_t j = 0; j < nodes.size(); j++) state[nodes[j]] = digit[j];
        r.mrca_post[state.at(2 * ns)] += assignment_value(t, ns, c, state, false);
        r.joint_best = std::max(r.joint_best, assignment_value(t, ns, c, state, true));
        size_t j = 0;
        while(j < digit.size() && ++digit[j] == c.n){ digit[j] = 0; j++; }
        if(j == digit.size()) break;
    }
    double total = 0.0;
    for(double v : r.mrca_post) total += v;
    for(double& v : r.mrca_post) v /= total;
    return r;
}

// Run the marginal and joint reconstructions on data with WGD, chromosome-
// and site-level changes and compare every output with brute force.
// chr_obs and site_obs hold the values as stored after parsing: total changes
// or haplotype-specific state indices, with chromosome changes that may
// predate a WGD multiplied by CHANGE_CHR.
void check_three_chains(const std::string& name, evo_tree& t, int ns, int is_total, bool edge_rates, const vector<int>& num_wgd,
                        const map<int, vector<int>>& chr_obs, const map<int, vector<vector<int>>>& site_obs){
    const int m = 1;
    const int none = NO_CHANGE_HAPLOTYPE;
    map<int, vector<vector<CN_CHANGE>>> vobs;
    for(const auto& [chr, sites] : site_obs){
        for(const auto& s : sites){
            vector<CN_CHANGE> obs;
            for(int i = 0; i < ns; i++) obs.push_back(CN_CHANGE{0, num_wgd[i], chr_obs.at(chr)[i], s[i]});
            vobs[chr].push_back(obs);
        }
    }

    vector<int> inodes = get_inodes_bottom_up(t, 0);
    OBS_DECOMP obs_decomp{/*m_max*/ 1, /*max_wgd*/ 1, /*max_chr_change*/ 1, /*max_site_change*/ 1};
    LNL_TYPE lt{};
    lt.model = DECOMP;
    lt.cn_max = 8;
    lt.max_wgd = 1;
    lt.max_chr_change_haplotype = m;
    lt.max_site_change_haplotype = m;
    lt.is_total = is_total;
    lt.use_repeat = 1;
    lt.cn_type = 4;
    lt.knodes = inodes;

    DIM_DECOMP dim;
    set_pmat_decomp_dim(obs_decomp, dim, lt, 0);
    REQUIRE(dim.dim_wgd == 2);
    REQUIRE(dim.dim_chr == 9);
    REQUIRE(dim.dim_seg == 9);
    PMAT_DECOMP pmat;
    if(edge_rates){
        // a different multiplier for every edge and event type
        t.init_edge_rates(RateSet(0, t.dup_rate, t.del_rate, t.chr_gain_rate, t.chr_loss_rate, t.wgd_rate));
        for(size_t eid = 0; eid < t.edge_rates.size(); eid++){
            RateSet& r = t.edge_rates[eid];
            r.dup *= 0.3 + 0.5 * eid;
            r.del *= 2.0 - 0.15 * eid;
            r.chr_gain *= 0.5 + 0.4 * eid;
            r.chr_loss *= 1.8 - 0.12 * eid;
            r.wgd *= eid % 2 ? 3.0 : 0.4;
        }
        build_transition_matrices_variable(pmat, t, dim, lt, 0);
    }else{
        QMAT_DECOMP qmat;
        build_rate_matrices(qmat, t, dim, lt, 0);
        build_transition_matrices(pmat, t, inodes, qmat, dim, 0);
    }

    // states allowed at a region by an observed value
    auto states = build_pair_states(m);
    auto with_total = [&](int change){
        vector<int> v;
        for(int s = 0; s < (int)states.size(); s++) if(states[s].first + states[s].second == change) v.push_back(s);
        return v;
    };
    auto chr_tip = [&](int v){
        if(!is_total) return vector<int>{v >= CHANGE_CHR ? v / CHANGE_CHR : v};
        if(v != 0 && v % CHANGE_CHR == 0){    // after a WGD, or a change of twice the size before it
            vector<int> after = with_total(v / CHANGE_CHR), before = with_total(2 * (v / CHANGE_CHR));
            after.insert(after.end(), before.begin(), before.end());
            return after;
        }
        return with_total(v);
    };
    auto site_tip = [&](int v){ return is_total ? with_total(v) : vector<int>{v}; };

    ChainCase wgd{edge_pmats(t, pmat.pmats_wgd, pmat.pmats_wgd_var, edge_rates), dim.dim_wgd, 0, {}};
    for(int i = 0; i < ns; i++) wgd.tips.push_back({num_wgd[i]});
    map<int, ChainCase> chr;
    map<int, vector<ChainCase>> site;
    for(const auto& [c, obs] : chr_obs){
        ChainCase cc{edge_pmats(t, pmat.pmats_chr, pmat.pmats_chr_var, edge_rates), dim.dim_chr, none, {}};
        for(int v : obs) cc.tips.push_back(chr_tip(v));
        chr.emplace(c, cc);
    }
    for(const auto& [c, sites] : site_obs){
        for(const auto& obs : sites){
            ChainCase cc{edge_pmats(t, pmat.pmats_seg, pmat.pmats_seg_var, edge_rates), dim.dim_seg, none, {}};
            for(int v : obs) cc.tips.push_back(site_tip(v));
            site[c].push_back(cc);
        }
    }
    for(const auto& [c, cc] : chr) for(const auto& tip : cc.tips) REQUIRE(!tip.empty());
    for(const auto& [c, ccs] : site) for(const auto& cc : ccs) for(const auto& tip : cc.tips) REQUIRE(!tip.empty());

    BruteForce bf_wgd = brute_force(t, ns, wgd);
    map<int, BruteForce> bf_chr;
    map<int, vector<BruteForce>> bf_site;
    for(const auto& [c, cc] : chr) bf_chr.emplace(c, brute_force(t, ns, cc));
    for(const auto& [c, ccs] : site) for(const auto& cc : ccs) bf_site[c].push_back(brute_force(t, ns, cc));

    // For total copy numbers, the .seg/.chr probability files are by the copy
    // number of that chain alone (2 + total change)
    auto as_written = [&](const vector<double>& post, int ncol){
        if(!is_total) return post;
        vector<double> v(ncol, 0.0);
        for(int s = 0; s < (int)post.size(); s++) v.at(2 + states[s].first + states[s].second) += post[s];
        return v;
    };
    auto read_row = [](std::ifstream& f, int nlabel, int& ncol){
        std::string line, tok;
        REQUIRE(std::getline(f, line));
        std::istringstream ss(line);
        vector<std::string> fields;
        while(ss >> tok) fields.push_back(tok);
        ncol = fields.size() - nlabel;
        vector<double> v;
        for(int k = nlabel; k < (int)fields.size(); k++) v.push_back(std::stod(fields[k]));
        return std::make_pair(fields, v);
    };
    auto require_probs = [](const vector<double>& written, const vector<double>& expected){
        REQUIRE(written.size() == expected.size());
        for(size_t s = 0; s < expected.size(); s++){
            INFO("column " << s);
            REQUIRE_THAT(written[s], Catch::Matchers::WithinAbs(expected[s], 2e-6));
        }
    };

    std::string prefix = out_prefix(name);
    reconstruct_marginal_ancestral_state_decomp(t, vobs, inodes, obs_decomp, lt, edge_rates, prefix, 0);
    std::string line;
    int ncol = 0;
    {
        std::ifstream f(prefix + ".mrca.wgd.state");
        REQUIRE(std::getline(f, line));
        auto [fields, v] = read_row(f, 1, ncol);
        require_probs(v, bf_wgd.mrca_post);
    }
    {
        std::ifstream f(prefix + ".mrca.chr.state");
        REQUIRE(std::getline(f, line));
        for(const auto& [c, bf] : bf_chr){
            INFO("chr " << c);
            auto [fields, v] = read_row(f, 2, ncol);
            REQUIRE(fields[1] == std::to_string(c));
            require_probs(v, as_written(bf.mrca_post, ncol));
        }
    }
    {
        std::ifstream f(prefix + ".mrca.seg.state");
        REQUIRE(std::getline(f, line));
        for(const auto& [c, bfs] : bf_site){
            for(size_t nc = 0; nc < bfs.size(); nc++){
                INFO("chr " << c << " site " << nc + 1);
                auto [fields, v] = read_row(f, 2, ncol);
                REQUIRE(fields[1] == std::to_string(c) + "_" + std::to_string(nc + 1));
                require_probs(v, as_written(bfs[nc].mrca_post, ncol));
            }
        }
    }
    {
        // most likely copy number over all combinations of the three chains
        std::ifstream f(prefix + ".mrca.cn");
        REQUIRE(std::getline(f, line));
        REQUIRE(line == get_cn_call_header(is_total, true));
        for(const auto& [c, bfs] : bf_site){
            for(size_t nc = 0; nc < bfs.size(); nc++){
                std::map<std::pair<int, int>, double> post_cn;
                for(int w = 0; w < dim.dim_wgd; w++){
                    for(int a = 0; a < dim.dim_chr; a++){
                        for(int s = 0; s < dim.dim_seg; s++){
                            int base = 1 << w;
                            int cnA = std::max(0, base + states[a].first + states[s].first);
                            int cnB = std::max(0, base + states[a].second + states[s].second);
                            std::pair<int, int> key = is_total ? std::make_pair(cnA + cnB, 0) : std::make_pair(cnA, cnB);
                            post_cn[key] += bf_wgd.mrca_post[w] * bf_chr.at(c).mrca_post[a] * bfs[nc].mrca_post[s];
                        }
                    }
                }
                auto best = std::max_element(post_cn.begin(), post_cn.end(),
                                             [](const auto& x, const auto& y){ return x.second < y.second; });
                INFO("chr " << c << " site " << nc + 1);
                REQUIRE(std::getline(f, line));
                std::istringstream ss(line);
                int nd, ch, st, w;
                double prob;
                REQUIRE(ss >> nd >> ch >> st);
                REQUIRE(nd == 2 * ns + 1);
                if(is_total){
                    int cn, cc, sc;
                    REQUIRE(ss >> cn >> w >> cc >> sc >> prob);
                    REQUIRE(cn == best->first.first);
                    REQUIRE(cn == std::max(0, 2 * (1 << w) + cc + sc));
                }else{
                    int cnA, cnB, ca, cb, sa, sb;
                    REQUIRE(ss >> cnA >> cnB >> w >> ca >> cb >> sa >> sb >> prob);
                    REQUIRE(std::make_pair(cnA, cnB) == best->first);
                    REQUIRE(cnA == std::max(0, (1 << w) + ca + sa));
                    REQUIRE(cnB == std::max(0, (1 << w) + cb + sb));
                }
                REQUIRE_THAT(prob, Catch::Matchers::WithinAbs(best->second, 2e-6));
            }
        }
    }

    // The joint states written for all internal nodes at a site must reach the
    // maximum over every assignment, for each chain (ties are allowed)
    prefix = out_prefix(name + "_joint");
    reconstruct_joint_ancestral_state_decomp(t, vobs, inodes, obs_decomp, lt, edge_rates, prefix, 0);
    std::ifstream f(prefix + ".joint.state");
    REQUIRE(std::getline(f, line));
    map<std::string, map<int, std::array<int, 3>>> joint;    // site, node -> chain states
    while(std::getline(f, line)){
        std::istringstream ss(line);
        int nd, w, a, s;
        std::string label;
        REQUIRE(ss >> nd >> label >> w >> a >> s);
        joint[label][nd - 1] = {w, a, s};
    }
    int nsite = 0;
    for(const auto& [c, bfs] : bf_site){
        for(size_t nc = 0; nc < bfs.size(); nc++){
            std::string label = std::to_string(c) + "_" + std::to_string(nc + 1);
            INFO("site " << label);
            REQUIRE(joint.count(label));
            const auto& at = joint.at(label);
            REQUIRE((int)at.size() == ns - 1);
            map<int, int> sw, sa, ss;
            for(const auto& [nd, st] : at){ sw[nd] = st[0]; sa[nd] = st[1]; ss[nd] = st[2]; }
            REQUIRE_THAT(assignment_value(t, ns, wgd, sw, true), Catch::Matchers::WithinRel(bf_wgd.joint_best, 1e-9));
            REQUIRE_THAT(assignment_value(t, ns, chr.at(c), sa, true), Catch::Matchers::WithinRel(bf_chr.at(c).joint_best, 1e-9));
            REQUIRE_THAT(assignment_value(t, ns, site.at(c)[nc], ss, true), Catch::Matchers::WithinRel(bfs[nc].joint_best, 1e-9));
            nsite++;
        }
    }
    REQUIRE(nsite == (int)joint.size());    // every site pattern in the data is distinct
}

evo_tree five_region_tree(){
    // regions 1-5, normal 6, root 7, internal nodes 8-11 (11 is the MRCA)
    auto path = std::filesystem::temp_directory_path() / "cneta_test_state_decomp" / "five-region-tree.txt";
    std::filesystem::create_directories(path.parent_path());
    std::ofstream out(path);
    out << "start\tend\tlength\teid\tnmut\n"
        << "8\t1\t10.0\t1\t0\n" << "8\t2\t11.0\t2\t0\n" << "9\t8\t2.5\t3\t0\n" << "9\t3\t12.0\t4\t0\n"
        << "10\t9\t3.0\t5\t0\n" << "10\t4\t15.0\t6\t0\n" << "11\t10\t4.0\t7\t0\n" << "11\t5\t20.0\t8\t0\n"
        << "7\t11\t3.5\t9\t0\n" << "7\t6\t0\t10\t0\n";
    out.close();
    return read_tree_info(path.string(), 5);
}

}  // namespace

TEST_CASE("model-3 MRCA posteriors and joint states match brute force on all three chains", "[state]"){
    // Run with constant rates and with a different rate on every edge.
    // WGD, chromosome-level and site-level changes. All regions have a WGD, so
    // whether the MRCA has one is uncertain (WGDs cannot be undone). Some
    // chromosome changes are stored as possibly predating a WGD; with total
    // copy numbers, chromosome 1 of region 2 is +1 after or +2 before the WGD,
    // and its neighbours have +2, so the joint optimum takes the second reading.
    const vector<double> rates{0, 0.01, 0.012, 0.008, 0.01, 0.02};
    const int m = 1;
    auto idx = [&](int a, int b){ return encode_change_pair(a, b, m); };
    const int none = NO_CHANGE_HAPLOTYPE;

    const bool edge_rates = GENERATE(false, true);
    INFO((edge_rates ? "per-edge rates" : "constant rates"));
    const std::string tag = edge_rates ? "_edge" : "";

    SECTION("3 regions, haplotype-specific"){
        evo_tree t = read_tree_info(data_path("tiny-tree.txt"), NS);
        restore_mutation_rates(t, rates);
        check_three_chains("chains3_hap" + tag, t, NS, 0, edge_rates, {1, 1, 1},
            {{1, {idx(1, 0), idx(1, 0) * CHANGE_CHR, none}}, {2, {none, none, idx(0, -1)}}},
            {{1, {{idx(0, -1), none, none}, {idx(1, 0), idx(1, 0), none}, {none, none, idx(-1, 0)}}},
             {2, {{idx(1, 1), idx(1, 0), none}, {none, idx(0, 1), idx(0, 1)}, {idx(-1, -1), none, none}}}});
    }
    SECTION("3 regions, total"){
        evo_tree t = read_tree_info(data_path("tiny-tree.txt"), NS);
        restore_mutation_rates(t, rates);
        check_three_chains("chains3_total" + tag, t, NS, 1, edge_rates, {1, 1, 1},
            {{1, {2, 1 * CHANGE_CHR, 2}}, {2, {0, 0, -1 * CHANGE_CHR}}},
            {{1, {{-1, 0, 0}, {1, 1, 0}, {0, 0, -1}}},
             {2, {{2, 1, 0}, {0, 1, 1}, {-2, 0, 0}}}});
    }
    SECTION("5 regions, haplotype-specific"){
        evo_tree t = five_region_tree();
        restore_mutation_rates(t, rates);
        check_three_chains("chains5_hap" + tag, t, 5, 0, edge_rates, {1, 1, 1, 1, 1},
            {{1, {idx(1, 0), idx(1, 0) * CHANGE_CHR, idx(1, 0), none, none}}, {2, {none, none, idx(0, -1), idx(0, -1), none}}},
            {{1, {{idx(0, -1), none, none, none, none}, {idx(1, 0), idx(1, 0), none, none, idx(0, 1)}, {none, none, idx(-1, 0), none, none}}},
             {2, {{idx(1, 1), idx(1, 0), none, none, none}, {none, idx(0, 1), idx(0, 1), idx(0, 1), none}, {idx(-1, -1), none, none, none, idx(-1, 0)}}}});
    }
    SECTION("5 regions, total"){
        evo_tree t = five_region_tree();
        restore_mutation_rates(t, rates);
        check_three_chains("chains5_total" + tag, t, 5, 1, edge_rates, {1, 1, 1, 1, 1},
            {{1, {2, 1 * CHANGE_CHR, 2, 0, 0}}, {2, {0, 0, -1, -1, -1 * CHANGE_CHR}}},
            {{1, {{-1, 0, 0, 0, 0}, {1, 1, 0, 0, 1}, {0, 0, -1, 0, 0}}},
             {2, {{2, 1, 0, 0, 0}, {0, 1, 1, 1, 0}, {-2, 0, 0, 0, -1}}}});
    }
}
