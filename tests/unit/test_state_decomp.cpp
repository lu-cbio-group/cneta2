// Unit tests for turning independent-chain (model 3) states back into copy
// numbers in ancestral state reconstruction (state.cpp).
//
// Chromosome- and site-level states are indices into the list built by
// build_pair_states() (common.hpp), which is ordered by total copy number and
// then by haplotype A. Decoding must use the same order, or reconstructed copy
// numbers are wrong even when the reconstructed states are right.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

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
    reconstruct_marginal_ancestral_state_decomp(t, in.vobs_change, in.inodes, in.obs_decomp, in.lnl_type, prefix, 0);

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

    // .mrca.cn reports the most likely copy number at each site and its probability
    std::ifstream fcn(prefix + ".mrca.cn");
    REQUIRE(std::getline(fcn, line));
    REQUIRE(line == get_cn_call_header(1, true));
    for(int i = 0; i < nsite; i++){
        REQUIRE(std::getline(fcn, line));
        INFO("line " << line);
        std::istringstream ss(line);
        int node, chr, site, cn, num_wgd, chr_change, site_change;
        double prob;
        REQUIRE(ss >> node >> chr >> site >> cn >> num_wgd >> chr_change >> site_change >> prob);
        REQUIRE(node == mrca + 1);
        REQUIRE(cn == best_cn[i]);
        REQUIRE(cn == 2 + chr_change + site_change);
        REQUIRE_THAT(prob, Catch::Matchers::WithinAbs(best_prob[i], 2e-6));
    }
    REQUIRE_FALSE(std::getline(fcn, line));
}

TEST_CASE("MRCA copy number posterior combines the independent chains", "[state]"){
    // Chromosome and site chains each: no change 0.4, +1 (1|0) 0.3, +2 (1|1) 0.3.
    // The most likely state of each chain is no change (copy number 2, 0.16),
    // but the most likely copy number is 4: 0.3 * 0.3 + 2 * 0.4 * 0.3 = 0.33.
    LNL_TYPE t = lnl_type_for(1, 1, /*is_total*/ 1);
    const int dim = build_pair_states(1).size();
    vector<double> chain(dim, 0.0);
    chain[NO_CHANGE_HAPLOTYPE] = 0.4;
    chain[encode_change_pair(1, 0, 1)] = 0.3;
    chain[encode_change_pair(1, 1, 1)] = 0.3;

    auto post_cn = get_allele_cn_posterior_decomp({}, chain, chain, t);
    std::map<int, double> by_total;
    for(const auto& [cn, p] : post_cn) by_total[cn.first + cn.second] += p;
    const std::map<int, double> expected{{2, 0.16}, {3, 0.24}, {4, 0.33}, {5, 0.18}, {6, 0.09}};
    REQUIRE(by_total.size() == expected.size());
    for(const auto& [cn, p] : expected){
        INFO("copy number " << cn);
        REQUIRE_THAT(by_total[cn], Catch::Matchers::WithinAbs(p, 1e-12));
    }

    CN_CALL best = get_most_likely_cn_call_decomp({}, chain, chain, t);
    REQUIRE(best.cn.first + best.cn.second == 4);
    REQUIRE_THAT(best.prob, Catch::Matchers::WithinAbs(0.33, 1e-12));
    // the chain changes reported are one of the two most likely ways to get 4 (0.12 each)
    REQUIRE(best.chr_change.first + best.chr_change.second + best.site_change.first + best.site_change.second == 2);
    REQUIRE(best.chr_change != best.site_change);

    // unnormalised chain posteriors give the same result
    vector<double> scaled = chain;
    for(double& p : scaled) p *= 7.0;
    REQUIRE_THAT(get_most_likely_cn_call_decomp({}, scaled, chain, t).prob, Catch::Matchers::WithinAbs(0.33, 1e-12));

    // haplotype-specific: (2,2) arises from chr 1|1 + site 0|0 or chr 0|0 + site 1|1
    REQUIRE_THAT(post_cn.at({2, 2}), Catch::Matchers::WithinAbs(0.24, 1e-12));
    REQUIRE_THAT(post_cn.at({3, 2}), Catch::Matchers::WithinAbs(0.18, 1e-12));
}

TEST_CASE("a WGD in the MRCA posterior doubles the baseline copy number", "[state]"){
    LNL_TYPE t = lnl_type_for(1, 1, /*is_total*/ 1);
    auto post_cn = get_allele_cn_posterior_decomp({0.25, 0.75}, {}, {}, t);
    REQUIRE(post_cn.size() == 2);
    REQUIRE_THAT(post_cn.at({1, 1}), Catch::Matchers::WithinAbs(0.25, 1e-12));
    REQUIRE_THAT(post_cn.at({2, 2}), Catch::Matchers::WithinAbs(0.75, 1e-12));
    CN_CALL best = get_most_likely_cn_call_decomp({0.25, 0.75}, {}, {}, t);
    REQUIRE(best.num_wgd == 1);
    REQUIRE(best.cn == std::make_pair(2, 2));
    REQUIRE_THAT(best.prob, Catch::Matchers::WithinAbs(0.75, 1e-12));
}

TEST_CASE("model-3 copy number lines have one value per column", "[state]"){
    CN_CALL call;
    call.num_wgd = 1;
    call.chr_change = {1, 0};
    call.site_change = {0, -1};
    call.cn = {3, 1};
    call.prob = 0.5;

    REQUIRE(get_cn_call_header(1, false) == "#node\tchr\tsite\tcn\tnum_wgd\tchr_change\tsite_change");
    REQUIRE(get_cn_call_line(10, 2, 4, call, 1, false) == "11\t2\t5\t4\t1\t1\t-1");
    REQUIRE(get_cn_call_header(0, true) == "#node\tchr\tsite\tcnA\tcnB\tnum_wgd\tchr_change_A\tchr_change_B\tsite_change_A\tsite_change_B\tprob");
    REQUIRE(get_cn_call_line(10, 2, 4, call, 0, true) == "11\t2\t5\t3\t1\t1\t1\t0\t0\t-1\t0.500000");
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

TEST_CASE("model-3 MRCA posteriors and joint states match brute force on all three chains", "[state]"){
    // Haplotype-specific changes on the 3-region tree, with WGD, chromosome-level
    // and site-level changes. Each chain is enumerated over every state of the
    // MRCA and of its internal child: the sum gives the marginal posterior at the
    // MRCA (including the trunk edge from the normal root), the maximum gives the
    // joint reconstruction. Region 2's chromosome-1 gain is written as possibly
    // before a WGD (state index times CHANGE_CHR).
    const int m = 1;
    auto idx = [&](int a, int b){ return encode_change_pair(a, b, m); };
    const int none = NO_CHANGE_HAPLOTYPE;

    const vector<int> num_wgd{1, 1, 1};    // all regions doubled, so the MRCA may or may not be
    const map<int, vector<int>> chr_tip{{1, {idx(1, 0), idx(1, 0), none}}, {2, {none, none, idx(0, -1)}}};
    const map<int, vector<vector<int>>> site_tip{
        {1, {{idx(0, -1), none, none}, {idx(1, 0), idx(1, 0), none}, {none, none, idx(-1, 0)}}},
        {2, {{idx(1, 1), idx(1, 0), none}, {none, idx(0, 1), idx(0, 1)}, {idx(-1, -1), none, none}}}};

    map<int, vector<vector<CN_CHANGE>>> vobs;
    for(const auto& [chr, sites] : site_tip){
        for(const auto& s : sites){
            vector<CN_CHANGE> obs;
            for(int i = 0; i < NS; i++){
                int c = chr_tip.at(chr)[i];
                if(chr == 1 && i == 1) c *= CHANGE_CHR;
                obs.push_back(CN_CHANGE{0, num_wgd[i], c, s[i]});
            }
            vobs[chr].push_back(obs);
        }
    }

    evo_tree t = read_tree_info(data_path("tiny-tree.txt"), NS);
    restore_mutation_rates(t, vector<double>{0, 0.01, 0.012, 0.008, 0.01, 0.02});
    vector<int> inodes = get_inodes_bottom_up(t, 0);
    OBS_DECOMP obs_decomp{/*m_max*/ 1, /*max_wgd*/ 1, /*max_chr_change*/ 1, /*max_site_change*/ 1};
    LNL_TYPE lt{};
    lt.model = DECOMP;
    lt.cn_max = 8;
    lt.max_wgd = 1;
    lt.max_chr_change_haplotype = m;
    lt.max_site_change_haplotype = m;
    lt.is_total = 0;
    lt.use_repeat = 1;
    lt.cn_type = 4;
    lt.knodes = inodes;

    DIM_DECOMP dim;
    set_pmat_decomp_dim(obs_decomp, dim, lt, 0);
    REQUIRE(dim.dim_wgd == 2);
    REQUIRE(dim.dim_chr == 9);
    REQUIRE(dim.dim_seg == 9);
    QMAT_DECOMP qmat;
    build_rate_matrices(qmat, t, dim, lt, 0);
    PMAT_DECOMP pmat;
    build_transition_matrices(pmat, t, inodes, qmat, dim, 0);

    const int mrca = 2 * NS;
    const int inner = NS + 2;
    REQUIRE(t.edges[t.nodes[inner].e_in].start == mrca);

    // posterior (sum) and best states (max) of one chain at the MRCA and its child
    struct ChainResult { vector<double> post; int best_mrca; int best_inner; };
    auto solve = [&](const map<double, double*>& pmats, int n, int root_state, const vector<int>& tips){
        auto p = [&](int child, int from, int to){ return pmats.at(t.edges[t.nodes[child].e_in].length)[from + to * n]; };
        ChainResult r{vector<double>(n, 0.0), -1, -1};
        double best = -1, second = -1;
        for(int a = 0; a < n; a++){
            for(int b = 0; b < n; b++){
                double v = p(mrca, root_state, a) * p(inner, a, b);
                for(int k : {mrca, inner}){
                    for(int e = 0; e < 2; e++){
                        int child = t.edges[t.nodes[k].e_ot[e]].end;
                        if(child < NS) v *= p(child, k == mrca ? a : b, tips[child]);
                    }
                }
                r.post[a] += v;
                if(v > best){ second = best; best = v; r.best_mrca = a; r.best_inner = b; }
                else if(v > second){ second = v; }
            }
        }
        REQUIRE(best > second * (1 + 1e-6));    // no tie in the joint optimum
        double total = 0.0;
        for(double v : r.post) total += v;
        for(double& v : r.post) v /= total;
        return r;
    };

    ChainResult wgd = solve(pmat.pmats_wgd, dim.dim_wgd, 0, num_wgd);
    map<int, ChainResult> chr;
    map<int, vector<ChainResult>> site;
    for(const auto& [c, tips] : chr_tip) chr.emplace(c, solve(pmat.pmats_chr, dim.dim_chr, none, tips));
    for(const auto& [c, sites] : site_tip){
        for(const auto& tips : sites) site[c].push_back(solve(pmat.pmats_seg, dim.dim_seg, none, tips));
    }

    auto read_probs = [](std::istringstream& ss, int n){
        vector<double> v(n);
        for(double& x : v) REQUIRE(ss >> x);
        return v;
    };
    auto require_probs = [](const vector<double>& written, const vector<double>& expected){
        for(size_t s = 0; s < expected.size(); s++){
            INFO("state " << s);
            REQUIRE_THAT(written[s], Catch::Matchers::WithinAbs(expected[s], 2e-6));
        }
    };

    std::string prefix = out_prefix("three_chains");
    reconstruct_marginal_ancestral_state_decomp(t, vobs, inodes, obs_decomp, lt, prefix, 0);
    std::string line, node, label;

    {
        std::ifstream f(prefix + ".mrca.wgd.state");
        REQUIRE(std::getline(f, line));
        REQUIRE(std::getline(f, line));
        std::istringstream ss(line);
        ss >> node;
        require_probs(read_probs(ss, dim.dim_wgd), wgd.post);
    }
    {
        std::ifstream f(prefix + ".mrca.chr.state");
        REQUIRE(std::getline(f, line));
        for(const auto& [c, r] : chr){
            INFO("chr " << c);
            REQUIRE(std::getline(f, line));
            std::istringstream ss(line);
            ss >> node >> label;
            REQUIRE(label == std::to_string(c));
            require_probs(read_probs(ss, dim.dim_chr), r.post);
        }
    }
    {
        std::ifstream f(prefix + ".mrca.seg.state");
        REQUIRE(std::getline(f, line));
        for(const auto& [c, rs] : site){
            for(size_t nc = 0; nc < rs.size(); nc++){
                INFO("chr " << c << " site " << nc + 1);
                REQUIRE(std::getline(f, line));
                std::istringstream ss(line);
                ss >> node >> label;
                REQUIRE(label == std::to_string(c) + "_" + std::to_string(nc + 1));
                require_probs(read_probs(ss, dim.dim_seg), rs[nc].post);
            }
        }
    }
    {
        // most likely (cnA, cnB) over all combinations of the three chains
        auto states = build_pair_states(m);
        std::ifstream f(prefix + ".mrca.cn");
        REQUIRE(std::getline(f, line));
        REQUIRE(line == get_cn_call_header(0, true));
        for(const auto& [c, rs] : site){
            for(size_t nc = 0; nc < rs.size(); nc++){
                std::map<std::pair<int, int>, double> post_cn;
                for(int w = 0; w < dim.dim_wgd; w++){
                    for(int a = 0; a < dim.dim_chr; a++){
                        for(int s = 0; s < dim.dim_seg; s++){
                            int base = 1 << w;
                            std::pair<int, int> cn{std::max(0, base + states[a].first + states[s].first),
                                                   std::max(0, base + states[a].second + states[s].second)};
                            post_cn[cn] += wgd.post[w] * chr.at(c).post[a] * rs[nc].post[s];
                        }
                    }
                }
                auto best = std::max_element(post_cn.begin(), post_cn.end(),
                                             [](const auto& x, const auto& y){ return x.second < y.second; });
                INFO("chr " << c << " site " << nc + 1);
                REQUIRE(std::getline(f, line));
                std::istringstream ss(line);
                int nd, ch, st, cnA, cnB, w, ca, cb, sa, sb;
                double prob;
                REQUIRE(ss >> nd >> ch >> st >> cnA >> cnB >> w >> ca >> cb >> sa >> sb >> prob);
                REQUIRE(std::make_pair(cnA, cnB) == best->first);
                REQUIRE_THAT(prob, Catch::Matchers::WithinAbs(best->second, 2e-6));
                REQUIRE(cnA == std::max(0, (1 << w) + ca + sa));
                REQUIRE(cnB == std::max(0, (1 << w) + cb + sb));
            }
        }
    }

    prefix = out_prefix("three_chains_joint");
    reconstruct_joint_ancestral_state_decomp(t, vobs, inodes, obs_decomp, lt, prefix, 0);
    std::ifstream f(prefix + ".joint.state");
    REQUIRE(std::getline(f, line));
    int nrow = 0;
    while(std::getline(f, line)){
        std::istringstream ss(line);
        int nd, w, a, s;
        REQUIRE(ss >> nd >> label >> w >> a >> s);
        int c = std::stoi(label.substr(0, label.find('_')));
        int nc = std::stoi(label.substr(label.find('_') + 1)) - 1;
        INFO("node " << nd << " site " << label);
        REQUIRE((nd == mrca + 1 || nd == inner + 1));
        bool at_mrca = nd == mrca + 1;
        REQUIRE(w == (at_mrca ? wgd.best_mrca : wgd.best_inner));
        REQUIRE(a == (at_mrca ? chr.at(c).best_mrca : chr.at(c).best_inner));
        REQUIRE(s == (at_mrca ? site.at(c)[nc].best_mrca : site.at(c)[nc].best_inner));
        nrow++;
    }
    REQUIRE(nrow == 2 * 6);    // two internal nodes, six distinct site patterns
}
