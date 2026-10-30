// Unit tests for ancestral state reconstruction under the independent-chain
// model (state.cpp), with and without branch-specific rates.
//
// The reconstruction recomputes the tree likelihood on the way, so its
// log-likelihood must agree with the likelihood code used during fitting.
// That agreement is what ties the per-edge transition matrices used here to
// the ones the optimiser sees.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <filesystem>
#include <fstream>
#include <iterator>
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
// number, segment-level changes only (the fixture has no chromosome-level
// events or WGD).
struct AsrInput {
    evo_tree tree;
    map<int, vector<vector<CN_CHANGE>>> vobs_change;
    OBS_DECOMP obs_decomp;
    LNL_TYPE lnl_type;
    vector<int> inodes;
};

AsrInput load_fixture(){
    AsrInput in;

    const int cn_max = 4, max_wgd = 1, max_chr_change_haplotype = 1, max_site_change_haplotype = 2;
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

RateSet global_rates(const evo_tree& tree){
    return RateSet(0, tree.dup_rate, tree.del_rate, tree.chr_gain_rate, tree.chr_loss_rate, tree.wgd_rate);
}

// A different rate on every edge, so a matrix taken from the wrong edge shows up.
void set_distinct_edge_rates(evo_tree& tree){
    tree.init_edge_rates(global_rates(tree));
    for(size_t eid = 0; eid < tree.edge_rates.size(); eid++){
        tree.edge_rates[eid] = tree.edge_rates[eid] * (1.0 + 0.7 * eid);
    }
}

std::string read_file(const std::string& path){
    std::ifstream in(path);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

// Removes earlier outputs with this prefix, so a run that writes nothing cannot
// pass by reading files left over from a previous run.
std::string out_prefix(const std::string& name){
    auto dir = std::filesystem::temp_directory_path() / "cneta_test_state";
    std::filesystem::create_directories(dir);
    for(const auto& entry : std::filesystem::directory_iterator(dir)){
        if(entry.path().filename().string().rfind(name + ".", 0) == 0){
            std::filesystem::remove(entry.path());
        }
    }
    return (dir / name).string();
}

}  // namespace

TEST_CASE("marginal reconstruction with per-edge rates reproduces the variable-rate likelihood", "[state]"){
    AsrInput in = load_fixture();
    set_distinct_edge_rates(in.tree);

    double lnl_asr = reconstruct_marginal_ancestral_state_decomp(
        in.tree, in.vobs_change, in.inodes, in.obs_decomp, in.lnl_type, true, out_prefix("var"), 0);

    DIM_DECOMP dim;
    set_pmat_decomp_dim(in.obs_decomp, dim, in.lnl_type, 0);
    PMAT_DECOMP pmat;
    build_transition_matrices_variable(pmat, in.tree, dim, in.lnl_type, 0);
    double lnl_fit = get_likelihood_chr_change_variable(in.tree, in.vobs_change, pmat, dim, in.lnl_type, 0);

    REQUIRE_THAT(lnl_asr, Catch::Matchers::WithinRel(lnl_fit, 1e-12));

    SECTION("and differs from the global-rate likelihood"){
        double lnl_global = reconstruct_marginal_ancestral_state_decomp(
            in.tree, in.vobs_change, in.inodes, in.obs_decomp, in.lnl_type, false, out_prefix("global"), 0);
        REQUIRE_FALSE(lnl_asr == lnl_global);
    }
}

TEST_CASE("reconstruction with the global rate on every edge matches the constant-rate path", "[state]"){
    AsrInput in = load_fixture();
    in.tree.init_edge_rates(global_rates(in.tree));

    std::string p_const = out_prefix("const");
    std::string p_edge = out_prefix("edge");

    double lnl_const = reconstruct_marginal_ancestral_state_decomp(
        in.tree, in.vobs_change, in.inodes, in.obs_decomp, in.lnl_type, false, p_const, 0);
    double lnl_edge = reconstruct_marginal_ancestral_state_decomp(
        in.tree, in.vobs_change, in.inodes, in.obs_decomp, in.lnl_type, true, p_edge, 0);
    REQUIRE(lnl_const == lnl_edge);

    reconstruct_joint_ancestral_state_decomp(
        in.tree, in.vobs_change, in.inodes, in.obs_decomp, in.lnl_type, false, p_const, 0);
    reconstruct_joint_ancestral_state_decomp(
        in.tree, in.vobs_change, in.inodes, in.obs_decomp, in.lnl_type, true, p_edge, 0);

    for(const std::string suffix : {".mrca.seg.state", ".mrca.cn", ".joint.state", ".joint.cn"}){
        INFO("output " << suffix);
        std::string a = read_file(p_const + suffix);
        REQUIRE_FALSE(a.empty());
        REQUIRE(a == read_file(p_edge + suffix));
    }
}

TEST_CASE("joint reconstruction with per-edge rates finds the most probable assignment", "[state]"){
    // On the 3-region tree there are only two internal nodes below the root, so
    // every assignment of segment-level states can be scored directly with the
    // per-edge matrices. The dynamic programme must reach the same maximum, and
    // the assignment it traces back must attain it. A matrix taken from the
    // wrong edge would give a different maximum.
    AsrInput in = load_fixture();
    set_distinct_edge_rates(in.tree);

    DIM_DECOMP dim;
    set_pmat_decomp_dim(in.obs_decomp, dim, in.lnl_type, 0);
    REQUIRE(dim.dim_seg > 0);
    REQUIRE(dim.dim_chr == 0);
    REQUIRE(dim.dim_wgd == 0);
    PMAT_DECOMP pmat;
    build_transition_matrices_variable(pmat, in.tree, dim, in.lnl_type, 0);

    const evo_tree& t = in.tree;
    const int n = dim.dim_seg;
    const int root = t.nleaf;            // fixed in the normal state
    const int mrca = 2 * NS;
    const int inner = NS + 2;            // the other internal node, child of the MRCA
    REQUIRE(t.edges[t.nodes[inner].e_in].start == mrca);

    auto p = [&](int child, int from, int to){
        return pmat.pmats_seg_var.at(t.nodes[child].e_in)[from + to * n];
    };

    int nsite = 0;
    for(const auto& [chr, sites] : in.vobs_change){
        for(const auto& obs : sites){
            nsite++;
            // best transition into each tip over the states its observation allows
            auto tip = [&](int leaf, int from){
                vector<int> chr_states, site_states;
                set_tip_states_decomp(obs[leaf], chr_states, site_states, in.lnl_type, 0);
                double best = 0.0;
                for(int s : site_states) best = max(best, p(leaf, from, s));
                return best;
            };
            auto score = [&](int s_mrca, int s_inner){
                double v = p(mrca, NO_CHANGE_HAPLOTYPE, s_mrca) * p(inner, s_mrca, s_inner);
                for(int k : {mrca, inner}){
                    for(int e = 0; e < 2; e++){
                        int child = t.edges[t.nodes[k].e_ot[e]].end;
                        if(child < NS) v *= tip(child, k == mrca ? s_mrca : s_inner);
                    }
                }
                return v;
            };

            double brute = 0.0;
            for(int a = 0; a < n; a++){
                for(int b = 0; b < n; b++) brute = max(brute, score(a, b));
            }

            LNL_TABLE lnl;
            STATE_TABLE state;
            initialize_lnl_state_table(lnl, state, 2 * t.nleaf - 1, dim);
            initialize_asr_table_decomp(obs, t, pmat, dim, lnl, state, in.lnl_type, true, 0);
            get_ancestral_states_site_decomp(lnl, state, t, in.inodes, pmat, dim, true, 0);

            INFO("chr " << chr << ", site " << nsite);
            REQUIRE_THAT(lnl.lnl_table_seg[mrca][NO_CHANGE_HAPLOTYPE], Catch::Matchers::WithinRel(brute, 1e-12));
            int s_mrca = state.state_table_seg[mrca][NO_CHANGE_HAPLOTYPE];
            int s_inner = state.state_table_seg[inner][s_mrca];
            REQUIRE_THAT(score(s_mrca, s_inner), Catch::Matchers::WithinRel(brute, 1e-12));
        }
    }
    REQUIRE(nsite > 0);
}

TEST_CASE("joint reconstruction uses the per-edge rates", "[state]"){
    AsrInput in = load_fixture();
    std::string p_global = out_prefix("joint_global");
    std::string p_var = out_prefix("joint_var");

    reconstruct_joint_ancestral_state_decomp(
        in.tree, in.vobs_change, in.inodes, in.obs_decomp, in.lnl_type, false, p_global, 0);

    // Near-frozen tip edges and fast internal edges: a change seen in a sample can
    // no longer arise on its own branch, so the ancestors must carry it instead.
    in.tree.init_edge_rates(global_rates(in.tree));
    for(size_t eid = 0; eid < in.tree.edges.size(); eid++){
        bool tip_edge = in.tree.edges[eid].end < NS;
        in.tree.edge_rates[eid] = in.tree.edge_rates[eid] * (tip_edge ? 1e-6 : 50.0);
    }
    reconstruct_joint_ancestral_state_decomp(
        in.tree, in.vobs_change, in.inodes, in.obs_decomp, in.lnl_type, true, p_var, 0);

    REQUIRE_FALSE(read_file(p_global + ".joint.cn") == read_file(p_var + ".joint.cn"));
}
