// Unit tests for pieces of the ML-RLC (bsr_mode=3) shift-edge search in
// optimization.hpp.
//
// The random-local-clock shift-edge search must never pick the root-LUCA
// edge as a shift edge when estimate_bsr0_first=false, because that edge's
// rate is instead estimated directly as mu0 (see optimization.cpp around
// luca_eid_bsr1). find_luca_eid identifies the edge, get_active_bsr_eids
// excludes it from the candidate pool when asked to, and
// pick_best_improving is the comparison step that actually decides which
// candidate (LUCA or otherwise) wins each round of the search.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <string>
#include <vector>

#include "optimization.hpp"
#include "tree_op.hpp"

namespace {

std::string data_path(const std::string& name){
    return std::string(CNETA_TEST_DATA_DIR) + "/" + name;
}

const int NS = 3;
const int NLEAF = NS + 1;

}  // namespace

TEST_CASE("find_luca_eid finds the length>0 edge leaving the root", "[optimization]"){
    evo_tree tree = read_tree_info(data_path("tiny-tree.txt"), NS);

    int luca_eid = find_luca_eid(tree);

    REQUIRE(luca_eid >= 0);
    REQUIRE(tree.edges[luca_eid].start == tree.root_node_id);
    REQUIRE(tree.edges[luca_eid].length > 0);

    // The root's other outgoing edge (the zero-length one, if it has one) is
    // not the LUCA edge, even though it also starts at the root.
    for(int eid = 0; eid < (int)tree.edges.size(); ++eid){
        if(tree.edges[eid].start == tree.root_node_id && eid != luca_eid){
            REQUIRE(tree.edges[eid].length == 0);
        }
    }
}

TEST_CASE("get_active_bsr_eids excludes the normal-sample edge", "[optimization]"){
    evo_tree tree = read_tree_info(data_path("tiny-tree.txt"), NS);

    std::vector<int> active = get_active_bsr_eids(tree);

    for(int eid : active){
        REQUIRE_FALSE(tree.is_normal_sample_edge(eid));
    }
    // One normal-sample edge is excluded from the full edge set.
    REQUIRE(active.size() == tree.edges.size() - 1);
}

TEST_CASE("get_active_bsr_eids also excludes the LUCA edge when asked", "[optimization]"){
    evo_tree tree = read_tree_info(data_path("tiny-tree.txt"), NS);
    int luca_eid = find_luca_eid(tree);

    std::vector<int> active_with_luca = get_active_bsr_eids(tree);
    std::vector<int> active_without_luca = get_active_bsr_eids(tree, luca_eid);

    // estimate_bsr0_first=true path: LUCA is a normal shift-edge candidate.
    REQUIRE(std::find(active_with_luca.begin(), active_with_luca.end(), luca_eid)
            != active_with_luca.end());

    // estimate_bsr0_first=false path: LUCA can never be picked as a shift
    // edge, because its rate is estimated directly as mu0 instead.
    REQUIRE(std::find(active_without_luca.begin(), active_without_luca.end(), luca_eid)
            == active_without_luca.end());
    REQUIRE(active_without_luca.size() == active_with_luca.size() - 1);
}

// pick_best_improving is the comparison step inside stepwise_search_shift_edges
// that actually decides which candidate edge (if any) gets added as a shift
// edge — see optimization.cpp:2424 (forward selection) and :2537 (backward
// removal). It takes already-computed scores rather than computing them
// itself, so it can be tested with hand-picked numbers instead of running a
// real likelihood optimization.

TEST_CASE("pick_best_improving requires beating the baseline by more than eps", "[optimization]"){
    double cur_score = 5.0;
    double eps = 1e-6;

    SECTION("a score exactly at cur_score + eps does not count as an improvement"){
        std::vector<char> evaluated{1};
        std::vector<double> scores{cur_score + eps};
        std::vector<int> ids{0};

        REQUIRE(pick_best_improving(1, evaluated, scores, ids, cur_score, eps) == -1);
    }

    SECTION("a score just above cur_score + eps is picked"){
        std::vector<char> evaluated{1};
        std::vector<double> scores{cur_score + eps * 2};
        std::vector<int> ids{0};

        REQUIRE(pick_best_improving(1, evaluated, scores, ids, cur_score, eps) == 0);
    }
}

TEST_CASE("pick_best_improving picks the highest-scoring candidate", "[optimization]"){
    // Candidate order in the arrays does not matter; index 1 has the best score.
    std::vector<char> evaluated{1, 1, 1};
    std::vector<double> scores{7.0, 9.0, 8.0};
    std::vector<int> ids{10, 11, 12};

    int best = pick_best_improving(3, evaluated, scores, ids, 5.0, 1e-6);

    REQUIRE(best == 1);
}

TEST_CASE("pick_best_improving skips candidates that were not evaluated", "[optimization]"){
    // Index 1 has the top score but was never evaluated (e.g. already selected
    // in a previous round), so index 2 should win instead.
    std::vector<char> evaluated{1, 0, 1};
    std::vector<double> scores{7.0, 9.0, 8.0};
    std::vector<int> ids{10, 11, 12};

    int best = pick_best_improving(3, evaluated, scores, ids, 5.0, 1e-6);

    REQUIRE(best == 2);
}

TEST_CASE("pick_best_improving returns -1 when nothing improves on the baseline", "[optimization]"){
    std::vector<char> evaluated{1, 1};
    std::vector<double> scores{3.0, 4.0};
    std::vector<int> ids{0, 1};

    int best = pick_best_improving(2, evaluated, scores, ids, 5.0, 1e-6);

    REQUIRE(best == -1);
}

TEST_CASE("pick_best_improving decides whether the LUCA edge is picked as a shift edge",
          "[optimization]"){
    // Mirrors the forward-selection call inside stepwise_search_shift_edges:
    // the LUCA edge is just another candidate id competing on score. Scores
    // here are made up by hand -- in the real search they would come from
    // compute_rlc_ic on a real likelihood run.
    const int luca_index = 0;
    const int luca_id = 4;  // matches find_luca_eid's result on tiny-tree.txt

    SECTION("LUCA scores clearly higher than every other candidate: LUCA is picked"){
        std::vector<char> evaluated{1, 1, 1};
        std::vector<double> scores{10.0, 6.0, 7.0};
        std::vector<int> ids{luca_id, 1, 2};

        int best = pick_best_improving(3, evaluated, scores, ids, 5.0, 1e-6);

        REQUIRE(best == luca_index);
    }

    SECTION("LUCA ties with another candidate but has the larger id: LUCA loses the tie"){
        std::vector<char> evaluated{1, 1};
        std::vector<int> ids{luca_id, 1};  // luca_id (4) > 1
        std::vector<double> scores{9.0, 9.0};

        int best = pick_best_improving(2, evaluated, scores, ids, 5.0, 1e-6);

        REQUIRE(best != luca_index);
        REQUIRE(ids[best] == 1);
    }

    SECTION("LUCA ties with another candidate but has the smaller id: LUCA wins the tie"){
        std::vector<char> evaluated{1, 1};
        std::vector<int> ids{luca_id, 9};  // luca_id (4) < 9
        std::vector<double> scores{9.0, 9.0};

        int best = pick_best_improving(2, evaluated, scores, ids, 5.0, 1e-6);

        REQUIRE(best == luca_index);
    }
}
