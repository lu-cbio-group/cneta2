// Unit tests for rescaling the model-3 likelihood tables (likelihood.cpp).
//
// With many tips the conditional likelihoods of a node can fall below the
// smallest double and become 0, which turns a valid tree into a penalty
// likelihood. Nodes are therefore rescaled by a power of two and the log of
// the scale is kept in the last row of the table; the root log likelihood
// must come out the same as without rescaling.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <cmath>
#include <vector>

#include "likelihood.hpp"

using Catch::Matchers::WithinRel;

namespace {

// Table for a tree with two leaves: rows for the nodes, then the scale row.
vector<vector<double>> make_table(int dim){
    return vector<vector<double>>(4, vector<double>(dim, 0.0));
}

}  // namespace


TEST_CASE("a small node is rescaled exactly and its log scale recorded", "[likelihood]"){
    auto table = make_table(3);
    const vector<double> orig{1e-200, 3e-203, 0.0};
    table[1] = orig;

    rescale_node_lnl(table, 1);

    double scale = get_lnl_log_scale(table);
    REQUIRE(scale < 0.0);
    REQUIRE(table[1][0] >= 0.5);
    REQUIRE(table[1][0] < 1.0);
    REQUIRE(table[1][2] == 0.0);
    for(int s = 0; s < 2; ++s){
        REQUIRE_THAT(log(table[1][s]) + scale, WithinRel(log(orig[s]), 1e-14));
    }
}


TEST_CASE("nodes above the threshold are left unchanged", "[likelihood]"){
    auto table = make_table(2);
    table[1] = {0.25, 1e-100};
    rescale_node_lnl(table, 1);
    REQUIRE(table[1][0] == 0.25);
    REQUIRE(table[1][1] == 1e-100);
    REQUIRE(get_lnl_log_scale(table) == 0.0);

    table[2] = {0.0, 0.0};    // all-zero rows stay zero, so impossible data keep the penalty
    rescale_node_lnl(table, 2);
    REQUIRE(get_lnl_log_scale(table) == 0.0);
}


TEST_CASE("a root probability below the smallest double stays finite", "[likelihood]"){
    // The root (row 2 for nleaf = 2) is the product of its two children.
    // 1e-200 * 1e-200 underflows to 0; with rescaling the product is kept.
    auto table = make_table(1);
    table[0] = {1e-200};
    table[1] = {1e-200};
    rescale_node_lnl(table, 0);
    rescale_node_lnl(table, 1);
    table[2] = {table[0][0] * table[1][0]};
    rescale_node_lnl(table, 2);

    REQUIRE_THAT(extract_tree_lnl_change(table, 2, 0, 0), WithinRel(2 * log(1e-200), 1e-14));
}
