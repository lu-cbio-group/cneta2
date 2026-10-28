// Unit tests for turning independent-chain (model 3) states back into copy
// numbers in ancestral state reconstruction (state.cpp).
//
// Chromosome- and site-level states are indices into the list built by
// build_pair_states() (common.hpp), which is ordered by total copy number and
// then by haplotype A. Decoding must use the same order, or reconstructed copy
// numbers are wrong even when the reconstructed states are right.

#include <catch2/catch_test_macros.hpp>

#include <vector>

#include "likelihood.hpp"
#include "parse_cn.hpp"
#include "state.hpp"

namespace {

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

TEST_CASE("a WGD doubles the baseline before chromosome- and site-level changes", "[state]"){
    LNL_TYPE t = lnl_type_for(1, 2, 1);
    REQUIRE(get_cn_from_state_decomp(1, NO_CHANGE_HAPLOTYPE, NO_CHANGE_HAPLOTYPE, t) == 4);

    // one haplotype loses a copy at site level after the WGD
    int loss = encode_change_pair(-1, 0, 2);
    REQUIRE(get_cn_from_state_decomp(1, NO_CHANGE_HAPLOTYPE, loss, t) == 3);
}
