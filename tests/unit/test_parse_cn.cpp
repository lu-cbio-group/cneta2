// Unit tests for the input parsing and state-encoding helpers in parse_cn.cpp.
//
// The haplotype-specific state encoding is the single most load-bearing
// convention in the codebase: states are ordered by total copy number, then
// by cnA, i.e.
//
//     state:  0    1    2    3    4    5    6   ...
//     cnA/cnB 0/0  0/1  1/0  0/2  1/1  2/0  0/3 ...
//
// Every likelihood calculation depends on it, so the round-trip below is the
// test most likely to catch an accidental change during the libcneta work.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <string>
#include <vector>

#include "parse_cn.hpp"

namespace {

std::string data_path(const std::string& name){
    return std::string(CNETA_TEST_DATA_DIR) + "/" + name;
}

}  // namespace

TEST_CASE("allele_cn_to_state matches the documented ordering", "[parse_cn]"){
    REQUIRE(allele_cn_to_state(0, 0) == 0);
    REQUIRE(allele_cn_to_state(0, 1) == 1);
    REQUIRE(allele_cn_to_state(1, 0) == 2);
    REQUIRE(allele_cn_to_state(0, 2) == 3);
    REQUIRE(allele_cn_to_state(1, 1) == 4);   // NORM_ALLElE_STATE
    REQUIRE(allele_cn_to_state(2, 0) == 5);
    REQUIRE(allele_cn_to_state(0, 3) == 6);
    REQUIRE(allele_cn_to_state(4, 0) == 14);
}

TEST_CASE("the normal diploid state is state 4", "[parse_cn]"){
    // common.hpp hard-codes NORM_ALLElE_STATE = 4 as "1/1". If the encoding
    // ever shifts, that constant becomes silently wrong.
    REQUIRE(allele_cn_to_state(1, 1) == NORM_ALLElE_STATE);
}

TEST_CASE("state and haplotype-specific copy number round-trip", "[parse_cn]"){
    const int cn_max = 4;

    for(int cnA = 0; cnA <= cn_max; ++cnA){
        for(int cnB = 0; cnB + cnA <= cn_max; ++cnB){
            const int state = allele_cn_to_state(cnA, cnB);

            INFO("cnA=" << cnA << " cnB=" << cnB << " state=" << state);

            // Round-trip back to the pair we started from.
            int back_a = 0;
            int back_b = 0;
            state_to_allele_cn(state, cn_max, back_a, back_b);
            REQUIRE(back_a == cnA);
            REQUIRE(back_b == cnB);

            // And the total must agree with the two halves.
            REQUIRE(state_to_total_cn(state, cn_max) == cnA + cnB);
        }
    }
}

TEST_CASE("state 0 decodes to a homozygous deletion", "[parse_cn]"){
    // state_to_allele_cn returns early for state 0 without writing to its
    // out-parameters, so callers must initialise them. Both conventions agree
    // that the answer is 0/0.
    int cnA = 0;
    int cnB = 0;
    state_to_allele_cn(0, 4, cnA, cnB);

    REQUIRE(cnA == 0);
    REQUIRE(cnB == 0);
    REQUIRE(state_to_total_cn(0, 4) == 0);
}

TEST_CASE("states are contiguous and unique up to cn_max", "[parse_cn]"){
    const int cn_max = 4;
    const int nstate = (cn_max + 1) * (cn_max + 2) / 2;

    std::vector<int> states;
    for(int cnA = 0; cnA <= cn_max; ++cnA){
        for(int cnB = 0; cnB + cnA <= cn_max; ++cnB){
            states.push_back(allele_cn_to_state(cnA, cnB));
        }
    }

    std::sort(states.begin(), states.end());

    // Exactly the numbers 0..nstate-1, each once. This is what lets the
    // rate matrices be indexed directly by state.
    REQUIRE(static_cast<int>(states.size()) == nstate);
    for(int i = 0; i < nstate; ++i){
        REQUIRE(states[i] == i);
    }
}

TEST_CASE("get_num_wgd reads WGD count off the ploidy", "[parse_cn]"){
    // WGD_CUTOFF is 3.0: above it means one doubling, above 2x it means two.
    std::vector<double> avg_cn{1.9, 2.0, 3.0, 3.5, 6.0, 6.5};
    std::vector<int> nwgd;

    get_num_wgd(avg_cn, nwgd, 0);

    REQUIRE(nwgd == std::vector<int>{0, 0, 0, 1, 1, 2});
}

TEST_CASE("read_cn parses a gzipped copy-number file", "[parse_cn]"){
    // The fixture holds 4 samples (3 tumour regions + normal) over 50 sites,
    // in the four-column total-copy-number layout described in the docs.
    const int Ns = 3;
    const int cn_max = 4;
    int num_total_bins = 0;

    std::vector<std::vector<std::vector<int>>> s_info;
    read_cn(s_info, data_path("tiny-cn.txt.gz"), Ns, num_total_bins, cn_max,
            /*is_total=*/1, /*is_rcn=*/0, /*debug=*/0);

    REQUIRE(num_total_bins == 50);

    // One entry per tumour sample; the normal is skipped on read.
    REQUIRE(s_info.size() == static_cast<size_t>(Ns));

    SECTION("every record is a well-formed (chr, site, cn) triple"){
        // s_info is indexed [sample][record][chr, site_id, cn] -- the inner
        // vector is a triple, not a list of copy numbers.
        for(const auto& sample : s_info){
            REQUIRE(sample.size() == static_cast<size_t>(num_total_bins));

            for(const auto& record : sample){
                REQUIRE(record.size() == 3);

                const int chr = record[0];
                const int site = record[1];
                const int cn = record[2];

                REQUIRE(chr >= 1);
                REQUIRE(chr <= NUM_CHR);
                REQUIRE(site >= 1);
                REQUIRE(cn >= 0);
                REQUIRE(cn <= cn_max);   // read_cn clamps anything above
            }
        }
    }
}

TEST_CASE("read_cn gives the same answer for plain and gzipped input", "[parse_cn]"){
    // gzstream is meant to be transparent; this is the test that says so.
    const int Ns = 3;
    const int cn_max = 4;

    int bins_gz = 0;
    std::vector<std::vector<std::vector<int>>> from_gz;
    read_cn(from_gz, data_path("tiny-cn.txt.gz"), Ns, bins_gz, cn_max, 1, 0, 0);

    int bins_plain = 0;
    std::vector<std::vector<std::vector<int>>> from_plain;
    read_cn(from_plain, data_path("tiny-cn.txt"), Ns, bins_plain, cn_max, 1, 0, 0);

    REQUIRE(bins_gz == bins_plain);
    REQUIRE(from_gz == from_plain);
}
