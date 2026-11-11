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

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "common.hpp"
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
    // Callers do not initialise the out-parameters, so state_to_allele_cn must
    // write them for state 0 too. Start from values it would never return.
    int cnA = -99;
    int cnB = -99;
    state_to_allele_cn(0, 4, cnA, cnB);

    REQUIRE(cnA == 0);
    REQUIRE(cnB == 0);
    REQUIRE(state_to_total_cn(0, 4) == 0);
}

TEST_CASE("every state decodes without relying on initialised outputs", "[parse_cn]"){
    for(int cn_max = 1; cn_max <= 8; ++cn_max){
        const int nstate = (cn_max + 1) * (cn_max + 2) / 2;
        for(int state = 0; state < nstate; ++state){
            INFO("cn_max=" << cn_max << " state=" << state);
            int cnA = -99;
            int cnB = -99;
            state_to_allele_cn(state, cn_max, cnA, cnB);
            REQUIRE(cnA >= 0);
            REQUIRE(cnB >= 0);
            REQUIRE(cnA + cnB <= cn_max);
            REQUIRE(allele_cn_to_state(cnA, cnB) == state);
            REQUIRE(state_to_total_cn(state, cn_max) == cnA + cnB);
        }
    }
}

TEST_CASE("haplotype ploidy counts 0/0 sites as zero copies", "[parse_cn]"){
    // rows are chr, site, state; 0/0, 1/1 and 2/0
    const std::vector<std::vector<int>> s_cn{{1, 1, allele_cn_to_state(0, 0)},
                                             {1, 2, allele_cn_to_state(1, 1)},
                                             {1, 3, allele_cn_to_state(2, 0)}};
    map<int, vector<int>> chr_cnA, chr_cnB;
    auto avg = compute_sample_avg_cn_haplotype(s_cn, chr_cnA, chr_cnB, 4);
    REQUIRE_THAT(avg.first, Catch::Matchers::WithinAbs(1.0, 1e-12));
    REQUIRE_THAT(avg.second, Catch::Matchers::WithinAbs(1.0 / 3, 1e-12));
    REQUIRE(chr_cnA.at(1) == std::vector<int>{0, 1, 2});
    REQUIRE(chr_cnB.at(1) == std::vector<int>{0, 1, 0});
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

// ===========================================================================
// Model 3 (DECOMP): converting copy numbers into multi-level changes
// ===========================================================================
//
// Model 3 does not see copy numbers directly: read_data_var_regions_by_chr_change()
// converts each observed copy number into changes at three levels, and the
// likelihood runs one Markov chain per level:
//
//   - WGD count w, shared by the whole sample;
//   - chromosome-level change c, shared by every segment on a chromosome
//     (including segments whose copy number is back at normal);
//   - segment-level change d, for one segment only.
//
// Per haplotype (normal copy number 1) and without WGD, a segment change may
// come after the chromosome change or before it, in which case a chromosome
// gain copies it:
//
//     CN = 1 + c + d          (segment change after the chromosome change)
//     CN = (1 + c)(1 + d)     (segment change before it)
//
// A lost copy cannot change again: if c = -1, then CN = 0 and d = 0. A total
// copy number (normal 2) is the sum over the two haplotypes, with changes on
// one or both. A WGD doubles everything that happened before it.
//
// The tests below cover segment-level changes only (w = 0, c = 0), where the
// conversion must give segment change = copy number - normal. Chromosome-level
// changes and WGD are being reworked (#50) and will be added later.
//
// To print the converted changes:
//
//     CNETA_PRINT_DECOMP=1 code/build/tests/cneta_unit_tests "[decomp]"

namespace {

const int CN_MAX = 8;

// One sample at one segment, decoded from the CN_CHANGE values.
struct Decomposition {
    std::vector<int> observed;    // {total} or {cnA, cnB}
    int num_wgd;
    std::vector<int> chr_change;
    std::vector<int> seg_change;
};

struct Segment {
    int chr;
    int start;                    // 0-based bin ids, as in the data map
    int end;
    std::vector<Decomposition> samples;
};

std::vector<Segment> convert(const std::string& cn_file, int Ns, bool is_total){
    // Ns, cn_max, max_wgd, max_chr_change_haplotype, max_site_change_haplotype,
    // model, is_total, is_rcn, is_bin, incl_all -- cnetml's defaults otherwise.
    const INPUT_PROPERTY property{Ns, CN_MAX, 2, 1, 2, DECOMP, is_total ? 1 : 0, 0, 1, 1};
    INPUT_DATA input{0, 0, 0, {}, {}, {}};
    std::map<int, std::vector<std::vector<int>>> data;

    read_data_var_regions_by_chr_change(data, cn_file, "", "", property, input, 0);

    const auto chr_states = build_pair_states(property.max_chr_change_haplotype);
    const auto site_states = build_pair_states(property.max_site_change_haplotype);

    std::vector<Segment> segments;
    for(const auto& [chr, rows] : data){
        for(const auto& row : rows){
            // chr, start, end, then (cn_state, num_wgd, chr change, site change) per sample
            REQUIRE(row.size() == 3 + 4 * static_cast<size_t>(Ns));
            Segment seg{row[0], row[1], row[2], {}};
            for(int j = 0; j < Ns; ++j){
                const int cn_state = row[3 + 4 * j];
                const int chr_code = row[5 + 4 * j];
                const int site_code = row[6 + 4 * j];

                Decomposition d;
                d.num_wgd = row[4 + 4 * j];
                if(is_total){
                    d.observed = {cn_state};
                    d.chr_change = {chr_code};
                    d.seg_change = {site_code};
                }else{
                    int a = 0, b = 0;
                    state_to_allele_cn(cn_state, CN_MAX, a, b);
                    d.observed = {a, b};
                    d.chr_change = {chr_states.at(chr_code).first, chr_states.at(chr_code).second};
                    d.seg_change = {site_states.at(site_code).first, site_states.at(site_code).second};
                }
                seg.samples.push_back(d);
            }
            segments.push_back(seg);
        }
    }
    return segments;
}

// The decomposition expected for a segment-only change: no WGD, no
// chromosome-level change, and segment change = copy number - normal.
bool is_segment_only(const Decomposition& d, bool is_total){
    const int normal = is_total ? 2 : 1;
    if(d.num_wgd != 0) return false;
    for(size_t h = 0; h < d.observed.size(); ++h){
        if(d.chr_change[h] != 0 || d.seg_change[h] != d.observed[h] - normal) return false;
    }
    return true;
}

std::string join(const std::vector<int>& v, bool sign){
    std::ostringstream os;
    for(size_t i = 0; i < v.size(); ++i){
        if(i) os << "/";
        if(sign && v[i] >= 0) os << "+";
        os << v[i];
    }
    return os.str();
}

// One table per sample.
std::string table(const std::vector<Segment>& segments, int Ns, bool is_total){
    std::ostringstream os;
    os << "Model 3 conversion (" << (is_total ? "total" : "haplotype-specific") << ")\n";
    for(int j = 0; j < Ns; ++j){
        os << "Sample " << j + 1 << "\n"
           << "  chr\tbins\tobserved\tchr change\tseg change\tcheck\n";
        for(const auto& seg : segments){
            const auto& d = seg.samples[j];
            os << "  " << seg.chr << "\t" << seg.start + 1 << "-" << seg.end + 1 << "\t"
               << join(d.observed, false) << "\t\t"
               << join(d.chr_change, true) << "\t\t"
               << join(d.seg_change, true) << "\t\t"
               << (is_segment_only(d, is_total) ? "ok" : "FAIL") << "\n";
        }
    }
    return os.str();
}

// Convert an input with segment-level changes only, print if asked, and require
// every sample at every segment to convert as a segment-only change.
void require_segment_only(const std::string& cn_file, int Ns, bool is_total){
    const auto segments = convert(cn_file, Ns, is_total);
    REQUIRE_FALSE(segments.empty());

    const std::string report = table(segments, Ns, is_total);
    if(std::getenv("CNETA_PRINT_DECOMP")) std::cout << "\n" << cn_file << "\n" << report;
    INFO(report);

    for(const auto& seg : segments){
        for(int j = 0; j < Ns; ++j){
            INFO("sample " << j + 1 << ", chromosome " << seg.chr << ", bins "
                 << seg.start + 1 << "-" << seg.end + 1);
            CHECK(is_segment_only(seg.samples[j], is_total));
        }
    }
}

// Synthetic input: 4 chromosomes of 10 bins, every bin 1/1 (total 2) unless a
// FocalChange sets it. Samples 1 and 2 carry focal changes of one or two bins
// on otherwise normal chromosomes; sample 3 is normal. The changes are small
// enough that no chromosome mean, and no sample ploidy, rounds away from normal.
struct FocalChange { int sample; int chr; int first_bin; int last_bin; int cnA; int cnB; };

const std::vector<FocalChange> SEGMENT_ONLY_CHANGES{
    {1, 1, 3, 4, 2, 1},   // gain (+1/0)
    {1, 2, 5, 5, 0, 1},   // loss (-1/0)
    {1, 3, 7, 7, 3, 1},   // amplification (+2/0)
    {1, 4, 2, 2, 0, 0},   // homozygous deletion (-1/-1)
    {2, 1, 8, 8, 2, 0},   // copy-neutral LOH (+1/-1, total unchanged)
    {2, 2, 1, 2, 1, 2},   // gain (0/+1)
    {2, 3, 4, 4, 1, 3},   // amplification (0/+2)
    {2, 4, 9, 9, 2, 2},   // gain of both haplotypes (+1/+1)
};

std::string write_segment_only_input(bool is_total){
    const auto dir = std::filesystem::temp_directory_path() / "cneta_test_parse_cn";
    std::filesystem::create_directories(dir);
    const auto path = dir / (is_total ? "segment-only-cn.txt" : "segment-only-haplotype-cn.txt");
    std::ofstream out(path);
    for(int sample = 1; sample <= 3; ++sample){
        for(int chr = 1; chr <= 4; ++chr){
            for(int bin = 1; bin <= 10; ++bin){
                int a = 1, b = 1;
                for(const auto& f : SEGMENT_ONLY_CHANGES){
                    if(f.sample == sample && f.chr == chr && bin >= f.first_bin && bin <= f.last_bin){
                        a = f.cnA;
                        b = f.cnB;
                    }
                }
                out << sample << "\t" << chr << "\t" << bin << "\t";
                if(is_total) out << a + b << "\n";
                else         out << a << "\t" << b << "\n";
            }
        }
    }
    return path.string();
}

}  // namespace

TEST_CASE("model 3 conversion recovers segment-only changes from total copy numbers", "[parse_cn][decomp]"){
    require_segment_only(write_segment_only_input(true), 3, true);
}

TEST_CASE("model 3 conversion recovers segment-only changes from haplotype copy numbers", "[parse_cn][decomp]"){
    require_segment_only(write_segment_only_input(false), 3, false);
}
