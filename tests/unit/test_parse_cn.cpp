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
#include <functional>
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
// The tests below cover inputs without WGD (w = 0), with segment-level changes
// only and with chromosome- and segment-level changes together. They use only
// cases with a single valid answer, so the expected segment change is the
// copy number minus 1 + c per haplotype (0 on a lost copy). Cases with two
// valid answers, such as 4/1 on a chromosome gained on A (segment +1 before the
// gain, or +2 after it), and WGD are left for the rework in #50.
//
// Each test writes an input as copy numbers (Regions), converts it, and
// compares every segment with an Expectation: a rule giving the expected WGD
// count, chromosome change and segment change. To add a test, e.g. with WGD,
// write its Regions and an Expectation for it; the rest is shared.
//
// To print the converted changes:
//
//     CNETA_PRINT_DECOMP=1 code/build/tests/cneta_unit_tests "[decomp]"

namespace {

const int CN_MAX = 8;
const int NS = 3;      // samples in each synthetic input: two tumours, then a normal sample
const int NCHR = 4;    // chromosomes in each synthetic input
const int NBIN = 10;   // bins per chromosome

// ---------------------------------------------------------------------------
// Writing inputs
// ---------------------------------------------------------------------------

// A run of bins at the given copy numbers. Regions are applied in order, so a
// later one overrides an earlier one; bins outside every region are 1/1.
struct Region { int sample; int chr; int first_bin; int last_bin; int cnA; int cnB; };

Region whole_chr(int sample, int chr, int cnA, int cnB){
    return {sample, chr, 1, NBIN, cnA, cnB};
}

// Writes the input as total (sample, chr, bin, cn) or haplotype-specific
// (sample, chr, bin, cnA, cnB) copy numbers and returns its path.
std::string write_input(const std::string& name, const std::vector<Region>& regions, bool is_total){
    const auto dir = std::filesystem::temp_directory_path() / "cneta_test_parse_cn";
    std::filesystem::create_directories(dir);
    const auto path = dir / (name + (is_total ? "-cn.txt" : "-haplotype-cn.txt"));
    std::ofstream out(path);
    for(int sample = 1; sample <= NS; ++sample){
        for(int chr = 1; chr <= NCHR; ++chr){
            for(int bin = 1; bin <= NBIN; ++bin){
                int a = 1, b = 1;
                for(const auto& r : regions){
                    if(r.sample == sample && r.chr == chr && bin >= r.first_bin && bin <= r.last_bin){
                        a = r.cnA;
                        b = r.cnB;
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

// ---------------------------------------------------------------------------
// Reading the conversion
// ---------------------------------------------------------------------------

// One sample at one segment, decoded from the CN_CHANGE values. Changes are
// {total} or {haplotype A, haplotype B}.
struct Decomposition {
    std::vector<int> observed;
    int num_wgd;
    std::vector<int> chr_change;
    bool chr_tagged;              // halved for a WGD (stored as a multiple of CHANGE_CHR)
    std::vector<int> seg_change;
};

struct Segment {
    int chr;
    int start;                    // 0-based bin ids, as in the data map
    int end;
    std::vector<Decomposition> samples;
};

std::vector<Segment> convert(const std::string& cn_file, bool is_total){
    // Ns, cn_max, max_wgd, max_chr_change_haplotype, max_site_change_haplotype,
    // model, is_total, is_rcn, is_bin, incl_all -- cnetml's defaults otherwise.
    const INPUT_PROPERTY property{NS, CN_MAX, 2, 1, 2, DECOMP, is_total ? 1 : 0, 0, 1, 1};
    INPUT_DATA input{0, 0, 0, {}, {}, {}};
    std::map<int, std::vector<std::vector<int>>> data;
    read_data_var_regions_by_chr_change(data, cn_file, "", "", property, input, 0);

    const auto chr_states = build_pair_states(property.max_chr_change_haplotype);
    const auto site_states = build_pair_states(property.max_site_change_haplotype);
    auto as_vector = [](const std::pair<int, int>& p){ return std::vector<int>{p.first, p.second}; };

    std::vector<Segment> segments;
    for(const auto& [chr, rows] : data){
        for(const auto& row : rows){
            // chr, start, end, then (cn_state, num_wgd, chr change, site change) per sample
            REQUIRE(row.size() == 3 + 4 * static_cast<size_t>(NS));
            Segment seg{row[0], row[1], row[2], {}};
            for(int j = 0; j < NS; ++j){
                const int cn_state = row[3 + 4 * j];
                const int chr_code = row[5 + 4 * j];
                const int site_code = row[6 + 4 * j];
                const bool tagged = chr_code != 0 && chr_code % CHANGE_CHR == 0;
                const int chr_value = tagged ? chr_code / CHANGE_CHR : chr_code;

                Decomposition d{{}, row[4 + 4 * j], {}, tagged, {}};
                if(is_total){
                    d.observed = {cn_state};
                    d.chr_change = {chr_value};
                    d.seg_change = {site_code};
                }else{
                    int a = 0, b = 0;
                    state_to_allele_cn(cn_state, CN_MAX, a, b);
                    d.observed = {a, b};
                    d.chr_change = as_vector(chr_states.at(chr_value));
                    d.seg_change = as_vector(site_states.at(site_code));
                }
                seg.samples.push_back(d);
            }
            segments.push_back(seg);
        }
    }
    return segments;
}

// ---------------------------------------------------------------------------
// Expected answers
// ---------------------------------------------------------------------------

struct Expected {
    int num_wgd;
    std::vector<int> chr_change;
    std::vector<int> seg_change;
};

// Gives the expected decomposition of one sample at one segment, from its
// sample number (1-based), chromosome and observed copy numbers.
using Expectation = std::function<Expected(int sample, int chr, const std::vector<int>& observed)>;

// Chromosome-level change of one sample on one chromosome, per haplotype.
struct ChrChange { int sample; int chr; int cA; int cB; };

// The answer without WGD for cases with a single valid decomposition:
// each chromosome's change as given (0 if not listed),
// and per haplotype a segment change of copy number - (1 + c), or 0 on a lost copy.
// For total copy numbers, both summed over the two haplotypes.
Expectation without_wgd(const std::vector<ChrChange>& chr_changes){
    return [chr_changes](int sample, int chr, const std::vector<int>& observed){
        int cA = 0, cB = 0;
        for(const auto& c : chr_changes){
            if(c.sample == sample && c.chr == chr){ cA = c.cA; cB = c.cB; }
        }
        if(observed.size() == 1) return Expected{0, {cA + cB}, {observed[0] - (2 + cA + cB)}};
        auto seg = [](int cn, int c){ return c == -1 ? 0 : cn - (1 + c); };
        return Expected{0, {cA, cB}, {seg(observed[0], cA), seg(observed[1], cB)}};
    };
}

// ---------------------------------------------------------------------------
// Checking
// ---------------------------------------------------------------------------

std::string join(const std::vector<int>& v, bool sign){
    std::ostringstream os;
    for(size_t i = 0; i < v.size(); ++i){
        if(i) os << "/";
        if(sign && v[i] >= 0) os << "+";
        os << v[i];
    }
    return os.str();
}

bool matches(const Decomposition& d, const Expected& e){
    return d.num_wgd == e.num_wgd && d.chr_change == e.chr_change && d.seg_change == e.seg_change;
}

// One table per sample. '*' marks a chromosome change halved for a WGD.
std::string table(const std::vector<Segment>& segments, const std::vector<std::vector<Expected>>& want){
    std::ostringstream os;
    for(int j = 0; j < NS; ++j){
        os << "Sample " << j + 1 << "\n"
           << "  chr\tbins\tobserved\tWGD\tchr change\tseg change\tcheck\n";
        for(size_t i = 0; i < segments.size(); ++i){
            const auto& seg = segments[i];
            const auto& d = seg.samples[j];
            os << "  " << seg.chr << "\t" << seg.start + 1 << "-" << seg.end + 1 << "\t"
               << join(d.observed, false) << "\t\t" << d.num_wgd << "\t"
               << join(d.chr_change, true) << (d.chr_tagged ? "*" : "") << "\t\t"
               << join(d.seg_change, true) << "\t\t"
               << (matches(d, want[i][j]) ? "ok" : "FAIL") << "\n";
        }
    }
    return os.str();
}

// Convert an input and require every sample at every segment to match expect.
void check_conversion(const std::string& cn_file, bool is_total, const Expectation& expect){
    const auto segments = convert(cn_file, is_total);
    REQUIRE_FALSE(segments.empty());

    std::vector<std::vector<Expected>> want(segments.size());
    for(size_t i = 0; i < segments.size(); ++i){
        for(int j = 0; j < NS; ++j){
            want[i].push_back(expect(j + 1, segments[i].chr, segments[i].samples[j].observed));
        }
    }

    const std::string report = table(segments, want);
    if(std::getenv("CNETA_PRINT_DECOMP")) std::cout << "\n" << cn_file << "\n" << report;
    INFO(report);

    for(size_t i = 0; i < segments.size(); ++i){
        const auto& seg = segments[i];
        for(int j = 0; j < NS; ++j){
            const auto& d = seg.samples[j];
            const auto& e = want[i][j];
            INFO("sample " << j + 1 << ", chromosome " << seg.chr << ", bins "
                 << seg.start + 1 << "-" << seg.end + 1);
            CHECK(d.num_wgd == e.num_wgd);
            CHECK(d.chr_change == e.chr_change);
            CHECK(d.seg_change == e.seg_change);
        }
    }
}

// ---------------------------------------------------------------------------
// Inputs
// ---------------------------------------------------------------------------

// Segment-level changes only: focal changes of one or two bins on otherwise
// normal chromosomes.
const std::vector<Region> SEGMENT_ONLY_INPUT{
    {1, 1, 3, 4, 2, 1},   // gain (+1/0)
    {1, 2, 5, 5, 0, 1},   // loss (-1/0)
    {1, 3, 7, 7, 3, 1},   // amplification (+2/0)
    {1, 4, 2, 2, 0, 0},   // homozygous deletion (-1/-1)
    {2, 1, 8, 8, 2, 0},   // copy-neutral LOH (+1/-1, total unchanged)
    {2, 2, 1, 2, 1, 2},   // gain (0/+1)
    {2, 3, 4, 4, 1, 3},   // amplification (0/+2)
    {2, 4, 9, 9, 2, 2},   // gain of both haplotypes (+1/+1)
};

// Chromosome- and segment-level changes together, each segment with a single
// valid decomposition. Changes are small enough that no chromosome mean, and
// no sample ploidy, rounds the wrong way.
const std::vector<ChrChange> CHR_CHANGES{
    {1, 1, 1, 0},    // chr 1 gained on A
    {1, 2, -1, 0},   // chr 2 lost on A
    {2, 3, 0, 1},    // chr 3 gained on B
    {2, 4, 0, -1},   // chr 4 lost on B
};
const std::vector<Region> CHR_AND_SEGMENT_INPUT{
    whole_chr(1, 1, 2, 1),
    whole_chr(1, 2, 0, 1),
    whole_chr(2, 3, 1, 2),
    whole_chr(2, 4, 1, 0),
    {1, 1, 3, 3, 3, 1},     // gain after the chromosome gain: seg +1/0
    {1, 1, 10, 10, 1, 1},   // back at normal: seg -1/0, cancelling the gain
    {2, 1, 10, 10, 2, 1},   // gain on an unchanged chromosome: seg +1/0 (keeps bin 10
                            // from being dropped as the same in every sample)
    {1, 2, 5, 5, 0, 2},     // gain on the remaining copy: seg 0/+1
    {2, 3, 7, 7, 2, 2},     // gain on the other haplotype: seg +1/0
    {2, 4, 9, 9, 2, 0},     // gain on the remaining copy: seg +1/0
};

}  // namespace

TEST_CASE("model 3 conversion recovers segment-only changes from total copy numbers", "[parse_cn][decomp]"){
    check_conversion(write_input("segment-only", SEGMENT_ONLY_INPUT, true), true, without_wgd({}));
}

TEST_CASE("model 3 conversion recovers segment-only changes from haplotype-specific copy numbers", "[parse_cn][decomp]"){
    check_conversion(write_input("segment-only", SEGMENT_ONLY_INPUT, false), false, without_wgd({}));
}

TEST_CASE("model 3 conversion recovers chromosome- and segment-level changes from total copy numbers", "[parse_cn][decomp]"){
    check_conversion(write_input("chr-and-segment", CHR_AND_SEGMENT_INPUT, true), true, without_wgd(CHR_CHANGES));
}

TEST_CASE("model 3 conversion recovers chromosome- and segment-level changes from haplotype-specific copy numbers", "[parse_cn][decomp]"){
    check_conversion(write_input("chr-and-segment", CHR_AND_SEGMENT_INPUT, false), false, without_wgd(CHR_CHANGES));
}
