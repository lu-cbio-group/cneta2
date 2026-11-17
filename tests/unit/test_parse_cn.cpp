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
//   - segment-level change s, for one segment only.
//
// Chromosome-level events are called gains and losses, segment-level events
// duplications and deletions.
//
// Per haplotype (normal copy number 1) and without WGD, a segment change may
// come after the chromosome change or before it, in which case a chromosome
// gain copies it:
//
//     CN = 1 + c + s          (segment change after the chromosome change)
//     CN = (1 + c)(1 + s)     (segment change before it)
//
// A lost copy cannot change again: if c = -1, then CN = 0 and s = 0. A total
// copy number (normal 2) is the sum over the two haplotypes, with changes on
// one or both. A WGD doubles everything that happened before it.
//
// The tests below cover inputs with no, one or two WGDs (w = 0, 1, 2), with
// segment-level changes, chromosome-level changes, or both.
//
// The conversion keeps only the net change at each level, not the order of
// the events behind it: real data cannot tell orders apart (a large
// duplication followed by a small deletion looks like a smaller duplication),
// and trying every order is infeasible. So at each level an increase is read
// as the fewest gains or duplications, and a decrease as the fewest losses or
// deletions, never a mix of both, where a WGD, or a chromosome gain that
// doubles a haplotype, doubles every event before it. E.g. 4/1 on a
// chromosome gained on A is segment +1/0 (one duplication before the
// chromosome gain), not +2/0 (two after it), and 5/1 is +2/0 (one before,
// one after). See expected_changes.
//
// A history mixing duplications and deletions can be just as short, e.g.
// 7/2 on a chromosome at 4/2 after one WGD is +2 +1 (expected) or +4 -1. Such
// histories are not considered, so every case has a single expected answer.
// The order of events written above each input is one history giving that
// answer, there to make it easy to check, not the only possible one.
// Tests the current code fails are tagged [!shouldfail] with a note on #50.
//
// Each test writes an input as copy numbers (Regions), converts it, and
// compares every segment with an Expectation: a rule giving the expected WGD
// count, chromosome change and segment change. To add a test, write its Regions
// and an Expectation for it (usually expected_changes); the rest is shared.
//
// To print the converted changes:
//
//     CNETA_PRINT_DECOMP=1 code/build/tests/cneta_unit_tests "[decomp]"

namespace {

const int CN_MAX = 10;   // as in typical runs; higher values slow the program
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
    const INPUT_PROPERTY property{.Ns = NS, .cn_max = CN_MAX, .max_wgd = 2,
                                 .max_chr_change_haplotype = 1, .max_site_change_haplotype = 2,
                                 .model = DECOMP, .is_total = is_total ? 1 : 0, .is_rcn = 0,
                                 .is_bin = 1, .incl_all = 1};
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

// Copies gained (+) or lost (-) on a whole chromosome of one sample, per
// haplotype, as seen in the observed copy numbers (so already doubled by a WGD
// if the change came before it).
struct ChrChange { int sample; int chr; int cA; int cB; };

// Number of WGDs in each sample (1-based sample = index + 1).
const std::vector<int> NO_WGD{0, 0, 0};
const std::vector<int> ONE_WGD{1, 1, 0};   // both tumours, not the normal sample

// The fewest events adding up to a change of copies, signed: all gains (or
// duplications) for an increase, all losses (or deletions) for a decrease, at
// the chromosome (or segment) level. An event adds (or removes) 1 copy,
// doubled by each doubling after it: sizes 1, 2, ..., 2^doublings.
// E.g. +3 with one doubling is two events (+1 before it, +1 after), not three.
// doublings is the most doublings an event can have had (see expected_changes),
// so it sets the largest event size.
int fewest_events(int change, int doublings){
    int left = std::abs(change), events = 0;
    // Take as many of the largest events as fit, then the next size down. As
    // the sizes are powers of 2, this gives the fewest events.
    for(int size = 1 << doublings; size >= 1; size /= 2){
        events += left / size;
        left %= size;
    }
    return change < 0 ? -events : events;
}

// The expected decomposition (net changes, read as above), given each sample's
// WGD count (0, 1 or 2) and the chromosome changes (0 if not listed).
//
// Per haplotype, with w WGDs:
//   normal copy number      = 2^w (1 without WGD, 2 after one, 4 after two)
//   chromosome copy number  = normal + chromosome change
//   segment change          = fewest duplications (or deletions) making up
//                             segment copy number - chromosome copy number,
//                             or 0 if the chromosome copy number is 0 (a lost
//                             copy cannot change again)
//
// Events are counted with fewest_events. A WGD doubles earlier events, and so
// does a chromosome gain that doubles the haplotype (chromosome at twice the
// normal), so the events of a segment change may come before or after them,
// but are all duplications or all deletions. E.g. without WGD, 5/1 on a
// chromosome at 2/1 is segment +2/0: one duplication before the chromosome
// gain (copied by it) and one after.
//
// For total copy numbers, the same is applied to the sums over the two
// haplotypes, so tests change one haplotype per chromosome or segment. The
// segment change is 0 only if the chromosome's total is 0 (both haplotypes
// lost): a total above 0 cannot show that one haplotype is at 0.
Expectation expected_changes(const std::vector<int>& num_wgd, const std::vector<ChrChange>& chr_changes){
    return [num_wgd, chr_changes](int sample, int chr, const std::vector<int>& observed){
        const int w = num_wgd.at(sample - 1);
        const int normal = 1 << w;

        int cA = 0, cB = 0;
        for(const auto& c : chr_changes){
            if(c.sample == sample && c.chr == chr){ cA = c.cA; cB = c.cB; }
        }
        // The most doublings a segment event can have had, on a haplotype whose
        // chromosome changed by c: the w WGDs, plus one if the chromosome gain
        // doubled the haplotype (chromosome at twice the normal, so gained
        // before every WGD). E.g. with one WGD, a chromosome at 4 copies on the
        // haplotype gives 2: a duplication then adds 4, 2 or 1 copies to the
        // bin, depending on how many doublings come after it. At 3 copies it
        // gives 1 (adds 2 or 1).
        // Simplification: a gain that does not double the haplotype is left
        // out, although it copies a duplication made before a WGD: 1 -> 2 ->
        // 4 -> 6 on a chromosome at 3 is one event, counted as two (+2 +1).
        // See COPIED_BY_CHR_GAIN_INPUT.
        auto doublings = [w, normal](int c){ return w + (normal + c == 2 * normal ? 1 : 0); };

        const bool is_total = observed.size() == 1;   // {total} or {cnA, cnB}
        if(is_total){
            const int chr_cn = 2 * normal + cA + cB;
            // Chromosome events are doubled only by WGDs, hence w. A total does
            // not show which haplotype a segment event was on, so allow the
            // larger number of doublings of the two.
            return Expected{w, {fewest_events(cA + cB, w)},
                            {chr_cn == 0 ? 0 : fewest_events(observed[0] - chr_cn, std::max(doublings(cA), doublings(cB)))}};
        }
        auto seg = [&](int cn, int c){
            const int chr_cn = normal + c;
            return chr_cn == 0 ? 0 : fewest_events(cn - chr_cn, doublings(c));
        };
        // Chromosome events are doubled only by WGDs, hence w.
        return Expected{w, {fewest_events(cA, w), fewest_events(cB, w)},
                        {seg(observed[0], cA), seg(observed[1], cB)}};
    };
}

Expectation without_wgd(const std::vector<ChrChange>& chr_changes){
    return expected_changes(NO_WGD, chr_changes);
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

// Each region is commented with its expected change and one history of events
// from 1/1 that gives it, with the fewest events, all increases or all
// decreases on each haplotype at each level (haplotype A/B):
//   -W->            a WGD
//   -chr gain A->   a chromosome-level gain (-chr loss A->: loss) on haplotype A
//   -dup A->        a segment-level duplication (-del A->: deletion) on haplotype A

// Segment-level changes only: focal changes of one or two bins on otherwise
// normal chromosomes.
const std::vector<Region> SEGMENT_ONLY_INPUT{
    // duplication: seg +1/0
    // 1/1 -dup A-> 2/1
    {1, 1, 3, 4, 2, 1},
    // deletion: seg -1/0
    // 1/1 -del A-> 0/1
    {1, 2, 5, 5, 0, 1},
    // amplification: seg +2/0
    // 1/1 -dup A-> 2/1 -dup A-> 3/1
    {1, 3, 7, 7, 3, 1},
    // homozygous deletion: seg -1/-1
    // 1/1 -del A-> 0/1 -del B-> 0/0
    {1, 4, 2, 2, 0, 0},
    // copy-neutral LOH: seg +1/-1, total unchanged
    // 1/1 -dup A-> 2/1 -del B-> 2/0
    {2, 1, 8, 8, 2, 0},
    // duplication: seg 0/+1
    // 1/1 -dup B-> 1/2
    {2, 2, 1, 2, 1, 2},
    // amplification: seg 0/+2
    // 1/1 -dup B-> 1/2 -dup B-> 1/3
    {2, 3, 4, 4, 1, 3},
    // duplication of both haplotypes: seg +1/+1
    // 1/1 -dup A-> 2/1 -dup B-> 2/2
    {2, 4, 9, 9, 2, 2},
};

// Chromosome- and segment-level changes together, each segment with a single
// expected decomposition. Changes are small enough that no chromosome mean, and
// no sample ploidy, rounds the wrong way.
const std::vector<ChrChange> CHR_CHANGES{
    {1, 1, 1, 0},    // chr 1 gained on A
    {1, 2, -1, 0},   // chr 2 lost on A
    {2, 3, 0, 1},    // chr 3 gained on B
    {2, 4, 0, -1},   // chr 4 lost on B
};
const std::vector<Region> CHR_AND_SEGMENT_INPUT{
    // 1/1 -chr gain A-> 2/1
    whole_chr(1, 1, 2, 1),
    // 1/1 -chr loss A-> 0/1
    whole_chr(1, 2, 0, 1),
    // 1/1 -chr gain B-> 1/2
    whole_chr(2, 3, 1, 2),
    // 1/1 -chr loss B-> 1/0
    whole_chr(2, 4, 1, 0),
    // duplication after the chromosome gain: seg +1/0
    // 1/1 -chr gain A-> 2/1 -dup A-> 3/1
    {1, 1, 3, 3, 3, 1},
    // back at normal: seg -1/0, cancelling the chromosome gain
    // 1/1 -chr gain A-> 2/1 -del A-> 1/1
    {1, 1, 10, 10, 1, 1},
    // duplication on an unchanged chromosome: seg +1/0 (keeps bin 10 from being
    // dropped as the same in every sample)
    // 1/1 -dup A-> 2/1
    {2, 1, 10, 10, 2, 1},
    // duplication on the remaining copy: seg 0/+1
    // 1/1 -chr loss A-> 0/1 -dup B-> 0/2
    {1, 2, 5, 5, 0, 2},
    // duplication on A, not on B, whose chromosome was gained: seg +1/0
    // 1/1 -chr gain B-> 1/2 -dup A-> 2/2
    {2, 3, 7, 7, 2, 2},
    // duplication on the remaining copy: seg +1/0
    // 1/1 -chr loss B-> 1/0 -dup A-> 2/0
    {2, 4, 9, 9, 2, 0},
};

// num_wgd WGDs in both tumours: every chromosome starts at 2^num_wgd copies
// per haplotype (2/2 after one, 4/4 after two), then the changes. Changes keep
// each haplotype's sample mean at that level, so the ploidy (about 4 or 8)
// gives the WGD count and no chromosome mean rounds the wrong way.
std::vector<Region> after_wgd(int num_wgd, const std::vector<Region>& changes){
    const int normal = 1 << num_wgd;
    std::vector<Region> regions;
    for(int sample : {1, 2}){
        for(int chr = 1; chr <= NCHR; ++chr) regions.push_back(whole_chr(sample, chr, normal, normal));
    }
    regions.insert(regions.end(), changes.begin(), changes.end());
    return regions;
}

// Segment-level changes and WGD.
const std::vector<Region> SEGMENT_AND_WGD_INPUT = after_wgd(1, {
    // duplication after WGD: seg +1/0
    // 1/1 -W-> 2/2 -dup A-> 3/2
    {1, 1, 3, 4, 3, 2},
    // deletion after WGD: seg -1/0
    // 1/1 -W-> 2/2 -del A-> 1/2
    {1, 2, 5, 5, 1, 2},
    // duplication before WGD: seg +1/0, not two duplications after it
    // 1/1 -dup A-> 2/1 -W-> 4/2
    {1, 3, 7, 7, 4, 2},
    // deletion before WGD: seg -1/0
    // 1/1 -del A-> 0/1 -W-> 0/2
    {1, 4, 2, 2, 0, 2},
    // duplication before WGD: seg 0/+1
    // 1/1 -dup B-> 1/2 -W-> 2/4
    {2, 1, 8, 8, 2, 4},
    // duplication after WGD: seg 0/+1
    // 1/1 -W-> 2/2 -dup B-> 2/3
    {2, 2, 1, 2, 2, 3},
    // deletion after WGD: seg 0/-1
    // 1/1 -W-> 2/2 -del B-> 2/1
    {2, 3, 4, 4, 2, 1},
    // amplification before WGD: seg +2/0, not four duplications after it
    // 1/1 -dup A-> 2/1 -dup A-> 3/1 -W-> 6/2
    {2, 4, 9, 9, 6, 2},
    // homozygous deletion before WGD: seg -1/-1
    // 1/1 -del A-> 0/1 -del B-> 0/0 -W-> 0/0
    {2, 4, 2, 2, 0, 0},
});

// Chromosome-level changes and WGD, given as observed copies per haplotype.
const std::vector<ChrChange> CHR_AND_WGD_CHANGES{
    {1, 1, 2, 0},    // chr 1 gained on A before WGD: chr +1/0
    {1, 2, -2, 0},   // chr 2 lost on A before WGD: chr -1/0
    {2, 3, 0, 1},    // chr 3 gained on B after WGD: chr 0/+1
    {2, 4, 0, -1},   // chr 4 lost on B after WGD: chr 0/-1
};
const std::vector<Region> CHR_AND_WGD_INPUT = after_wgd(1, {
    // 1/1 -chr gain A-> 2/1 -W-> 4/2
    whole_chr(1, 1, 4, 2),
    // 1/1 -chr loss A-> 0/1 -W-> 0/2
    whole_chr(1, 2, 0, 2),
    // 1/1 -W-> 2/2 -chr gain B-> 2/3
    whole_chr(2, 3, 2, 3),
    // 1/1 -W-> 2/2 -chr loss B-> 2/1
    whole_chr(2, 4, 2, 1),
});

// Segment- and chromosome-level changes and WGD, on the chromosomes above.
std::vector<Region> with_chr_and_wgd(const std::vector<Region>& changes){
    std::vector<Region> regions = CHR_AND_WGD_INPUT;
    regions.insert(regions.end(), changes.begin(), changes.end());
    return regions;
}
const std::vector<Region> SEGMENT_CHR_AND_WGD_INPUT = with_chr_and_wgd({
    // duplication after WGD on a chromosome gained before it: seg +1/0
    // 1/1 -chr gain A-> 2/1 -W-> 4/2 -dup A-> 5/2
    {1, 1, 3, 3, 5, 2},
    // duplication between the chromosome gain and WGD: seg +1/0
    // 1/1 -chr gain A-> 2/1 -dup A-> 3/1 -W-> 6/2
    {1, 1, 5, 5, 6, 2},
    // back at 2/2: seg -1/0, cancelling the chromosome gain
    // 1/1 -chr gain A-> 2/1 -del A-> 1/1 -W-> 2/2
    {1, 1, 10, 10, 2, 2},
    // duplication after WGD on the remaining copy: seg 0/+1
    // 1/1 -chr loss A-> 0/1 -W-> 0/2 -dup B-> 0/3
    {1, 2, 5, 5, 0, 3},
    // duplication after WGD on A, not on B, whose chromosome was gained: seg +1/0
    // 1/1 -W-> 2/2 -chr gain B-> 2/3 -dup A-> 3/3
    {2, 3, 7, 7, 3, 3},
    // duplication on one copy of B after WGD, then copied by the chromosome
    // gain: seg 0/+1. The gain copies one of the two B copies (chromosome at
    // 3, not 4), so it adds 2, not 3.
    // 1/1 -W-> 2/2 -dup B-> 2/3 -chr gain B-> 2/5
    {2, 3, 2, 2, 2, 5},
    // duplication before WGD on A, not on B, whose chromosome was lost: seg +1/0
    // 1/1 -dup A-> 2/1 -W-> 4/2 -chr loss B-> 4/1
    {2, 4, 9, 9, 4, 1},
});

// A duplication before WGD that a later chromosome gain also copies, on the
// chromosomes above. It is one event, but counted as two under the
// simplification in doublings, so the expected value is seg 0/+2.
const std::vector<Region> COPIED_BY_CHR_GAIN_INPUT = with_chr_and_wgd({
    // duplication before WGD, copied by the WGD and then by the chromosome
    // gain (chromosome at 3): seg 0/+2, two events (+2 +1) here
    // 1/1 -dup B-> 1/2 -W-> 2/4 -chr gain B-> 2/6
    {2, 3, 5, 5, 2, 6},
});

// Segment changes made of events before and after a doubling (chromosome gain
// or WGD). Sample 1 has no WGD; sample 2 has one.
const std::vector<int> MIXED_TIMING_WGD{0, 1, 0};
const std::vector<ChrChange> MIXED_TIMING_CHR_CHANGES{
    {1, 1, 1, 0},    // chr 1 gained on A: chr +1/0
    {1, 2, -1, 0},   // chr 2 lost on A: chr -1/0
    {2, 1, 2, 0},    // chr 1 gained on A before WGD: chr +1/0
    {2, 2, -2, 0},   // chr 2 lost on A before WGD: chr -1/0
};
const std::vector<Region> MIXED_TIMING_INPUT = [](){
    std::vector<Region> regions;
    for(int chr = 1; chr <= NCHR; ++chr) regions.push_back(whole_chr(2, chr, 2, 2));
    const std::vector<Region> changes{
        // 1/1 -chr gain A-> 2/1
        whole_chr(1, 1, 2, 1),
        // 1/1 -chr loss A-> 0/1
        whole_chr(1, 2, 0, 1),
        // 1/1 -chr gain A-> 2/1 -W-> 4/2
        whole_chr(2, 1, 4, 2),
        // 1/1 -chr loss A-> 0/1 -W-> 0/2
        whole_chr(2, 2, 0, 2),
        // duplication before the chromosome gain: seg +1/0, not two duplications after it
        // 1/1 -dup A-> 2/1 -chr gain A-> 4/1
        {1, 1, 3, 3, 4, 1},
        // duplications before and after the chromosome gain: seg +2/0, not three after it
        // 1/1 -dup A-> 2/1 -chr gain A-> 4/1 -dup A-> 5/1
        {1, 1, 5, 5, 5, 1},
        // duplication between the chromosome gain and WGD: seg +1/0
        // 1/1 -chr gain A-> 2/1 -dup A-> 3/1 -W-> 6/2
        {2, 1, 3, 3, 6, 2},
        // duplications before and after the WGD: seg +2/0, not three after it
        // 1/1 -chr gain A-> 2/1 -dup A-> 3/1 -W-> 6/2 -dup A-> 7/2
        {2, 1, 8, 8, 7, 2},
        // duplications before and after the WGD: seg +2/0, not three after it
        // 1/1 -dup A-> 2/1 -W-> 4/2 -dup A-> 5/2
        {2, 3, 4, 4, 5, 2},
    };
    regions.insert(regions.end(), changes.begin(), changes.end());
    return regions;
}();

// Two WGDs in both tumours: every chromosome starts at 4/4 (total 8). A
// segment event adds 4 copies before both WGDs, 2 between them and 1 after
// them. With totals at most CN_MAX (10), only deletions can come before both
// WGDs.
const std::vector<int> TWO_WGD{2, 2, 0};

// Segment-level changes and two WGDs.
const std::vector<Region> SEGMENT_AND_TWO_WGD_INPUT = after_wgd(2, {
    // duplication after both WGDs: seg +1/0
    // 1/1 -W-> 2/2 -W-> 4/4 -dup A-> 5/4
    {1, 1, 3, 4, 5, 4},
    // duplication between the WGDs: seg +1/0, not two duplications after both
    // 1/1 -W-> 2/2 -dup A-> 3/2 -W-> 6/4
    {1, 2, 5, 5, 6, 4},
    // deletion between the WGDs: seg -1/0
    // 1/1 -W-> 2/2 -del A-> 1/2 -W-> 2/4
    {1, 3, 7, 7, 2, 4},
    // deletion after both WGDs: seg -1/0
    // 1/1 -W-> 2/2 -W-> 4/4 -del A-> 3/4
    {1, 4, 2, 2, 3, 4},
    // duplication between the WGDs: seg 0/+1
    // 1/1 -W-> 2/2 -dup B-> 2/3 -W-> 4/6
    {2, 1, 8, 8, 4, 6},
    // deletion before both WGDs: seg 0/-1; B is gone, so nothing later changes it
    // 1/1 -del B-> 1/0 -W-> 2/0 -W-> 4/0
    {2, 2, 1, 1, 4, 0},
    // deletions between and after the WGDs: seg -2/0, not three deletions after both
    // 1/1 -W-> 2/2 -del A-> 1/2 -W-> 2/4 -del A-> 1/4
    {2, 3, 4, 4, 1, 4},
    // duplication after both WGDs: seg 0/+1
    // 1/1 -W-> 2/2 -W-> 4/4 -dup B-> 4/5
    {2, 4, 9, 9, 4, 5},
});

// Chromosome-level changes and two WGDs, given as observed copies per haplotype.
const std::vector<ChrChange> CHR_AND_TWO_WGD_CHANGES{
    {1, 1, 2, 0},    // chr 1 gained on A between the WGDs: chr +1/0
    {1, 2, -2, 0},   // chr 2 lost on A between the WGDs: chr -1/0
    {1, 3, 1, 0},    // chr 3 gained on A after both WGDs: chr +1/0
    {1, 4, -1, 0},   // chr 4 lost on A after both WGDs: chr -1/0
    {2, 1, 2, 0},    // chr 1 gained on A between the WGDs: chr +1/0
    {2, 2, -4, 0},   // chr 2 lost on A before both WGDs: chr -1/0
    {2, 3, 2, 0},    // chr 3 gained on A between the WGDs: chr +1/0
};
const std::vector<Region> CHR_AND_TWO_WGD_INPUT = after_wgd(2, {
    // 1/1 -W-> 2/2 -chr gain A-> 3/2 -W-> 6/4
    whole_chr(1, 1, 6, 4),
    // 1/1 -W-> 2/2 -chr loss A-> 1/2 -W-> 2/4
    whole_chr(1, 2, 2, 4),
    // 1/1 -W-> 2/2 -W-> 4/4 -chr gain A-> 5/4
    whole_chr(1, 3, 5, 4),
    // 1/1 -W-> 2/2 -W-> 4/4 -chr loss A-> 3/4
    whole_chr(1, 4, 3, 4),
    // 1/1 -W-> 2/2 -chr gain A-> 3/2 -W-> 6/4
    whole_chr(2, 1, 6, 4),
    // 1/1 -chr loss A-> 0/1 -W-> 0/2 -W-> 0/4
    whole_chr(2, 2, 0, 4),
    // 1/1 -W-> 2/2 -chr gain A-> 3/2 -W-> 6/4
    whole_chr(2, 3, 6, 4),
});

// Segment- and chromosome-level changes and two WGDs, on the chromosomes above.
const std::vector<Region> SEGMENT_CHR_AND_TWO_WGD_INPUT = [](){
    std::vector<Region> regions = CHR_AND_TWO_WGD_INPUT;
    const std::vector<Region> changes{
        // deletion after both WGDs on a chromosome at 6/4: seg -1/0
        // 1/1 -W-> 2/2 -chr gain A-> 3/2 -W-> 6/4 -del A-> 5/4
        {1, 1, 3, 3, 5, 4},
        // back at 4/4, a deletion between the WGDs: seg -1/0
        // 1/1 -W-> 2/2 -chr gain A-> 3/2 -del A-> 2/2 -W-> 4/4
        {1, 1, 5, 5, 4, 4},
        // deletion after both WGDs on B, not on A, whose chromosome was gained: seg 0/-1
        // 1/1 -W-> 2/2 -chr gain A-> 3/2 -W-> 6/4 -del B-> 6/3
        {1, 1, 10, 10, 6, 3},
        // duplication after both WGDs on a chromosome at 2/4: seg +1/0
        // 1/1 -W-> 2/2 -chr loss A-> 1/2 -W-> 2/4 -dup A-> 3/4
        {1, 2, 5, 5, 3, 4},
        // deletion between the WGDs on a chromosome at 2/4: seg -1/0
        // 1/1 -W-> 2/2 -chr loss A-> 1/2 -del A-> 0/2 -W-> 0/4
        {1, 2, 7, 7, 0, 4},
        // duplication after both WGDs on a chromosome at 5/4: seg +1/0
        // 1/1 -W-> 2/2 -W-> 4/4 -chr gain A-> 5/4 -dup A-> 6/4
        {1, 3, 7, 7, 6, 4},
        // duplication after both WGDs on B, not on A, whose chromosome was lost: seg 0/+1
        // 1/1 -W-> 2/2 -W-> 4/4 -chr loss A-> 3/4 -dup B-> 3/5
        {1, 4, 9, 9, 3, 5},
        // deletion between the WGDs on B, not on A, whose chromosome was gained: seg 0/-1
        // 1/1 -W-> 2/2 -chr gain A-> 3/2 -del B-> 3/1 -W-> 6/2
        {2, 1, 2, 2, 6, 2},
        // duplication after both WGDs on the remaining copy: seg 0/+1
        // 1/1 -chr loss A-> 0/1 -W-> 0/2 -W-> 0/4 -dup B-> 0/5
        {2, 2, 5, 5, 0, 5},
        // duplication between the WGDs on the remaining copy: seg 0/+1
        // 1/1 -chr loss A-> 0/1 -W-> 0/2 -dup B-> 0/3 -W-> 0/6
        {2, 2, 6, 6, 0, 6},
        // deletion after both WGDs on B, not on A, whose chromosome was gained: seg 0/-1
        // 1/1 -W-> 2/2 -chr gain A-> 3/2 -W-> 6/4 -del B-> 6/3
        {2, 3, 7, 7, 6, 3},
    };
    regions.insert(regions.end(), changes.begin(), changes.end());
    return regions;
}();

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

TEST_CASE("model 3 conversion recovers segment-level changes and WGD from total copy numbers", "[parse_cn][decomp][wgd]"){
    check_conversion(write_input("segment-wgd", SEGMENT_AND_WGD_INPUT, true), true, expected_changes(ONE_WGD, {}));
}

TEST_CASE("model 3 conversion recovers segment-level changes and WGD from haplotype-specific copy numbers", "[parse_cn][decomp][wgd]"){
    check_conversion(write_input("segment-wgd", SEGMENT_AND_WGD_INPUT, false), false, expected_changes(ONE_WGD, {}));
}

TEST_CASE("model 3 conversion recovers chromosome-level changes and WGD from total copy numbers", "[parse_cn][decomp][wgd]"){
    check_conversion(write_input("chr-wgd", CHR_AND_WGD_INPUT, true), true, expected_changes(ONE_WGD, CHR_AND_WGD_CHANGES));
}

TEST_CASE("model 3 conversion recovers chromosome-level changes and WGD from haplotype-specific copy numbers", "[parse_cn][decomp][wgd]"){
    check_conversion(write_input("chr-wgd", CHR_AND_WGD_INPUT, false), false, expected_changes(ONE_WGD, CHR_AND_WGD_CHANGES));
}

// Known to fail (#50): the segment change is halved for the WGD before the chromosome change is subtracted,
// so a segment whose change and chromosome change differ in parity is wrong, e.g. 7 on a chromosome at 6 gives +2, not +1.
// Remove [!shouldfail] once fixed; Catch2 then reports this test if it passes.
TEST_CASE("model 3 conversion recovers segment- and chromosome-level changes and WGD from total copy numbers", "[parse_cn][decomp][wgd][!shouldfail]"){
    check_conversion(write_input("segment-chr-wgd", SEGMENT_CHR_AND_WGD_INPUT, true), true, expected_changes(ONE_WGD, CHR_AND_WGD_CHANGES));
}

TEST_CASE("model 3 conversion recovers segment- and chromosome-level changes and WGD from haplotype-specific copy numbers", "[parse_cn][decomp][wgd]"){
    check_conversion(write_input("segment-chr-wgd", SEGMENT_CHR_AND_WGD_INPUT, false), false, expected_changes(ONE_WGD, CHR_AND_WGD_CHANGES));
}

// Known to fail (#50): the expected value follows the simplification in
// doublings (seg +2), while the current code gives +1.
// Remove [!shouldfail] once fixed; Catch2 then reports these tests if they pass.
TEST_CASE("model 3 conversion counts a duplication copied by a later chromosome gain as two events from total copy numbers", "[parse_cn][decomp][wgd][!shouldfail]"){
    check_conversion(write_input("copied-by-chr-gain", COPIED_BY_CHR_GAIN_INPUT, true), true, expected_changes(ONE_WGD, CHR_AND_WGD_CHANGES));
}

TEST_CASE("model 3 conversion counts a duplication copied by a later chromosome gain as two events from haplotype-specific copy numbers", "[parse_cn][decomp][wgd][!shouldfail]"){
    check_conversion(write_input("copied-by-chr-gain", COPIED_BY_CHR_GAIN_INPUT, false), false, expected_changes(ONE_WGD, CHR_AND_WGD_CHANGES));
}

// Known to fail (#50): the code takes every segment event as coming after the
// doublings, so it counts too many events, e.g. 6 on a chromosome at 3 (5/1 on
// 2/1) gives +3, not +2. With haplotype-specific copy numbers, clamping the
// change to +-2 hides some of this, but 5/1 on 2/1 gives +1, not +2.
// Remove [!shouldfail] once fixed; Catch2 then reports these tests if they pass.
TEST_CASE("model 3 conversion recovers segment changes with events before and after a doubling from total copy numbers", "[parse_cn][decomp][wgd][!shouldfail]"){
    check_conversion(write_input("mixed-timing", MIXED_TIMING_INPUT, true), true, expected_changes(MIXED_TIMING_WGD, MIXED_TIMING_CHR_CHANGES));
}

TEST_CASE("model 3 conversion recovers segment changes with events before and after a doubling from haplotype-specific copy numbers", "[parse_cn][decomp][wgd][!shouldfail]"){
    check_conversion(write_input("mixed-timing", MIXED_TIMING_INPUT, false), false, expected_changes(MIXED_TIMING_WGD, MIXED_TIMING_CHR_CHANGES));
}

// Known to fail (#50): the code halves a change for WGD only once, whatever the
// WGD count, so after two WGDs a change before both is counted twice, e.g. 4/0
// on a chromosome at 4/4 gives segment -2 from total copy numbers, not -1, and
// a chromosome at 0/4 gives chromosome -2, not -1. The parity bug above also
// shows up here. Remove [!shouldfail] once fixed.
TEST_CASE("model 3 conversion recovers segment-level changes and two WGDs from total copy numbers", "[parse_cn][decomp][wgd][!shouldfail]"){
    check_conversion(write_input("segment-two-wgd", SEGMENT_AND_TWO_WGD_INPUT, true), true, expected_changes(TWO_WGD, {}));
}

TEST_CASE("model 3 conversion recovers segment-level changes and two WGDs from haplotype-specific copy numbers", "[parse_cn][decomp][wgd][!shouldfail]"){
    check_conversion(write_input("segment-two-wgd", SEGMENT_AND_TWO_WGD_INPUT, false), false, expected_changes(TWO_WGD, {}));
}

TEST_CASE("model 3 conversion recovers chromosome-level changes and two WGDs from total copy numbers", "[parse_cn][decomp][wgd][!shouldfail]"){
    check_conversion(write_input("chr-two-wgd", CHR_AND_TWO_WGD_INPUT, true), true, expected_changes(TWO_WGD, CHR_AND_TWO_WGD_CHANGES));
}

TEST_CASE("model 3 conversion recovers chromosome-level changes and two WGDs from haplotype-specific copy numbers", "[parse_cn][decomp][wgd]"){
    check_conversion(write_input("chr-two-wgd", CHR_AND_TWO_WGD_INPUT, false), false, expected_changes(TWO_WGD, CHR_AND_TWO_WGD_CHANGES));
}

TEST_CASE("model 3 conversion recovers segment- and chromosome-level changes and two WGDs from total copy numbers", "[parse_cn][decomp][wgd][!shouldfail]"){
    check_conversion(write_input("segment-chr-two-wgd", SEGMENT_CHR_AND_TWO_WGD_INPUT, true), true, expected_changes(TWO_WGD, CHR_AND_TWO_WGD_CHANGES));
}

TEST_CASE("model 3 conversion recovers segment- and chromosome-level changes and two WGDs from haplotype-specific copy numbers", "[parse_cn][decomp][wgd][!shouldfail]"){
    check_conversion(write_input("segment-chr-two-wgd", SEGMENT_CHR_AND_TWO_WGD_INPUT, false), false, expected_changes(TWO_WGD, CHR_AND_TWO_WGD_CHANGES));
}
