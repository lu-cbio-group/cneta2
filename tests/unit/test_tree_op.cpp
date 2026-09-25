// Unit tests for the tree helpers in tree_op.cpp / evo_tree.cpp.
//
// The tree-string functions deserve the most attention: tree search uses them
// to decide whether two candidate topologies are the same, so a change in
// their behaviour silently changes which trees get explored.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include <fstream>
#include <string>
#include <vector>

#include "tree_op.hpp"

namespace {

std::string data_path(const std::string& name){
    return std::string(CNETA_TEST_DATA_DIR) + "/" + name;
}

const int NS = 3;        // tumour regions in the fixture
const int NLEAF = NS + 1;  // plus the normal sample

}  // namespace

TEST_CASE("get_tree_height returns the largest node time", "[tree_op]"){
    std::vector<double> node_times{0.0, 3.5, 1.25, 7.0, 2.0};

    REQUIRE_THAT(get_tree_height(node_times),
                 Catch::Matchers::WithinRel(7.0, 1e-12));
}

TEST_CASE("get_total_time discounts the last sampling time", "[tree_op]"){
    // Height measured to the last sample, minus how long after the first
    // sample that last one was taken.
    std::vector<double> node_times{0.0, 3.5, 7.0};

    REQUIRE_THAT(get_total_time(node_times, 2.0),
                 Catch::Matchers::WithinRel(5.0, 1e-12));
    REQUIRE_THAT(get_total_time(node_times, 0.0),
                 Catch::Matchers::WithinRel(7.0, 1e-12));
}

TEST_CASE("get_first_sample finds the earliest sampled region", "[tree_op]"){
    // tobs holds times relative to the first sample, so the first sample is
    // the entry equal to zero.
    REQUIRE(get_first_sample(std::vector<double>{0.0, 4.0, 4.0}) == 0);
    REQUIRE(get_first_sample(std::vector<double>{4.0, 0.0, 2.0}) == 1);

    // No zero entry means no first sample could be identified.
    REQUIRE(get_first_sample(std::vector<double>{1.0, 4.0, 2.0}) == -1);
}

TEST_CASE("read_tree_info loads the fixture tree", "[tree_op]"){
    evo_tree tree = read_tree_info(data_path("tiny-tree.txt"), NS);

    REQUIRE(tree.nleaf == NLEAF);
    // A rooted binary tree on n leaves has 2n-2 edges.
    REQUIRE(tree.edges.size() == static_cast<size_t>(2 * NLEAF - 2));
    // ... and 2n-1 nodes.
    REQUIRE(tree.nodes.size() == static_cast<size_t>(2 * NLEAF - 1));
}

TEST_CASE("is_blen_valid rejects negative branch lengths", "[tree_op]"){
    evo_tree tree = read_tree_info(data_path("tiny-tree.txt"), NS);

    SECTION("the fixture tree is valid as read"){
        REQUIRE(is_blen_valid(tree));
    }

    SECTION("flipping one branch negative makes it invalid"){
        tree.edges[0].length = -1.0;
        REQUIRE_FALSE(is_blen_valid(tree));
    }
}

TEST_CASE("edge lengths agree with node times", "[tree_op]"){
    // Every edge must span exactly the time between the nodes it connects.
    // Optimisation updates times and lengths separately, so this consistency
    // is worth stating explicitly.
    evo_tree tree = read_tree_info(data_path("tiny-tree.txt"), NS);
    std::vector<double> times = tree.get_node_times();

    for(const auto& e : tree.edges){
        const double span = times[e.end] - times[e.start];
        INFO("edge " << e.start << " -> " << e.end);
        REQUIRE_THAT(e.length, Catch::Matchers::WithinAbs(span, 1e-9));
    }
}

TEST_CASE("tree strings ignore the order edges are listed in", "[tree_op]"){
    // tiny-tree-shuffled.txt holds exactly the same topology and branch
    // lengths as tiny-tree.txt, with the edge rows in a different order.
    //
    // Note it is order_tree_string_uniq() that supplies the invariance:
    // create_tree_string_uniq() emits the root's two children in whatever
    // order they were stored, so the two must be used together. That pairing
    // is the contract tree search relies on.
    evo_tree tree = read_tree_info(data_path("tiny-tree.txt"), NS);
    evo_tree shuffled = read_tree_info(data_path("tiny-tree-shuffled.txt"), NS);

    REQUIRE(order_tree_string_uniq(create_tree_string_uniq(tree))
            == order_tree_string_uniq(create_tree_string_uniq(shuffled)));
}

TEST_CASE("different topologies give different tree strings", "[tree_op]"){
    // Without this, tree search could treat two genuinely distinct
    // topologies as already-visited and skip one.
    evo_tree tree = read_tree_info(data_path("tiny-tree.txt"), NS);
    evo_tree other = read_tree_info(data_path("tiny-tree-alt.txt"), NS);

    REQUIRE(order_tree_string_uniq(create_tree_string_uniq(tree))
            != order_tree_string_uniq(create_tree_string_uniq(other)));
}

// KNOWN FAILURE - tagged [!mayfail] so it reports without breaking the build.
//
// order_tree_string_uniq() splits on ':' and iterates to split1.size()-1,
// intending to drop the empty string that the trailing ':' produces. But it
// sorts first, which moves that empty string to the *front*, so the loop
// instead drops the lexicographically largest component - a real node. One
// node is therefore missing from every canonical tree string, and applying
// the function twice loses another.
//
// Two topologies that differ only in their largest component would hash to
// the same string and tree search would treat one as already visited. Remove
// the tag once the off-by-one is fixed; the assertion below should then pass
// as written.
TEST_CASE("order_tree_string_uniq is idempotent", "[tree_op][!mayfail]"){
    evo_tree tree = read_tree_info(data_path("tiny-tree.txt"), NS);

    const std::string once = order_tree_string_uniq(create_tree_string_uniq(tree));
    REQUIRE(order_tree_string_uniq(once) == once);
}

TEST_CASE("canonical tree strings contain no empty component", "[tree_op][!mayfail]"){
    // Same root cause as above, stated so that it is visible directly. The
    // sorted-empty-string-first bug shows up as a leading ':' -- an empty
    // component standing in for the node that was dropped off the end.
    evo_tree tree = read_tree_info(data_path("tiny-tree.txt"), NS);

    const std::string canonical =
        order_tree_string_uniq(create_tree_string_uniq(tree));

    INFO("canonical string: " << canonical);
    REQUIRE_FALSE(canonical.empty());
    REQUIRE(canonical.front() != ':');
    REQUIRE(canonical.find("::") == std::string::npos);
}
