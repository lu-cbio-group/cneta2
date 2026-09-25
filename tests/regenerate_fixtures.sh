#!/usr/bin/env bash
#
# Regenerate the committed test fixtures in tests/data.
#
# You should rarely need this. The fixtures are committed precisely so that
# the regression tests do not depend on cnets continuing to produce the same
# simulation. Regenerate only when:
#
#   * the input file format changes, or
#   * the likelihood calculation changes deliberately and the reference
#     log-likelihoods in expected.json need to move with it.
#
# Either way, commit the new fixtures on their own, with a message saying
# what changed and why. A silent fixture update makes the regression tests
# worthless.
#
# Usage:
#   cd code && ./build.sh local 4     # fixtures come from the built binaries
#   tests/regenerate_fixtures.sh

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DATA_DIR="$REPO_ROOT/tests/data"
CNETA_BIN="${CNETA_BIN:-$REPO_ROOT/bin}"

# Single-threaded, so the reference values are reproducible.
export OMP_NUM_THREADS=1

SEED=12345
NS=3
SEG_MAX=50
CN_MAX=4

for prog in cnets cnetml; do
  if [[ ! -x "$CNETA_BIN/$prog" ]]; then
    echo "Error: $CNETA_BIN/$prog not found. Build first: cd code && ./build.sh local 4" >&2
    exit 1
  fi
done

WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

echo "Simulating with seed $SEED into $WORK"
"$CNETA_BIN/cnets" \
  -o "$WORK/" -p "" -r "$NS" -n 1 \
  --mode 1 --method 1 --fix_nseg 1 --seg_max "$SEG_MAX" --cn_max "$CN_MAX" \
  --model 2 --dup_rate 0.001 --del_rate 0.001 \
  --chr_gain 0 --chr_loss 0 --wgd 0 --dup_size 5 --del_size 5 \
  -e 9000000 -b 1.563e-3 --gtime 0.002739726 -t 2 \
  --age 60 --constrained 1 --print_relative 1 --verbose 0 --seed "$SEED" \
  > "$WORK/cnets.log"

echo "Copying fixtures into $DATA_DIR"
cp "$WORK/sim-data-1-cn.txt.gz"            "$DATA_DIR/tiny-cn.txt.gz"
cp "$WORK/sim-data-1-haplotype-cn.txt.gz"  "$DATA_DIR/tiny-haplotype-cn.txt.gz"
cp "$WORK/sim-data-1-tree.txt"             "$DATA_DIR/tiny-tree.txt"
cp "$WORK/sim-data-1-rel-times.txt"        "$DATA_DIR/tiny-rel-times.txt"
gunzip -c "$DATA_DIR/tiny-cn.txt.gz"     > "$DATA_DIR/tiny-cn.txt"

# tiny-tree-shuffled.txt and tiny-tree-alt.txt are hand-written variants of
# tiny-tree.txt (same topology with rows reordered, and a different topology
# respectively). If the tree above changes, update those two by hand as well.
echo "Note: tiny-tree-shuffled.txt and tiny-tree-alt.txt are hand-maintained." >&2

score_tree() {
  "$CNETA_BIN/cnetml" \
    -c "$DATA_DIR/tiny-cn.txt.gz" \
    -t "$DATA_DIR/tiny-rel-times.txt" \
    --tree_file "$1" \
    -s "$NS" --is_total 1 --is_bin 0 --incl_all 1 --m_max 1 \
    -d 2 --cn_max "$CN_MAX" --cn_type 0 --constrained 1 --correct_bias 1 \
    -x 0 --dup_rate 0.001 --del_rate 0.001 \
    --chr_gain_rate 0 --chr_loss_rate 0 --wgd_rate 0 \
    --mode 2 --verbose 0 --seed "$SEED" \
  | sed -n 's/.*log likelihood of the input tree is  *//p'
}

echo
echo "Reference log-likelihoods for tests/data/expected.json:"
for tree in tiny-tree.txt tiny-tree-alt.txt; do
  printf '    "%s": %s\n' "$tree" "$(score_tree "$DATA_DIR/$tree")"
done
echo
echo "Update the mode2_logl block in tests/data/expected.json with the values"
echo "above, along with the _generated metadata."
