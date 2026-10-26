#!/usr/bin/env bash

# This script is used to run program cnetmcmc.

# This script lives in <repo>/bin and finds the programs, config/ and util/
# relative to the repository, so it can be run from any directory. Output
# (CNETA_OUT) is still relative to the directory you run it from.
REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"

# Where the compiled executables live. The CMake build puts them in <repo>/bin,
# next to this script; override with CNETA_BIN to point at a different build.
CNETA_BIN="${CNETA_BIN:-$REPO_ROOT/bin}"

# Where the configuration files live; override with CNETA_CONFIG_DIR.
CNETA_CONFIG_DIR="${CNETA_CONFIG_DIR:-$REPO_ROOT/config}"

# Settings shared by all run scripts (output directory, seed). Anything already
# set in the environment takes precedence over this file.
if [[ -f "$CNETA_CONFIG_DIR/common.conf" ]]; then
  source "$CNETA_CONFIG_DIR/common.conf"
fi

# Use MCMC to infer tree
# With configuration file; override with CNETA_CONFIG to use a different one.
config_file="${CNETA_CONFIG:-$CNETA_CONFIG_DIR/cnet_mcmc.cfg}"
# cnetmcmc stops on a missing config file too; checking here gives the
# message before anything else runs.
if [[ ! -f $config_file ]]; then
  echo "Configuration file $config_file does not exist!" >&2
  exit 1
fi

# Seed for the random number generator. Empty or 0 (the default) means draw a
# fresh 32-bit seed for this run. Each chain's seed is written at the top of
# its log and appended to seed-history.txt in the output directory, so any run
# can be repeated with   CNETA_SEED=<that number> bin/run-cnetmcmc.sh
if [[ -z "${CNETA_SEED:-}" || "$CNETA_SEED" == "0" ]]; then
  seed=$(od -An -N4 -tu4 /dev/urandom | tr -d ' ')
  if [[ $seed -eq 0 ]]; then seed=1; fi
else
  seed=$CNETA_SEED
fi

# record_seed <output dir> <seed> <log file>: append one line to
# <output dir>/seed-history.txt. Re-running with the same settings overwrites
# the log, so this file is the record that survives.
record_seed () {
  printf '%s\t%s\tseed=%s\tlog=%s\n' "$(date '+%Y-%m-%dT%H:%M:%S')" "$(basename "$0")" "$2" "$3" >> "$1/seed-history.txt"
}

# derived_seed <seed> <n>: the seed for chain n + 1, seed + n, kept inside the
# 32-bit range and never 0 (which would mean "pick one").
derived_seed () {
  echo $(( ($1 + $2 - 1) % 4294967295 + 1 ))
}

####################### Parameters related to input  ###########################
Ns=3  # The number of regions

# input parameters
dir="${CNETA_OUT:-./example}"   # The output directory
# prefix=sim-data-N5-cons1-model1-1
prefix=sim-data-1
is_total=1 # If yes, the input is total copy number
input=$dir/"$prefix"-cn.txt.gz
# input=$dir/"$prefix"-allele-cn.txt.gz

# The input file of sample timings (optional, required for estimating mutation rates)
times=$dir/"$prefix"-rel-times.txt
if [[ ! -f $times ]]; then
  times=""
fi

# Optional reference real tree
# rtree=$dir/"$prefix"-tree.txt
rtree=""

init_tree=0     # Types of starting tree for MCMC. 0: random tree, 1: provided tree, 2: random tree with the same topology as real tree
file_itree=$dir/"$prefix"-tree.txt
################################################################################


####################### Set output directory and run the program ###############
# output parameters
odir=$dir/mcmc
if [[ ! -d $odir ]]; then
  mkdir -p $odir
fi
suffix=$prefix

echo "Start running cnetmcmc"

# Run multiple chains
for i in {1..2}
do
  # Each chain gets its own seed (seed, seed + 1, ...): with one shared seed
  # the chains would be identical copies rather than independent runs.
  chain_seed=$(derived_seed $seed $((i - 1)))
  echo "chain $i"
  # Use the suffix of MrBayes for compatibility with RWTY
  trace_param_file=$odir/mcmc_"$suffix"_${i}.p
  trace_tree_file=$odir/mcmc_"$suffix"_${i}.t
  sum_tree_file=$odir/sum-mcmc_"$suffix"_${i}.txt

  log=$odir/std_mcmc_"$suffix"_${i}
  echo "seed $chain_seed" > $log
  record_seed "$dir" "$chain_seed" "mcmc/$(basename "$log")"

  "$CNETA_BIN/cnetmcmc" -s $Ns --is_total $is_total -c $input -t "$times" --rtree "$rtree" --trace_param_file $trace_param_file --trace_tree_file $trace_tree_file --config_file "$config_file" --init_tree $init_tree --file_itree $file_itree --seed $chain_seed >> $log
  chain_status=$?
  if [[ $chain_status -ne 0 ]]; then
    echo "cnetmcmc chain $i FAILED (exit code ${chain_status}), see $log"
    exit "$chain_status"
  fi
  echo "cnetmcmc chain $i SUCCEEDED (exit code 0)"

  # Summarize the sampled trees into a maximum credibility tree with median heights
  # treeannotator -burnin 10 -heights median $trace_tree_file $sum_tree_file
done

echo "Finish running cnetmcmc"
