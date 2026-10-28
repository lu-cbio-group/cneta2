#!/usr/bin/env bash

# This script is used to run program cnetml.

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

# Seed for the random number generator. Empty or 0 (the default) means draw a
# fresh 32-bit seed for this run. The seed in use is written at the top of each
# log and appended to seed-history.txt in the output directory, so any run can
# be repeated with   CNETA_SEED=<that number> bin/run-cnetml.sh
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

# derived_seed <seed> <n>: the seed for the n-th bootstrap replicate, seed + n,
# kept inside the 32-bit range and never 0 (which would mean "pick one").
derived_seed () {
  echo $(( ($1 + $2 - 1) % 4294967295 + 1 ))
}

verbose=0  # Whether or not to print debug information. 0: default, 1: standard debug, 2: debug with details on tree

# Running mode. 0: build maximum likelihood tree; 1: test; 2: compute likelihood; 3: compute maximum likelihood; 4: infer ancestral state
mode=0
plot=0 # Whether or not to plot the inferred trees

# Branch-specific rate mode. 0: constant rate; 1: shared branch multiplier; 2: independent per-type; 3: ML-RLC
bsr_mode=0
# Search method for bsr_mode=3 (ML-RLC) shift-edge selection. 0: exhaustive search (all shift-edge subsets, small trees only); 1: stepwise greedy search (default); 2: genetic algorithm (not implemented yet)
rlc_search_method=1
# Model-selection criterion for bsr_mode=3 shift-edge search. 0: AIC (default), 1: BIC
rlc_criterion=0
# bsr_mode>0 only. How the reference rate for the branch multipliers is set.
# 0 (default): skip calibration; the trunk edge (root -> MRCA) keeps multiplier 1 and is never a shift edge, its own rate (mu0)
#    is estimated jointly with the other edges' multipliers/shifts, which are relative to it.
# 1: calibrate bsr0 rates first (constant-rate fit), then estimate per-edge multipliers relative to that reference.
estimate_bsr0_first=0

####################### Parameters related to input  ###########################
idir="${CNETA_OUT:-./example}"   # The input directory
prefix=sim-data-1   # The prefix of input files
# prefix=sim-data-"$1"

# The input parameters related to copy numbers
input=$idir/"$prefix"-cn.txt.gz
# input=$idir/"$prefix"-allele-cn.txt.gz
if [[ ! -f $input ]]; then
  echo "Input file $input does not exit!"
  exit 0
fi
seg_file=$idir/"$prefix"-segs.txt

# Whether or not the input copy number is relative copy number (to normal)
is_rcn=0

is_total=1    # If yes (1), the input is total copy number; or else (0), the input is haplotype-specific
# input=$idir/"$prefix"-allele-cn.txt.gz
# whether or not the input copy number is for each bin. If not, the input copy number is read as it is. Or else, consecutive bins will be merged
is_bin=0
# whether or not to include all the input copy numbers without further propressing
incl_all=1
cn_max=4  # Maximum copy number allowed

# The input parameters that should be consistent with the simulation or real data
Ns=3 # The number of samples for a patient
cons=1  # Whether or not the tree is constrained by patient age

# The input file of sample timings (optional, required for estimating mutation rates)
times=$idir/"$prefix"-rel-times.txt
if [[ ! -f $times ]]; then
  times=""
fi

model=2  # Model of evolution.  0: Mk model, 1: bounded model for total copy number, 2: bounded model for haplotype-specific copy number, 3: independent Markov chain model
# Parameters for independent Markov chain model
m_max=1

estmu=1  # Whether or not to estimate mutation rate
cn_type=0  # 0-ONLY_SEG: only use segments level CN changes (use only r1/r2); 1-EXCLUDE_SEG: exclude segment level CN changes (use only r3/r4/r5); 2-EXCLUDE_CHR: exclude chromosome level changes (use only r1/r2/r5) ; 3-EXCLUDE_WGD: exclude whole genome doubling (use only r1/r2/r3/r4); 4-ALL: three types of mutations (use r1 to r5)

if [[ $times == "" ]]; then
  estmu=0
  cons=0
fi

# With bsr_mode>0 the rates are fixed as the reference for the branch multipliers, so cnetml requires estmu=0.
# estmu_mode0 keeps the value before this guard; mode 3 uses it to find the tree built in mode 0.
estmu_mode0=$estmu
if [[ $bsr_mode -gt 0 && $estmu -eq 1 ]]; then
  echo "Warning: estmu=1 is incompatible with bsr_mode>0, setting estmu=0"
  estmu=0
fi
################################################################################


######### Parameters for optimization and tree search algorithms  ##############
correct_bias=1  # Whether or not to correct acquisition bias in likelihood computation
opt=1   # 1: L-BFGS-B method; 0: simplex method
tolerance=1e-2
miter=2000

# Bounds on the number of events per haplotype; they set the dimension of the rate matrices.
# Values below are the program defaults. They must be >= the observed max for the data:
# for simulated data, take the column max of simulation_cn_summary.txt in the cnets output
# directory (columns: prefix, max_wgd, max_chr_change, max_site_change).
max_chr_change_haplotype=1
max_site_change_haplotype=2

tree_search=2  # 0: genetic algorithm, 1: hill climbing, 2: exhaustive search

# Guard: cnetml supports bsr_mode>0 only with mode=3 (fixed topology). Build the tree with
# bsr_mode=0 in mode 0 first, then run mode 3 with bsr_mode>0 on that tree.
if [[ $bsr_mode -gt 0 && $mode -ne 3 ]]; then
  echo "Error: bsr_mode=$bsr_mode requires mode=3 (got mode=$mode)."
  echo "Build the tree with bsr_mode=0 and mode=0 first, then rerun with mode=3."
  exit 1
fi
# Parameters used in genetic algorithm
Npop=100  # size of initial trees for hill climbing
Ngen=20
Nstop=5
# whether or not to do reduced NNI while doing hill climbing NNIs
speed_nni=1

# Parameters about starting tree during search
init_tree=0  # 0: Random coalescence tree, 1: Maximum parsimony tree
dir_itrees=$idir/itrees  # The directory containing the initial trees (in TXT format)
tree_file=$idir/"$prefix"-tree.txt   # Only required when the tree is to be given as input
# tree_file=""
# effective population size for random initial coalescence tree
Ne=90000
beta=1.563e-3  # exponential growth rate
gtime=0.002739726 # generation time in year


# Specify the (initial) mutation rates when they are fixed
r1=0.001
r2=0.001
r3=0
r4=0
r5=0
mu=0  # Used in model 0

# When rates are estimated (estmu=1) or used as the reference for branch multipliers (bsr_mode>0),
# every rate used by cn_type must be non-zero. A used rate left at 0 above is set to 0.0005.
# cn_type=0: r1 r2; 1: r3 r4 r5; 2: r1 r2 r5; 3: r1 r2 r3 r4; 4: r1-r5
if [[ $estmu -eq 1 || $bsr_mode -gt 0 ]]; then
  case $cn_type in
    0) used_rates="r1 r2" ;;
    1) used_rates="r3 r4 r5" ;;
    2) used_rates="r1 r2 r5" ;;
    3) used_rates="r1 r2 r3 r4" ;;
    4) used_rates="r1 r2 r3 r4 r5" ;;
    *) echo "Invalid value for cn_type: $cn_type"; exit 1 ;;
  esac
  for r in $used_rates; do
    if awk -v x="${!r}" 'BEGIN { exit !(x == 0) }'; then
      echo "Warning: $r is used by cn_type=$cn_type but set to 0, using 0.0005 instead"
      printf -v "$r" '%s' 0.0005
    fi
  done
fi
################################################################################

# Options appended to every cnetml call below, so the main run and
# its bootstrap replicates always use the same settings
extra_opts="--is_rcn $is_rcn --max_chr_change_haplotype $max_chr_change_haplotype --max_site_change_haplotype $max_site_change_haplotype --bsr_mode $bsr_mode --rlc_search_method $rlc_search_method --rlc_criterion $rlc_criterion --estimate_bsr0_first $estimate_bsr0_first"


####################### Set output directory and run the program ###############
echo "Start running cnetml"

dir=$idir
if [[ ! -d $dir ]]; then
  mkdir -p $dir
fi
bsdir1=$dir/bootstrap
bsdir2=$dir/bootstrap_fixt

if [[ $mode -eq 0 ]]; then
  suffix=m$model-o"$opt"-s"$tree_search"-cons${cons}-estmu${estmu}-"$prefix"
  mltree=$dir/MaxL-"$suffix".txt

  log=$dir/std_cnetml_"$suffix"
  echo "seed $seed" > $log
  record_seed "$dir" "$seed" "$(basename "$log")"
  # valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes -v
  time "$CNETA_BIN/cnetml" -c $input --seg_file $seg_file -t "$times" --tree_file "$tree_file" --is_total $is_total --m_max $m_max --is_bin $is_bin --incl_all $incl_all -s $Ns -p $Npop -g $Ngen -e $Nstop -r $tolerance -o $mltree -d $model --cn_max $cn_max --cn_type $cn_type --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --estmu $estmu --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode --speed_nni $speed_nni $extra_opts --seed $seed 2>&1 >> $log
  cnetml_status=$?
  # #
  #
  if [[ $cnetml_status -ne 0 ]]; then
    echo "cnetml main run FAILED (exit code ${cnetml_status})"
  else
    echo "cnetml main run SUCCEEDED (exit code 0)"
  fi

  echo "Finish running cnetml"

  if [[ $cnetml_status -ne 0 ]]; then
    exit "$cnetml_status"
  fi

  if [[ $plot -eq 1 ]]; then
    Rscript "$REPO_ROOT/util/plot-trees-all.R" -f $mltree -b 0 -t "single" -l "xlim" --time_file "$times"  #>& /dev/null
  fi
  # Evaluate the estimation error
  # cmp_plot=$dir/cmp_plot-"$suffix".pdf
  # cmp_dist=$dir/cmp_dist-"$suffix".txt
  # Rscript test/tree_comparison/compare_trees.R -r "$tree_file" -i $mltree -o $cmp_plot -s $cmp_dist -p 1

  bootstrap=0
  if [[ $bootstrap -eq 1 ]]; then
    echo "Run bootstrapping"
    # Do bootstrapping
    mkdir -p $bsdir1
    for i in {1..10}
    do
      echo $i
      ofile=$bsdir1/MaxL-"$suffix"-btree-$i.txt
      # Each replicate resamples the sites at random, so each needs its own seed.
      bseed=$(derived_seed $seed $i)
      log=$bsdir1/std_cnetml_"$suffix"-btree-$i
      echo "seed $bseed" > $log
      record_seed "$dir" "$bseed" "bootstrap/$(basename "$log")"

      "$CNETA_BIN/cnetml" -c $input -t "$times" --tree_file "$tree_file" --is_total $is_total --m_max $m_max --is_bin $is_bin --incl_all $incl_all -s $Ns -p $Npop -g $Ngen -e $Nstop -r $tolerance -o $ofile -d $model --cn_max $cn_max --cn_type $cn_type --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --estmu $estmu  --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2  --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode -b 1 $extra_opts --seed $bseed >> $log
    done
    # Draw the ML tree with bootstrapping support
    if [[ $plot -eq 1 ]]; then
      Rscript "$REPO_ROOT/util/plot-trees-all.R" -s $bsdir1 -f $mltree --tree_file_nex "$mltree".nex -o $dir/MaxL-tree-"$suffix"-bootstrap.pdf -t "bootstrap" -l "age" --time_file "$times" -p "MaxL-"$suffix"-btree-*txt"
    fi
  fi

elif [[ $mode -eq 1 ]]; then
  suffix=sim1-m$model-test-"$prefix"
  mltree=$dir/MaxL-"$suffix".txt

  log=$dir/std_cnetml_"$suffix"
  echo "seed $seed" > $log
  record_seed "$dir" "$seed" "$(basename "$log")"
  "$CNETA_BIN/cnetml" -c $input -t "$times" --is_bin $is_bin --incl_all $incl_all --tree_file "$tree_file"  --is_total $is_total --m_max $m_max -s $Ns -p $Npop -g $Ngen -e $Nstop -r $tolerance -o $mltree -d $model --cn_max $cn_max --cn_type $cn_type --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --estmu $estmu --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose $extra_opts --seed $seed >> $log

elif [[ $mode -eq 2 ]]; then
  # In this mode, the mutation rates have to be specified
  suffix=m$model-"$cons""$estmu"-mode"$mode"-"$prefix"

  log=$dir/std_cnetml_"$suffix"
  echo "seed $seed" > $log
  record_seed "$dir" "$seed" "$(basename "$log")"
  "$CNETA_BIN/cnetml" -c $input -t "$times" --is_bin $is_bin --incl_all $incl_all -s $Ns --tree_file "$tree_file" --is_total $is_total --m_max $m_max -d $model --cn_max $cn_max --cn_type $cn_type --constrained $cons --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode  $extra_opts --seed $seed >> $log

elif [[ $mode -eq 3 ]]; then
  # In this mode, the mutation rates can be specified or not
  # Fixed topology to use:
  #   0: the tree built in mode 0 (always bsr_mode=0, named with the estmu set before the bsr guard); for real data
  #   1: the simulated true tree from cnets; for simulated data
  fixt_tree=0
  if [[ $fixt_tree -eq 1 ]]; then
    mltree=$idir/"$prefix"-tree.txt
  else
    mltree=$dir/MaxL-m$model-o"$opt"-s"$tree_search"-cons${cons}-estmu${estmu_mode0}-"$prefix".txt
  fi
  # Outputs carry the bsr settings so runs with different settings don't overwrite each other.
  suffix=m$model-o"$opt"-s"$tree_search"-cons${cons}-estmu${estmu}-bsr${bsr_mode}-rsm${rlc_search_method}-ric${rlc_criterion}-ebf${estimate_bsr0_first}-"$prefix"
  mltree2=$dir/MaxL-"$suffix"-fixt.txt
  echo "Estimating the branch lengths (and mutation rates) of the given tree $mltree"
  if [[ ! -f $mltree ]]; then
  echo "$mltree does not exist!"
  exit
  fi

  log=$dir/std_cnetml_"$suffix"-mode"$mode"
  echo "seed $seed" > $log
  record_seed "$dir" "$seed" "$(basename "$log")"
  # valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes -v
  time "$CNETA_BIN/cnetml" --tree_file "$mltree" -o $mltree2 -c $input -t "$times" --is_bin $is_bin --incl_all $incl_all -s $Ns --is_total $is_total --m_max $m_max -d $model --cn_max $cn_max --cn_type $cn_type --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --estmu $estmu --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode $extra_opts --seed $seed >> $log

  # plot the new tree
  if [[ $plot -eq 1 ]]; then
    Rscript "$REPO_ROOT/util/plot-trees-all.R" -f $mltree2 -b 0 -t "single" -l "xlim" --time_file "$times"  #>& /dev/null
  fi

  bootstrap=1
  if [[ $bootstrap -eq 1 ]]; then
    echo "Run bootstrapping"
    # Do bootstrapping
    mkdir -p $bsdir2

    for i in {1..10}
    do
      echo $i
      ofile=$bsdir2/MaxL-"$suffix"-btree-$i.txt
      # Each replicate resamples the sites at random, so each needs its own seed.
      bseed=$(derived_seed $seed $i)
      log=$bsdir2/std_cnetml_"$suffix"-btree-$i
      echo "seed $bseed" > $log
      record_seed "$dir" "$bseed" "bootstrap_fixt/$(basename "$log")"

      # valgrind --tool=memcheck --leak-check=full --show-leak-kinds=all --track-origins=yes -v /usr/bin/time
      "$CNETA_BIN/cnetml" --tree_file "$mltree" -o $ofile -c $input -t "$times" --is_bin $is_bin --incl_all $incl_all -s $Ns  --is_total $is_total --m_max $m_max -d $model --cn_max $cn_max --cn_type $cn_type --tree_search $tree_search --init_tree $init_tree --epop $Ne --beta $beta --gtime $gtime --dir_itrees $dir_itrees --optim $opt --constrained $cons --estmu $estmu --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode -b 1 $extra_opts --seed $bseed >> $log
    done

    # Collect the rates of all replicates from their summary files into one table:
    # one row per replicate, one column per rate (<type>_rate_reference and <type>_rate_final_global_weighted)
    frates=$bsdir2/bootstrap_rates-"$suffix".tsv
    for i in {1..10}; do
      fsum=$bsdir2/MaxL-"$suffix"-btree-$i.txt.summary.txt
      if [[ -f $fsum ]]; then
        awk -v rep=$i -F'\t' '$1 ~ /_rate_(reference|final_global_weighted)$/ {print rep "\t" $1 "\t" $2}' "$fsum"
      else
        echo "Warning: $fsum not found, replicate $i skipped" >&2
      fi
    done | awk -F'\t' '
      { if (!($2 in seen)) { seen[$2] = 1; cols[++nc] = $2 }
        if (!($1 in rseen)) { rseen[$1] = 1; reps[++nr] = $1 }
        val[$1, $2] = $3 }
      END { printf "replicate"; for (c = 1; c <= nc; c++) printf "\t%s", cols[c]; printf "\n"
            for (r = 1; r <= nr; r++) { printf "%s", reps[r]; for (c = 1; c <= nc; c++) printf "\t%s", val[reps[r], cols[c]]; printf "\n" } }' > "$frates"

    # Which rate the CIs are computed from: reference or final_global_weighted.
    # With bsr_mode=0 the two are the same. With bsr_mode>0, reference is the rate the branch
    # multipliers scale, final_global_weighted is the tree-wide average after the multipliers.
    ci_rate=reference
    # rate_column <column name> <output file>: write "<replicate> <rate>" lines, the input format of compute_ci.R
    rate_column () {
      awk -F'\t' -v col="$1" 'NR == 1 { for (c = 1; c <= NF; c++) if ($c == col) k = c; next } k { print $1 " " $k }' "$frates" > "$2"
    }
    fdup=$bsdir2/duprate-"$suffix".txt
    fdel=$bsdir2/delrate-"$suffix".txt
    rate_column dup_rate_$ci_rate "$fdup"
    rate_column del_rate_$ci_rate "$fdel"

    # Draw the ML tree with confidence intervals
    if [[ $plot -eq 1 ]]; then
      # N: CI from the 2.5th and 97.5th percentiles (Y: assume normal distribution)
      ci_dup=`Rscript "$REPO_ROOT/util/compute_ci.R" $fdup 4 N`
      ci_del=`Rscript "$REPO_ROOT/util/compute_ci.R" $fdel 4 N`
      echo $ci_dup
      echo $ci_del
      # Node-age CIs come from the replicates above (time trees only, not the .nmut.nex files).
      # If the tree was built in mode 0 with bootstrapping, its branch support is added as well.
      nex_pattern='btree-[0-9]+\.txt\.nex$'
      tree_name=$(basename "$mltree" .txt)
      if [[ $fixt_tree -eq 0 ]] && ls $bsdir1/"$tree_name"-btree-*.txt >/dev/null 2>&1; then
        Rscript "$REPO_ROOT/util/plot-trees-all.R" -s $bsdir1 -f $mltree2 --tree_file_nex "$mltree".nex -o $dir/MaxL-tree-"$suffix"-fixt-bootstrap.pdf -t "bootstrap" -l "ci" --time_file "$times" -p "$tree_name-btree-*txt" --bstrap_dir2 $bsdir2 --nex_pattern "$nex_pattern" --title "duplication rate $ci_dup; deletion rate $ci_del"
      else
        Rscript "$REPO_ROOT/util/plot-trees-all.R" -f $mltree2 --tree_file_nex "$mltree2".nex -o $dir/MaxL-tree-"$suffix"-fixt-ci.pdf -t "single" -l "ci" --time_file "$times" --bstrap_dir2 $bsdir2 --nex_pattern "$nex_pattern" --title "duplication rate $ci_dup; deletion rate $ci_del"
      fi
    fi
  fi

elif [[ $mode -eq 4 ]]; then  # Inferring ancestral states of a given tree from copy number profile
  # In this mode, the mutation rates have to be specified
  suffix=m$model-o"$opt"-s"$tree_search"-cons${cons}-estmu${estmu}-"$prefix"
  mltree=$dir/MaxL-"$suffix".txt
  echo "Inferring ancestral states of the given tree $mltree"

  log=$dir/std_cnetml_state_"$suffix"
  echo "seed $seed" > $log
  record_seed "$dir" "$seed" "$(basename "$log")"
  "$CNETA_BIN/cnetml" -c $input -t "$times" --tree_file "$mltree" -o "$mltree" --is_bin $is_bin --incl_all $incl_all -s $Ns --is_total $is_total --m_max $m_max -d $model --cn_max $cn_max --cn_type $cn_type --constrained $cons --correct_bias $correct_bias -x $mu --dup_rate $r1 --del_rate $r2 --chr_gain_rate $r3 --chr_loss_rate $r4 --wgd_rate $r5 --verbose $verbose --mode $mode  $extra_opts --seed $seed >> $log
fi
