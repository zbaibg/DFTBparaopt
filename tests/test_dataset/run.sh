#!/bin/bash
#
# Test: DFTB dataset save/load round-trip.
#
# Step 1 ("save") runs repopt normally and lets it cache every molecule's
#   DFTB eel/fel into dftb_cache.dat (high-precision, round-trip-safe).
# Step 2 ("load") runs the same GA problem with identical seed but uses
#   `dataset_load` instead of calling DFTB. If the cache is faithful the
#   whole optimisation -- final score, spline coefficients, initial and
#   final populations -- must reproduce bit-for-bit.
#
# The test passes iff
#   (a) step 2 never invokes DFTB (the dftb.tmp/ working dir must be
#       absent after the load run), and
#   (b) rep1_save.out and rep1_load.out are identical once the non-
#       deterministic timing line is stripped.
#
# This test lives next to test_no_bounds and does NOT touch any files
# from the other test directories. The geom/, grids/, dftb_inp/ folders
# are symlinks into tests/test_no_bounds/ so no data is duplicated.

set -e
cd "$(dirname "$0")"

REPOPT=../../repopt/repopt

if [ ! -x "$REPOPT" ]; then
  echo "FAIL: repopt binary not found or not executable: $REPOPT"
  echo "      (build it with: cd repopt && make)"
  exit 1
fi

# Optional launcher. Default to direct execution; export REPOPT_LAUNCHER
# (e.g. REPOPT_LAUNCHER="srun --partition=pre") to mirror the style of
# the other test scripts on a SLURM system.
LAUNCHER=${REPOPT_LAUNCHER:-}

run_repopt() {
  local inp=$1
  local out=$2
  if [ -n "$LAUNCHER" ]; then
    $LAUNCHER "$REPOPT" "$inp" >& "$out"
  else
    "$REPOPT" "$inp" >& "$out"
  fi
}

# --- Step 1: populate the cache ----------------------------------------
rm -rf dftb.tmp dftb_cache.dat rep1_save.out rep1_load.out rep1.diff \
       score.dat pop.initial.dat pop.final.dat

echo "[1/2] running repopt with dataset_save (DFTB will be executed)..."
run_repopt rep1_save.in rep1_save.out

if [ ! -s dftb_cache.dat ]; then
  echo "FAIL: dataset_save did not produce dftb_cache.dat"
  exit 1
fi

# --- Step 2: reuse the cache (DFTB must NOT be invoked) ----------------
# Remove dftb.tmp/ so we can detect whether the load run re-creates it.
# (Molecule::init unconditionally rm -rf's + mkdir's that directory at
# the start of a DFTB call, so its reappearance signals a cache miss.)
rm -rf dftb.tmp

echo "[2/2] running repopt with dataset_load (DFTB must be skipped)..."
run_repopt rep1_load.in rep1_load.out

if [ -d dftb.tmp ]; then
  echo "FAIL: dftb.tmp/ reappeared during the dataset_load run --"
  echo "      DFTB was invoked even though the cache should have covered"
  echo "      every molecule."
  exit 1
fi

# --- Compare the two runs ----------------------------------------------
# Ignore the non-deterministic timing line and the "#dataset:" bookkeeping
# lines (one says "wrote N", the other "loaded N"). Everything else --
# final score, spline coefficients, population dumps -- must match.
FILTER='/total run time/d;/^#dataset:/d;/^srun:/d'
diff <(sed "$FILTER" rep1_save.out) <(sed "$FILTER" rep1_load.out) > rep1.diff || true

if [ -s rep1.diff ]; then
  echo "FAIL: rep1_save.out and rep1_load.out differ (see rep1.diff)"
  head -40 rep1.diff
  exit 1
fi

echo "PASS: dataset save/load round-trip reproduces the save run exactly"
echo "      (wrote $(grep -c '^&molecule ' dftb_cache.dat) molecules to dftb_cache.dat)"
