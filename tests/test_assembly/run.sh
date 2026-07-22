#!/bin/bash
#
# Test: pair-table assembly cache vs legacy assemblers (bit-exact).
#
# Step 1: run once with dataset_save so DFTB results are cached.
# Step 2: run with dataset_load + verify_assembly N. That mode builds
#         eqmat/vref with the cached pair-table path and with the original
#         findcolind path for N perturbed knot vectors and requires every
#         entry of eqmat, vref and vweight to match bit-for-bit.
#
# geom/grids/dftb_inp are symlinks into tests/test_no_bounds/.
#
set -euo pipefail
cd "$(dirname "$0")"

REPOPT=../../repopt/repopt
if [ ! -x "$REPOPT" ]; then
  echo "FAIL: repopt binary not found: $REPOPT"
  exit 1
fi
LAUNCHER=${REPOPT_LAUNCHER:-}

run_repopt() {
  local inp=$1 out=$2
  if [ -n "$LAUNCHER" ]; then
    $LAUNCHER "$REPOPT" "$inp" >"$out" 2>&1
  else
    "$REPOPT" "$inp" >"$out" 2>&1
  fi
}

rm -rf dftb.tmp dftb_cache.dat rep_save.out rep_verify.out

echo "[1/2] dataset_save (may call DFTB)..."
run_repopt rep_save.in rep_save.out
if [ ! -s dftb_cache.dat ]; then
  echo "FAIL: no dftb_cache.dat after save"
  tail -30 rep_save.out || true
  exit 1
fi

echo "[2/2] verify_assembly against legacy ..."
run_repopt rep_verify.in rep_verify.out

if ! grep -q 'verify_assembly: all trials matched (bit-exact)' rep_verify.out; then
  echo "FAIL: verify_assembly did not report bit-exact match"
  tail -50 rep_verify.out || true
  exit 1
fi
if grep -q 'verify_assembly FAIL' rep_verify.out; then
  echo "FAIL: mismatch reported"
  exit 1
fi

echo "PASS: assembly cache matches legacy eqmat/vref/vweight"
