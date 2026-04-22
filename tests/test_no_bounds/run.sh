#!/bin/bash
# Test: grid files contain only knot values (no UP/LOW bounds).
# Expected result: identical to examples/rep1.out_ref_new, since
# apply_knot_bounds is a no-op when no UP/LOW is declared. We ignore the
# timing line which is non-deterministic.

set -e
cd "$(dirname "$0")"

REPOPT=../../repopt/repopt
REF=rep1.out_ref_new

srun --partition=pre "$REPOPT" rep1.in >& rep1.out

# Ignore the "total run time" line (varies per run) and the srun queue /
# allocation messages that are injected before repopt's own stdout.
FILTER='/total run time/d;/^srun:/d'
diff <(sed "$FILTER" rep1.out) <(sed "$FILTER" "$REF") > rep1.diff || true
if [ -s rep1.diff ]; then
  echo "FAIL: rep1.out differs from $REF (see rep1.diff)"
  head -40 rep1.diff
  exit 1
else
  echo "PASS: rep1.out matches $REF (ignoring run-time line)"
fi
