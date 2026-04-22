#!/bin/bash
# Test: grids/cc.grdx declares UP/LOW bounds, including one knot with
# UP == LOW (forced to 1.905 A). After the GA finishes the c_c_ inner
# knots must be inside the declared bounds.

set -e
cd "$(dirname "$0")"

REPOPT=../../repopt/repopt

srun --partition=pre "$REPOPT" rep1.in >& rep1.out

python3 check_bounds.py rep1.out
