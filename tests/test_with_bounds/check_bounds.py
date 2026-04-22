#!/usr/bin/env python3
"""
Verify that the c_c_ inner knots ended up inside the UP/LOW bounds that
were declared in grids/cc.grdx. Knot 2 had UP == LOW, so it should be
frozen at 1.905 A. Output spline intervals are printed in Bohr, so we
convert back to Angstrom for readability.
"""
import re
import sys

AA_BOHR = 0.5291772085936
TOL = 1e-3  # A

# expected bounds in Angstrom: (lo, hi) for inner knots 1, 2, 3.
# -inf / +inf is used where the corresponding bound is intentionally
# omitted in grids/cc.grdx, to exercise the "UP only" and "LOW only"
# code paths.
expected = {
    1: (float("-inf"), 1.80),    # only UP declared
    2: (1.905, 1.905),            # UP == LOW  (frozen)
    3: (2.20, float("inf")),      # only LOW declared
}

def parse_cc_knots(path):
    with open(path) as f:
        lines = f.readlines()
    # Header layout produced by writeout() for each potential:
    #   c_c_   r1,r2,Allequations ( 4 th order splines)
    #   Spline
    #   <nknots-1>  <cutoff>
    #   <exponential expA expB expC>
    #   <knot_i  knot_{i+1}  coeff0..coeff4>    <-- one row per spline piece
    for i, line in enumerate(lines):
        if line.lstrip().startswith("c_c_") and "r1,r2" in line:
            start = i + 4
            break
    else:
        sys.exit("Could not locate c_c_ spline block in " + path)
    knots_bohr = []
    for line in lines[start:]:
        parts = line.split()
        if len(parts) < 2:
            break
        try:
            r1 = float(parts[0]); r2 = float(parts[1])
        except ValueError:
            break
        if not knots_bohr:
            knots_bohr.append(r1)
        knots_bohr.append(r2)
    return knots_bohr

def main():
    path = sys.argv[1] if len(sys.argv) > 1 else "rep1.out"
    knots_bohr = parse_cc_knots(path)
    knots_aa = [k * AA_BOHR for k in knots_bohr]
    print("c_c_ knots (A):", ["%.4f" % k for k in knots_aa])
    ok = True
    for idx, (lo, hi) in expected.items():
        v = knots_aa[idx]
        if v < lo - TOL or v > hi + TOL:
            print("FAIL: knot %d = %.4f A not in [%.4f, %.4f]" % (idx, v, lo, hi))
            ok = False
        else:
            print("OK:   knot %d = %.4f A in [%.4f, %.4f]" % (idx, v, lo, hi))
    if not ok:
        sys.exit(1)
    print("PASS: all c_c_ inner knots respect the UP/LOW bounds")

if __name__ == "__main__":
    main()
