#!/bin/bash
#
# Unit tests for erepopt sparse genomes + pair-scoped SKF suffixes.
# No DFTB/skgen and no external chemistry datasets.
#
set -euo pipefail
cd "$(dirname "$0")"

ROOT="$(cd ../.. && pwd)"
CXX="${CXX:-g++}"
INCLUDES="-I${ROOT}/eigen3.37 -I${ROOT}/galib247/include -I${ROOT}/erepopt"
LDFLAGS="-L${ROOT}/galib247/lib"
LIBS="-lga"

# Ensure genes.o exists (built with the erepopt makefile / same flags is fine).
if [ ! -f "${ROOT}/erepopt/genes.o" ]; then
  echo "Building erepopt/genes.o ..."
  make -C "${ROOT}/erepopt" genes.o
fi

echo "Building test_genes ..."
${CXX} -std=c++11 -O2 ${INCLUDES} -o test_genes test_genes.cpp "${ROOT}/erepopt/genes.o" ${LDFLAGS} ${LIBS}

echo "Running test_genes ..."
./test_genes
echo "PASS"
