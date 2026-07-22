// Unit tests for sparse genomes and pair-scoped SKF onsite suffixes.
// Does not run DFTB/skgen and does not use any external chemistry datasets.
//
// Build/run: ./run.sh   (from this directory)

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <ga/ga.h>
#include "../../erepopt/genes.hpp"
#include "../../erepopt/erepobj.hpp"
#include "../../erepopt/ddh.hpp"
#include "../../erepopt/element.hpp"

using namespace std;

// Stub only the default ctor so we can build an in-memory Erepobj without
// linking the full erepobj.o dependency graph (DFTB I/O, molecules, etc.).
Erepobj::Erepobj() {}

static int g_fails = 0;

static void expect(bool ok, const string& msg) {
  if (!ok) {
    cerr << "FAIL: " << msg << endl;
    g_fails++;
  } else {
    cout << "OK: " << msg << endl;
  }
}

static Radius make_radius(double r, double amin, double amax) {
  Radius x;
  x.r = r;
  x.minr = amin;
  x.maxr = amax;
  x.delta = 0.0;
  x.precision = 2;
  return x;
}

static sdatapoint make_dp(const string& name, const string& type,
                          double amin, double val, double amax, int prec = 2) {
  sdatapoint p;
  p.name = name;
  p.type = type;
  p.precision = prec;
  p.min = amin;
  p.value = val;
  p.max = amax;
  p.delta = 0.0;
  return p;
}

static Element make_elem(const string& name, int lmax, const vector<Radius>& rad) {
  Element e;
  e.name = name;
  e.optype = 0;
  e.lmax = lmax;
  e.radius = rad;
  return e;
}

int main() {
  // -------------------------------------------------------------------------
  // 1) Sparse vs dense genome length with locked radii / onsite
  // -------------------------------------------------------------------------
  Erepobj erep;
  // dens + s + p  => lmax=1 => nvars = lmax+2 = 3
  // dens free, s locked (min==max), p free
  erep.velem.push_back(make_elem("C", 1, {
      make_radius(3.0, 2.5, 3.5),  // dens free
      make_radius(2.0, 2.0, 2.0),  // s locked
      make_radius(2.5, 2.0, 3.0),  // p free
  }));
  // dens free, s free  => lmax=0 => nvars=2
  erep.velem.push_back(make_elem("H", 0, {
      make_radius(1.5, 1.0, 2.0),
      make_radius(1.2, 1.0, 1.5),
  }));

  sddh ddh;
  ddh.td3 = ddh.tdamph = ddh.thubbardderivs = ddh.tdamphver2 = false;
  ddh.thubbards = true;
  ddh.tvorbes = true;
  ddh.thirdorderfull = false;
  // C Hubbard free, H Hubbard locked, O Hubbard free (third element for suffix test)
  ddh.hubbards.push_back(make_dp("C", "s", 0.3, 0.40, 0.5));
  ddh.hubbards.push_back(make_dp("H", "s", 0.4, 0.40, 0.4)); // locked
  ddh.hubbards.push_back(make_dp("O", "s", 0.2, 0.25, 0.3));
  ddh.vorbes.push_back(make_dp("C", "2S", -0.3, -0.25, -0.2));
  ddh.vorbes.push_back(make_dp("H", "1S", -0.2, -0.2, -0.2)); // locked
  ddh.vorbes.push_back(make_dp("O", "2S", -0.4, -0.35, -0.3));

  const int dense = dense_genome_length(erep, ddh);
  const int sparse = genome_length(erep, ddh);
  // dense: C(3)+H(2)+hub(3)+vor(3) = 11
  expect(dense == 11, "dense_genome_length == 11");
  // sparse free: C dens, C p, H dens, H s, C hub, O hub, C vor, O vor = 8
  // locked excluded: C s, H hub, H vor
  expect(sparse == 8, "genome_length excludes min==max (8 free genes)");
  expect(sparse < dense, "sparse length strictly less than dense");

  // -------------------------------------------------------------------------
  // 2) apply_genome_to_params leaves locked values alone
  // -------------------------------------------------------------------------
  GA1DArrayGenome<double> genome(sparse);
  // Fill free genes with distinctive values (within bounds where possible).
  // Order: C dens, C p, H dens, H s, C hub, O hub, C vor, O vor
  genome.gene(0, 3.1);
  genome.gene(1, 2.7);
  genome.gene(2, 1.6);
  genome.gene(3, 1.3);
  genome.gene(4, 0.45);
  genome.gene(5, 0.28);
  genome.gene(6, -0.22);
  genome.gene(7, -0.32);

  // Poison locked slots with wrong numbers; apply must not overwrite them.
  erep.velem[0].radius[1].r = 2.0;      // locked C s
  ddh.hubbards[1].value = 0.40;         // locked H hub
  ddh.vorbes[1].value = -0.20;          // locked H vor

  apply_genome_to_params(genome, erep, ddh);

  expect(std::fabs(erep.velem[0].radius[0].r - 3.1) < 1e-14, "free C dens applied");
  expect(std::fabs(erep.velem[0].radius[1].r - 2.0) < 1e-14, "locked C s unchanged");
  expect(std::fabs(erep.velem[0].radius[2].r - 2.7) < 1e-14, "free C p applied");
  expect(std::fabs(erep.velem[1].radius[0].r - 1.6) < 1e-14, "free H dens applied");
  expect(std::fabs(ddh.hubbards[0].value - 0.45) < 1e-14, "free C hub applied");
  expect(std::fabs(ddh.hubbards[1].value - 0.40) < 1e-14, "locked H hub unchanged");
  expect(std::fabs(ddh.hubbards[2].value - 0.28) < 1e-14, "free O hub applied");
  expect(std::fabs(ddh.vorbes[1].value - (-0.20)) < 1e-14, "locked H vor unchanged");

  // -------------------------------------------------------------------------
  // 3) pair_onsite_suffix only includes the pair's elements
  // -------------------------------------------------------------------------
  const string ch = pair_onsite_suffix(ddh, "C", "H");
  const string co = pair_onsite_suffix(ddh, "C", "O");
  const string hh = pair_onsite_suffix(ddh, "H", "H");

  expect(ch.find("_hub") != string::npos, "C-H suffix has _hub");
  expect(ch.find("_vor") != string::npos, "C-H suffix has _vor");
  expect(ch.find("_C") != string::npos, "C-H suffix includes C");
  expect(ch.find("_H") != string::npos, "C-H suffix includes H");
  expect(ch.find("_O") == string::npos, "C-H suffix excludes O");

  expect(co.find("_O") != string::npos, "C-O suffix includes O");
  expect(co.find("_H") == string::npos, "C-O suffix excludes H");

  expect(hh.find("_C") == string::npos, "H-H suffix excludes C");
  expect(hh.find("_O") == string::npos, "H-H suffix excludes O");
  expect(hh.find("_H") != string::npos, "H-H suffix includes H");

  // Changing an unrelated onsite (O) must not change the C-H SKF suffix.
  const string ch_before = ch;
  ddh.hubbards[2].value = 0.299;
  ddh.vorbes[2].value = -0.301;
  const string ch_after = pair_onsite_suffix(ddh, "C", "H");
  expect(ch_before == ch_after,
         "C-H SKF suffix unchanged when only O onsite values move");

  // Changing C onsite must change the C-H suffix.
  ddh.hubbards[0].value = 0.499;
  const string ch_changed = pair_onsite_suffix(ddh, "C", "H");
  expect(ch_changed != ch_before,
         "C-H SKF suffix changes when C hubbard moves");

  // -------------------------------------------------------------------------
  // 4) write_free_genes emits exactly sparse length tokens
  // -------------------------------------------------------------------------
  ostringstream oss;
  write_free_genes(oss, genome, erep, ddh);
  istringstream iss(oss.str());
  int ntok = 0;
  string tok;
  while (iss >> tok) ntok++;
  expect(ntok == sparse, "write_free_genes token count == genome_length");

  if (g_fails == 0) {
    cout << "PASS: all erepopt genes tests" << endl;
    return 0;
  }
  cerr << "FAILED: " << g_fails << " check(s)" << endl;
  return 1;
}
