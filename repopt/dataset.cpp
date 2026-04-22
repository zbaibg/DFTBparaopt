#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <cstdlib>
#include <string>
#include <vector>
#include <map>

#include "allequations.hpp"
#include "dataset.hpp"

using namespace std;

string dataset_load_path = "";
string dataset_save_path = "";
map<string, DftbMolData> dftb_dataset_cache;
bool   dftb_dataset_loaded = false;

// The dataset is written with `setprecision(max_digits10)` (17 for IEEE
// double) in scientific notation. This is the shortest decimal form that
// round-trips every finite double back to the exact same bit pattern, so
// reloading the cache introduces no additional fitting error relative to
// computing eel/fel directly from DFTB.
static const int DATASET_PRECISION = numeric_limits<double>::max_digits10;

void load_dftb_dataset(const string& path) {
  if (path.empty()) { dftb_dataset_loaded = true; return; }

  ifstream fin(path.c_str());
  if (!fin) {
    cerr << endl << "ERROR: cannot open dataset_load file \"" << path << "\"" << endl
         << "exit repopt" << endl << endl;
    exit(1);
  }

  dftb_dataset_cache.clear();

  string line;
  while (getline(fin, line)) {
    // strip comments and blanks
    size_t hpos = line.find('#');
    if (hpos != string::npos) line = line.substr(0, hpos);
    {
      size_t first = line.find_first_not_of(" \t\r\n");
      if (first == string::npos) continue;
    }
    istringstream iss(line);
    string tok;
    iss >> tok;
    if (tok != "&molecule") continue;

    string molname;
    iss >> molname;
    if (molname.empty()) {
      cerr << endl << "ERROR: malformed dataset file \"" << path
           << "\": &molecule entry without name" << endl
           << "exit repopt" << endl << endl;
      exit(1);
    }

    DftbMolData entry;
    entry.eel = 0.0;
    int natom = -1;
    bool have_eel = false, have_natom = false, have_fel = false;
    bool end_seen = false;

    while (getline(fin, line)) {
      size_t hp = line.find('#');
      if (hp != string::npos) line = line.substr(0, hp);
      {
        size_t first = line.find_first_not_of(" \t\r\n");
        if (first == string::npos) continue;
      }
      istringstream iss2(line);
      string key;
      iss2 >> key;
      if (key == "&end") { end_seen = true; break; }
      else if (key == "natom") {
        iss2 >> natom;
        if (natom < 0) {
          cerr << endl << "ERROR: malformed natom for \"" << molname
               << "\" in dataset file \"" << path << "\"" << endl
               << "exit repopt" << endl << endl;
          exit(1);
        }
        entry.fel.resize(natom, 3);
        for (int i = 0; i < natom; i++)
          for (int j = 0; j < 3; j++)
            entry.fel(i, j) = 0.0;
        have_natom = true;
      } else if (key == "eel") {
        iss2 >> entry.eel;
        have_eel = true;
      } else if (key == "fel") {
        if (!have_natom) {
          cerr << endl << "ERROR: dataset file \"" << path
               << "\" has `fel` before `natom` for molecule \"" << molname << "\"" << endl
               << "exit repopt" << endl << endl;
          exit(1);
        }
        for (int i = 0; i < natom; i++) {
          if (!getline(fin, line)) {
            cerr << endl << "ERROR: unexpected EOF while reading fel for \""
                 << molname << "\" in dataset file \"" << path << "\"" << endl
                 << "exit repopt" << endl << endl;
            exit(1);
          }
          size_t hp2 = line.find('#');
          if (hp2 != string::npos) line = line.substr(0, hp2);
          istringstream fiss(line);
          if (!(fiss >> entry.fel(i, 0) >> entry.fel(i, 1) >> entry.fel(i, 2))) {
            cerr << endl << "ERROR: malformed fel row " << i << " for \""
                 << molname << "\" in dataset file \"" << path << "\"" << endl
                 << "exit repopt" << endl << endl;
            exit(1);
          }
        }
        have_fel = true;
      }
    }

    if (!end_seen || !have_eel || !have_natom || !have_fel) {
      cerr << endl << "ERROR: incomplete entry for \"" << molname
           << "\" in dataset file \"" << path << "\"" << endl
           << "  (need natom, eel, fel, &end)" << endl
           << "exit repopt" << endl << endl;
      exit(1);
    }

    dftb_dataset_cache[molname] = entry;
  }

  fin.close();
  dftb_dataset_loaded = true;
  cout << "#dataset: loaded " << dftb_dataset_cache.size()
       << " DFTB entries from \"" << path << "\"" << endl;
}

static void write_one(ofstream& fout, const Molecule& m) {
  fout << "&molecule " << m.name << "\n";
  fout << "  natom " << m.natom << "\n";
  fout << "  eel   " << scientific << setprecision(DATASET_PRECISION)
       << m.eel << "\n";
  fout << "  fel\n";
  for (int i = 0; i < m.natom; i++) {
    fout << "    "
         << scientific << setprecision(DATASET_PRECISION) << setw(26) << m.fel(i, 0) << " "
         << scientific << setprecision(DATASET_PRECISION) << setw(26) << m.fel(i, 1) << " "
         << scientific << setprecision(DATASET_PRECISION) << setw(26) << m.fel(i, 2) << "\n";
  }
  fout << "&end\n\n";
}

void save_dftb_dataset(const string& path,
                       const vector<Molecule>& vmol,
                       const vector<Molecule>& vreamol) {
  if (path.empty()) return;

  ofstream fout(path.c_str());
  if (!fout) {
    cerr << endl << "ERROR: cannot open dataset_save file \"" << path
         << "\" for writing" << endl
         << "exit repopt" << endl << endl;
    exit(1);
  }

  fout << "# DFTBparaopt DFTB dataset (repopt)\n"
       << "# Per-molecule cache of the DFTB total Mermin free energy (eel)\n"
       << "# and force matrix (fel, natom x 3) in atomic units, as parsed\n"
       << "# from detailed.out. Numbers are written in scientific notation\n"
       << "# with " << DATASET_PRECISION
       << " significant digits, which round-trips IEEE doubles exactly.\n"
       << "#\n"
       << "# NOTE: the cache is keyed by the xyz filename used in the\n"
       << "# `$compounds:` / `$definition_reactions:` blocks. If you change\n"
       << "# the DFTB Hamiltonian settings (dftb_in.hsd template, d3/damph/\n"
       << "# hubbardderivs/damphver2, SKF directory, ...), regenerate the\n"
       << "# dataset instead of reusing this one.\n\n";

  // Deduplicate against entries that are trivially copies of vmol entries
  // (vreamol may share its xyz with vmol, in which case eel/fel are the
  // same object-by-value). Tracking by name keeps the written file
  // self-consistent.
  map<string, bool> seen;
  int nwritten = 0;
  for (size_t i = 0; i < vmol.size(); i++) {
    if (seen.count(vmol[i].name)) continue;
    seen[vmol[i].name] = true;
    write_one(fout, vmol[i]);
    nwritten++;
  }
  for (size_t i = 0; i < vreamol.size(); i++) {
    if (seen.count(vreamol[i].name)) continue;
    seen[vreamol[i].name] = true;
    write_one(fout, vreamol[i]);
    nwritten++;
  }

  fout.close();
  cout << "#dataset: wrote " << nwritten
       << " DFTB entries to \"" << path << "\"" << endl;
}
