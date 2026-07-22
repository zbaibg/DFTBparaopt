#include <iomanip>
#include "genes.hpp"

#include <sstream>

using namespace std;

int dense_genome_length(const Erepobj& erep, const sddh& ddh) {
  int length = 0;
  for (size_t i = 0; i < erep.velem.size(); i++)
    length += erep.velem[i].lmax + 2;
  if (ddh.td3) length += ddh.d3.size();
  if (ddh.tdamph) length += 1;
  if (ddh.thubbardderivs) length += ddh.hubbardderivs.size();
  if (ddh.tdamphver2) length += ddh.damphver2.size();
  if (ddh.thubbards) length += ddh.hubbards.size();
  if (ddh.tvorbes) length += ddh.vorbes.size();
  return length;
}

int genome_length(const Erepobj& erep, const sddh& ddh) {
  int length = 0;
  for (size_t i = 0; i < erep.velem.size(); i++) {
    int nvars = erep.velem[i].lmax + 2;
    for (int j = 0; j < nvars; j++) {
      if (gene_is_free(erep.velem[i].radius[j].minr, erep.velem[i].radius[j].maxr))
        length++;
    }
  }
  if (ddh.td3) {
    for (size_t i = 0; i < ddh.d3.size(); i++)
      if (gene_is_free(ddh.d3[i].min, ddh.d3[i].max)) length++;
  }
  if (ddh.tdamph) {
    if (gene_is_free(ddh.damph.min, ddh.damph.max)) length++;
  }
  if (ddh.thubbardderivs) {
    for (size_t i = 0; i < ddh.hubbardderivs.size(); i++)
      if (gene_is_free(ddh.hubbardderivs[i].min, ddh.hubbardderivs[i].max)) length++;
  }
  if (ddh.tdamphver2) {
    for (size_t i = 0; i < ddh.damphver2.size(); i++)
      if (gene_is_free(ddh.damphver2[i].min, ddh.damphver2[i].max)) length++;
  }
  if (ddh.thubbards) {
    for (size_t i = 0; i < ddh.hubbards.size(); i++)
      if (gene_is_free(ddh.hubbards[i].min, ddh.hubbards[i].max)) length++;
  }
  if (ddh.tvorbes) {
    for (size_t i = 0; i < ddh.vorbes.size(); i++)
      if (gene_is_free(ddh.vorbes[i].min, ddh.vorbes[i].max)) length++;
  }
  return length;
}

void clamp_free_genome(GAGenome& g, Erepobj& erep, sddh& ddh) {
  GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
  int idx = 0;
  for (size_t i = 0; i < erep.velem.size(); i++) {
    int nvars = erep.velem[i].lmax + 2;
    for (int j = 0; j < nvars; j++) {
      if (gene_is_free(erep.velem[i].radius[j].minr, erep.velem[i].radius[j].maxr)) {
        double value = genome.gene(idx);
        value = GAMax(erep.velem[i].radius[j].minr, value);
        value = GAMin(erep.velem[i].radius[j].maxr, value);
        genome.gene(idx++, value);
      }
    }
  }
  if (ddh.td3) {
    for (size_t i = 0; i < ddh.d3.size(); i++) {
      if (gene_is_free(ddh.d3[i].min, ddh.d3[i].max)) {
        double value = genome.gene(idx);
        value = GAMax(ddh.d3[i].min, value);
        value = GAMin(ddh.d3[i].max, value);
        genome.gene(idx++, value);
      }
    }
  }
  if (ddh.tdamph) {
    if (gene_is_free(ddh.damph.min, ddh.damph.max)) {
      double value = genome.gene(idx);
      value = GAMax(ddh.damph.min, value);
      value = GAMin(ddh.damph.max, value);
      genome.gene(idx++, value);
    }
  }
  if (ddh.thubbardderivs) {
    for (size_t i = 0; i < ddh.hubbardderivs.size(); i++) {
      if (gene_is_free(ddh.hubbardderivs[i].min, ddh.hubbardderivs[i].max)) {
        double value = genome.gene(idx);
        value = GAMax(ddh.hubbardderivs[i].min, value);
        value = GAMin(ddh.hubbardderivs[i].max, value);
        genome.gene(idx++, value);
      }
    }
  }
  if (ddh.tdamphver2) {
    for (size_t i = 0; i < ddh.damphver2.size(); i++) {
      if (gene_is_free(ddh.damphver2[i].min, ddh.damphver2[i].max)) {
        double value = genome.gene(idx);
        value = GAMax(ddh.damphver2[i].min, value);
        value = GAMin(ddh.damphver2[i].max, value);
        genome.gene(idx++, value);
      }
    }
  }
  if (ddh.thubbards) {
    for (size_t i = 0; i < ddh.hubbards.size(); i++) {
      if (gene_is_free(ddh.hubbards[i].min, ddh.hubbards[i].max)) {
        double value = genome.gene(idx);
        value = GAMax(ddh.hubbards[i].min, value);
        value = GAMin(ddh.hubbards[i].max, value);
        genome.gene(idx++, value);
      }
    }
  }
  if (ddh.tvorbes) {
    for (size_t i = 0; i < ddh.vorbes.size(); i++) {
      if (gene_is_free(ddh.vorbes[i].min, ddh.vorbes[i].max)) {
        double value = genome.gene(idx);
        value = GAMax(ddh.vorbes[i].min, value);
        value = GAMin(ddh.vorbes[i].max, value);
        genome.gene(idx++, value);
      }
    }
  }
}

void write_free_genes(ostream& out, const GAGenome& g,
                      const Erepobj& erep, const sddh& ddh) {
  const GA1DArrayGenome<double>& genome =
      (const GA1DArrayGenome<double>&)g;
  int idx = 0;
  out << std::fixed;
  for (size_t j = 0; j < erep.velem.size(); j++) {
    for (int k = 0; k < erep.velem[j].lmax + 2; k++) {
      if (gene_is_free(erep.velem[j].radius[k].minr, erep.velem[j].radius[k].maxr)) {
        out.precision(erep.velem[j].radius[k].precision);
        out.width(5);
        out << genome.gene(idx++) << "\t";
      }
    }
  }
  if (ddh.td3) {
    for (size_t j = 0; j < ddh.d3.size(); j++) {
      if (gene_is_free(ddh.d3[j].min, ddh.d3[j].max)) {
        out.precision(ddh.d3[j].precision);
        out.width(5);
        out << genome.gene(idx++) << "\t";
      }
    }
  }
  if (ddh.tdamph) {
    if (gene_is_free(ddh.damph.min, ddh.damph.max)) {
      out.precision(ddh.damph.precision);
      out.width(5);
      out << genome.gene(idx++) << "\t";
    }
  }
  if (ddh.thubbardderivs) {
    for (size_t j = 0; j < ddh.hubbardderivs.size(); j++) {
      if (gene_is_free(ddh.hubbardderivs[j].min, ddh.hubbardderivs[j].max)) {
        out.precision(ddh.hubbardderivs[j].precision);
        out.width(5);
        out << genome.gene(idx++) << "\t";
      }
    }
  }
  if (ddh.tdamphver2) {
    for (size_t j = 0; j < ddh.damphver2.size(); j++) {
      if (gene_is_free(ddh.damphver2[j].min, ddh.damphver2[j].max)) {
        out.precision(ddh.damphver2[j].precision);
        out.width(5);
        out << genome.gene(idx++) << "\t";
      }
    }
  }
  if (ddh.thubbards) {
    for (size_t j = 0; j < ddh.hubbards.size(); j++) {
      if (gene_is_free(ddh.hubbards[j].min, ddh.hubbards[j].max)) {
        out.precision(ddh.hubbards[j].precision);
        out.width(5);
        out << genome.gene(idx++) << "\t";
      }
    }
  }
  if (ddh.tvorbes) {
    for (size_t j = 0; j < ddh.vorbes.size(); j++) {
      if (gene_is_free(ddh.vorbes[j].min, ddh.vorbes[j].max)) {
        out.precision(ddh.vorbes[j].precision);
        out.width(5);
        out << genome.gene(idx++) << "\t";
      }
    }
  }
}

void write_individual_line(ostream& out, const GAGenome& g,
                           const Erepobj& erep, const sddh& ddh) {
  write_free_genes(out, g, erep, ddh);
  out << std::scientific;
  out.precision(12);
  out.width(20);
  out << g.score() << "\n";
  out << std::fixed;
}

void apply_genome_to_params(const GAGenome& g, Erepobj& erep, sddh& ddh) {
  const GA1DArrayGenome<double>& genome =
      (const GA1DArrayGenome<double>&)g;
  int idx = 0;
  for (size_t i = 0; i < erep.velem.size(); i++) {
    int nvars = erep.velem[i].lmax + 2;
    for (int j = 0; j < nvars; j++) {
      if (gene_is_free(erep.velem[i].radius[j].minr, erep.velem[i].radius[j].maxr)) {
        double value = genome.gene(idx++);
        value = GAMax(erep.velem[i].radius[j].minr, value);
        value = GAMin(erep.velem[i].radius[j].maxr, value);
        erep.velem[i].radius[j].r = value;
      }
      // else keep radius[j].r as loaded from input (min==max==r)
    }
  }
  if (ddh.td3) {
    for (size_t i = 0; i < ddh.d3.size(); i++) {
      if (gene_is_free(ddh.d3[i].min, ddh.d3[i].max)) {
        double value = genome.gene(idx++);
        value = GAMax(ddh.d3[i].min, value);
        value = GAMin(ddh.d3[i].max, value);
        ddh.d3[i].value = value;
      }
    }
  }
  if (ddh.tdamph) {
    if (gene_is_free(ddh.damph.min, ddh.damph.max)) {
      double value = genome.gene(idx++);
      value = GAMax(ddh.damph.min, value);
      value = GAMin(ddh.damph.max, value);
      ddh.damph.value = value;
    }
  }
  if (ddh.thubbardderivs) {
    for (size_t i = 0; i < ddh.hubbardderivs.size(); i++) {
      if (gene_is_free(ddh.hubbardderivs[i].min, ddh.hubbardderivs[i].max)) {
        double value = genome.gene(idx++);
        value = GAMax(ddh.hubbardderivs[i].min, value);
        value = GAMin(ddh.hubbardderivs[i].max, value);
        ddh.hubbardderivs[i].value = value;
      }
    }
  }
  if (ddh.tdamphver2) {
    for (size_t i = 0; i < ddh.damphver2.size(); i++) {
      if (gene_is_free(ddh.damphver2[i].min, ddh.damphver2[i].max)) {
        double value = genome.gene(idx++);
        value = GAMax(ddh.damphver2[i].min, value);
        value = GAMin(ddh.damphver2[i].max, value);
        ddh.damphver2[i].value = value;
      }
    }
  }
  if (ddh.thubbards) {
    for (size_t i = 0; i < ddh.hubbards.size(); i++) {
      if (gene_is_free(ddh.hubbards[i].min, ddh.hubbards[i].max)) {
        double value = genome.gene(idx++);
        value = GAMax(ddh.hubbards[i].min, value);
        value = GAMin(ddh.hubbards[i].max, value);
        ddh.hubbards[i].value = value;
      }
    }
  }
  if (ddh.tvorbes) {
    for (size_t i = 0; i < ddh.vorbes.size(); i++) {
      if (gene_is_free(ddh.vorbes[i].min, ddh.vorbes[i].max)) {
        double value = genome.gene(idx++);
        value = GAMax(ddh.vorbes[i].min, value);
        value = GAMin(ddh.vorbes[i].max, value);
        ddh.vorbes[i].value = value;
      }
    }
  }
}

string pair_onsite_suffix(const sddh& ddh, const string& e1, const string& e2) {
  ostringstream ss;
  string part;
  if (ddh.thubbards) {
    bool any = false;
    for (size_t i = 0; i < ddh.hubbards.size(); i++) {
      if (ddh.hubbards[i].name != e1 && ddh.hubbards[i].name != e2) continue;
      if (!any) {
        part += "_hub";
        any = true;
      }
      ss.str(string());
      ss.precision(ddh.hubbards[i].precision);
      ss << std::fixed << ddh.hubbards[i].value;
      part += "_" + ddh.hubbards[i].name + ss.str();
    }
  }
  if (ddh.tvorbes) {
    bool any = false;
    for (size_t i = 0; i < ddh.vorbes.size(); i++) {
      if (ddh.vorbes[i].name != e1 && ddh.vorbes[i].name != e2) continue;
      if (!any) {
        part += "_vor";
        any = true;
      }
      ss.str(string());
      ss.precision(ddh.vorbes[i].precision);
      ss << std::fixed << ddh.vorbes[i].value;
      part += "_" + ddh.vorbes[i].name + ss.str();
    }
  }
  return part;
}
