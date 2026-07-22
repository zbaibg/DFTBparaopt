#ifndef GENES_HPP
#define GENES_HPP

#include <cmath>
#include <iostream>
#include <ga/ga.h>
#include "erepobj.hpp"
#include "ddh.hpp"

// Helpers for sparse genomes: parameters with min==max are excluded from the GA
// chromosome and kept at their input/default values on erepobj / ddh.

inline bool gene_is_free(double amin, double amax) {
  return std::fabs(amax - amin) > 1e-14;
}

int genome_length(const Erepobj& erep, const sddh& ddh);

// Dense (legacy) chromosome length including fixed parameters.
int dense_genome_length(const Erepobj& erep, const sddh& ddh);

// Clamp only free gene slots in the sparse chromosome.
void clamp_free_genome(GAGenome& g, Erepobj& erep, sddh& ddh);

// Write free genes with per-parameter precision (pop.*.dat format).
void write_free_genes(std::ostream& out, const GAGenome& g,
                      const Erepobj& erep, const sddh& ddh);

// Write one individual: free genes + score (scientific, 12 digits).
void write_individual_line(std::ostream& out, const GAGenome& g,
                           const Erepobj& erep, const sddh& ddh);

// Copy free genes from genome into velem[].radius[].r and ddh.*.value.
// Fixed parameters are left unchanged (as loaded from input).
void apply_genome_to_params(const GAGenome& g, Erepobj& erep, sddh& ddh);

// Build SKF filename suffix from hubbard/vorbes that belong to either element
// of the pair (not the full global list).
std::string pair_onsite_suffix(const sddh& ddh,
                               const std::string& e1,
                               const std::string& e2);

#endif /* GENES_HPP */
