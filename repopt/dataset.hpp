#ifndef DATASET_HPP_INCLUDED
#define DATASET_HPP_INCLUDED

#include <string>
#include <map>
#include <vector>
#include <Eigen/Dense>

// DFTB-level per-molecule data that is expensive to (re)compute but does
// not depend on the spline parameters optimized by the GA. Caching these
// quantities lets a user rerun the GA with different hyper-parameters
// without calling the DFTB driver again. Only quantities that come out of
// the DFTB run are stored: `eel` (Total Mermin free energy) and the
// `fel(natom,3)` force matrix. Geometry, reference energies/forces are
// always re-read from the input files.
struct DftbMolData {
  double          eel;
  Eigen::MatrixXd fel;   // natom x 3, in atomic units (consistent with detailed.out)
};

// Paths controlled by the `$system:` block:
//   dataset_load <path>   - if non-empty, skip DFTB for any molecule whose
//                           xyz filename appears in <path> and reuse the
//                           cached eel/fel from there. Missing molecules
//                           fall back to a real DFTB call, so incremental
//                           caches work.
//   dataset_save <path>   - if non-empty, write the current (post-init)
//                           DFTB dataset to <path> with full double
//                           precision so it can be reused by a later run.
extern std::string dataset_load_path;
extern std::string dataset_save_path;

// In-memory dataset. Keyed by the xyz filename (Molecule::name) exactly as
// given in the input file, so the same string is used for lookup.
extern std::map<std::string, DftbMolData> dftb_dataset_cache;

// One-shot guard so that Molecule::init loads the cache lazily on its
// first invocation, regardless of which section ($compounds: or
// $definition_reactions:) triggers it.
extern bool dftb_dataset_loaded;

class Molecule;

void load_dftb_dataset(const std::string& path);
void save_dftb_dataset(const std::string& path,
                       const std::vector<Molecule>& vmol,
                       const std::vector<Molecule>& vreamol);

#endif /* DATASET_HPP_INCLUDED */
