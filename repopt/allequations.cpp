#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <chrono>
#include <algorithm>
#include <Eigen/Dense>
#include "allequations.hpp"
#include "auxiliary.hpp"
#include "tools.hpp"

using namespace std;

extern double min_step01;
using namespace Eigen;

extern bool runtest;
extern int ilmsfit,idecompose,nreplicate,dftbout_new;

Allequations::Allequations()
  : assembly_cache_ready(false), use_assembly_cache(true) {}

Allequations::Allequations(const string inputfile)
  : assembly_cache_ready(false), use_assembly_cache(true) {
  calculate(inputfile);    
}

void Allequations::calculate(const string inputfile) {
  int ieq=0;

  readinp(inputfile); 
  sort_inputdist();
  get_autogrids();
  sizeeqsys();    
  build_assembly_cache();
  ieq = include_splineeq(0);
  ieq = include_energyeq(ieq);
  ieq = include_forceeq(ieq);
  ieq = include_addeq(ieq);
  ieq = include_reaeq(ieq);
  ieq = include_smootheq(ieq);
  weighteqsys();
  svd_fulfill_spleq();
  get_residual();

}

void Allequations::get_minRbond(){
  string potname;
  double dist;
  int ipot,imol,irea,ind,at1,at2; 

  for (ipot=0; ipot < vpot.size() ;ipot++) {
    vpot[ipot].minRbond = 999.9;
  }

  for (imol=0; imol < vmol.size(); imol++){
    for (at1=0; at1 < vmol[imol].natom; at1++){
      for (at2=at1+1; at2 < vmol[imol].natom; at2++){
        dist = vmol[imol].dist(at1,at2); 
        potname = findpotname( vmol[imol].atomname[at1], vmol[imol].atomname[at2] ); 
        for (ipot=0; ipot < vpot.size() ;ipot++) {
          if ( vpot[ipot].potname == potname ) {
            if( vpot[ipot].minRbond>dist) vpot[ipot].minRbond = dist; 
            break;        
          }
        }
      }
    }
  }

  for (irea=0; irea < vrea.size(); irea++ ){
    for (imol=0; imol < vrea[irea].vmol.size(); imol++){
      ind=vrea[irea].vmol[imol];
      for (at1=0; at1 < vreamol[ind].natom; at1++){
        for (at2=at1+1; at2 < vreamol[ind].natom; at2++){
          dist = vreamol[ind].dist(at1,at2); 
          potname = findpotname( vreamol[ind].atomname[at1], vreamol[ind].atomname[at2] ); 
          for (ipot=0; ipot < vpot.size() ;ipot++) {
            if ( vpot[ipot].potname == potname ) {
              if( vpot[ipot].minRbond>dist) vpot[ipot].minRbond = dist; 
              break;        
            }
          }
        }
      }
    }
  }

  for (ipot=0; ipot < vpot.size() ;ipot++) {
    vpot[ipot].minRbond = vpot[ipot].minRbond - 0.0001 ; 
  }
}

void Allequations::validate_grid_knot_spacing(){
  // Feasibility check for the GA corridor of each potential's inner knots.
  //
  // Inner knots j = 1 .. nknots-2 must lie in
  //   [ max(minRbond+min_step01, has_low[j] ? vr_low[j] : -inf),
  //     min(cutoff-max_step,     has_up [j] ? vr_up [j] : +inf) ]
  // and consecutive inner knots must be separated by at least max_step
  // (enforced by the push-apart loops in MyObjective).
  //
  // We build per-knot effective lo/hi, then propagate the max_step spacing
  // in both directions so each knot's corridor also reflects the user-
  // supplied UP/LOW bounds of ALL other knots. If any knot ends up with
  // lo > hi, the grid is infeasible and we report which knot (and which
  // upstream/downstream UP or LOW caused the squeeze) before exiting.
  //
  // (min_step01 is read in Angstrom; convert to Bohr like main() does.)
  const double min_step01_bohr = min_step01 / AA_Bohr;
  const double eps   = 1e-9;
  const double NEGINF = -1e300;
  const double POSINF =  1e300;

  for (int ipot = 0; ipot < (int)vpot.size(); ipot++){
    if (vpot[ipot].gridname.size() >= 8 &&
        vpot[ipot].gridname.substr(0, 8) == "autogrid") continue;
    const int nk = vpot[ipot].nknots;
    if (nk < 2) continue;

    const double cutoff   = vpot[ipot].vr[nk - 1];
    const double max_step = vpot[ipot].max_step;
    const double minRbond = vpot[ipot].minRbond;

    // First check the base corridor (no UP/LOW): MyObjective always
    // clamps inner knots into [minRbond+min_step01, cutoff-max_step],
    // so this corridor must be non-empty even before any user bound.
    const double lo_base = minRbond + min_step01_bohr;
    const double hi_base = cutoff   - max_step;
    if (hi_base < lo_base - eps){
      cerr << endl << "ERROR:" << endl
           << "  Grid knot spacing (potential \"" << vpot[ipot].potname << "\", file \""
           << vpot[ipot].gridname << "\"):" << endl
           << "  (cutoff - max_step) < (minRbond + min_step01); no valid inner-knot interval." << endl
           << "  cutoff = " << cutoff * AA_Bohr << " A, max_step = " << max_step * AA_Bohr << " A  =>  cutoff-max_step = "
           << hi_base * AA_Bohr << " A," << endl
           << "  minRbond = " << minRbond * AA_Bohr << " A, min_step01 = " << min_step01 << " A  =>  minRbond+min_step01 = "
           << lo_base * AA_Bohr << " A." << endl
           << "  Increase cutoff, decrease max_step, or decrease minRbond / min_step01." << endl
           << "exit repopt" << endl << endl;
      exit(1);
    }

    // No inner knots -> only the base check is meaningful.
    if (nk <= 2) continue;

    // Build per-knot effective lo/hi, seeded with base + any user UP/LOW.
    // Index j ranges over inner knots 1..nk-2; we size the vectors to nk
    // for index convenience and leave endpoints unused.
    vector<double> eff_lo(nk, NEGINF);
    vector<double> eff_hi(nk, POSINF);
    // Track which "source" (base / LOW@j' / UP@j') is currently the
    // tightest on each side, so error messages pinpoint the cause.
    // src_lo[j] / src_hi[j] = -1 means the base (minRbond+min_step01 or
    // cutoff-max_step); >= 0 means the knot index whose user bound is
    // responsible. Initialized to -1 (base).
    vector<int>    src_lo(nk, -1);
    vector<int>    src_hi(nk, -1);

    for (int j = 1; j <= nk - 2; j++){
      eff_lo[j] = lo_base;
      eff_hi[j] = hi_base;
      if ((int)vpot[ipot].has_low.size() > j && vpot[ipot].has_low[j]){
        if (vpot[ipot].vr_low[j] > eff_lo[j]){
          eff_lo[j] = vpot[ipot].vr_low[j];
          src_lo[j] = j;
        }
      }
      if ((int)vpot[ipot].has_up.size() > j && vpot[ipot].has_up[j]){
        if (vpot[ipot].vr_up[j] < eff_hi[j]){
          eff_hi[j] = vpot[ipot].vr_up[j];
          src_hi[j] = j;
        }
      }
    }

    // Propagate the max_step spacing: knot j must sit at least
    // max_step above knot j-1's minimum, and at least max_step below
    // knot j+1's maximum. Iterated linearly over inner knots in each
    // direction (this is sufficient because max_step is constant).
    for (int j = 2; j <= nk - 2; j++){
      const double cand = eff_lo[j - 1] + max_step;
      if (cand > eff_lo[j]){
        eff_lo[j] = cand;
        src_lo[j] = src_lo[j - 1];  // propagate the original source
      }
    }
    for (int j = nk - 3; j >= 1; j--){
      const double cand = eff_hi[j + 1] - max_step;
      if (cand < eff_hi[j]){
        eff_hi[j] = cand;
        src_hi[j] = src_hi[j + 1];
      }
    }

    // Feasibility: every inner knot must have a non-empty corridor.
    for (int j = 1; j <= nk - 2; j++){
      if (eff_lo[j] > eff_hi[j] + eps){
        cerr << endl << "ERROR:" << endl
             << "  Grid knot spacing (potential \"" << vpot[ipot].potname << "\", file \""
             << vpot[ipot].gridname << "\"):" << endl
             << "  inner knot " << j << " has empty corridor under the current UP/LOW + max_step constraints." << endl
             << "  effective lower bound = " << eff_lo[j] * AA_Bohr << " A";
        if (src_lo[j] < 0){
          cerr << "  (from minRbond+min_step01 = " << lo_base * AA_Bohr
               << " A, possibly +i*max_step)";
        } else {
          cerr << "  (from LOW on knot " << src_lo[j] << " = "
               << vpot[ipot].vr_low[src_lo[j]] * AA_Bohr
               << " A, propagated via max_step)";
        }
        cerr << "," << endl
             << "  effective upper bound = " << eff_hi[j] * AA_Bohr << " A";
        if (src_hi[j] < 0){
          cerr << "  (from cutoff-max_step = " << hi_base * AA_Bohr
               << " A, possibly -i*max_step)";
        } else {
          cerr << "  (from UP on knot " << src_hi[j] << " = "
               << vpot[ipot].vr_up[src_hi[j]] * AA_Bohr
               << " A, propagated via max_step)";
        }
        cerr << "." << endl
             << "  cutoff = " << cutoff   * AA_Bohr << " A, max_step = " << max_step * AA_Bohr << " A," << endl
             << "  minRbond = " << minRbond * AA_Bohr << " A, min_step01 = " << min_step01 << " A, nknots = " << nk << "." << endl
             << "  Loosen the offending UP/LOW, reduce nknots, increase cutoff / minRbond margin," << endl
             << "  decrease max_step, or adjust min_step01." << endl
             << "exit repopt" << endl << endl;
        exit(1);
      }
    }
  }
}

void Allequations::prepare(const string inputfile) {
  int ieq=0;

  readinp(inputfile); 
  sort_inputdist();
  get_autogrids();
  sizeeqsys();    
  get_minRbond();
  validate_grid_knot_spacing();
  build_assembly_cache();
}

void Allequations::reset() {
  int i,j;
  for (i=0; i < nrows; i++){
    vref[i]    = 0.0;
    vweight[i] = 1.0;
    vres[i]    = 0.0;
    for (j=0; j < ncols; j++){
      eqmat(i,j) = 0.0;
    }
  }
  for (j=0; j < ncols; j++) {
    vunknown[j] = 0.0;
    sigma[j]    = 0.0;
  }
  for (i=0; i < nfree+natfit; i++){
    sigma_e_f[i] = 0.0;
  }
}

double Allequations::score() {

  Timer t;

  int i,ieq=0;

  //cout << "T0 taken: " << t.elapsed() << " seconds\n"; t.reset();
  ieq = include_splineeq(0);
  //cout << "T1 taken: " << t.elapsed() << " seconds\n"; t.reset();
  ieq = include_energyeq(ieq);
  //cout << "T2 taken: " << t.elapsed() << " seconds\n"; t.reset();
  ieq = include_forceeq(ieq);
  //cout << "T3 taken: " << t.elapsed() << " seconds\n"; t.reset();
  ieq = include_addeq(ieq);
  //cout << "T4 taken: " << t.elapsed() << " seconds\n"; t.reset();
  ieq = include_reaeq(ieq);
  //cout << "T5 taken: " << t.elapsed() << " seconds\n"; t.reset();
  ieq = include_smootheq(ieq);
  //cout << "T6 taken: " << t.elapsed() << " seconds\n"; t.reset();
  weighteqsys();
  //cout << "T7 taken: " << t.elapsed() << " seconds\n"; t.reset();
  if(!runtest){
    for (i=0; i < nreplicate; i++){
      svd_fulfill_spleq();
    }
    //cout << "T8 taken: " << t.elapsed() << " seconds\n"; t.reset();
    get_residual();
    //cout << "T9 taken: " << t.elapsed() << " seconds\n"; t.reset();
  }  
  return 1.0;
}

void Allequations::sizeeqsys() {                
    
  nrows=0;ncols=0;nspleq=0;neeq=0;nfeq=0;naddeq=0;nreaeq=0;
  nsmootheq=0;nallequations=0;nfree=0;nelem=0;natom=0;natfit=0;nzl=0;

  // find out numbers of spline-coefficients
  // f and derivatives up to (order-1) are assumed to be continuous
  for (int ipot=0; ipot < vpot.size(); ipot++){
    nspleq    += (vpot[ipot].nknots - 1) * (vpot[ipot].ordspl); 
    nallequations += (vpot[ipot].nknots - 1) * (vpot[ipot].ordspl + 1);
    nfree     += (vpot[ipot].nknots - 1);
  }
  
  // find out number of elements 
  for (int imol=0; imol < vmol.size(); imol++){
    for (int ielem=0; ielem < vmol[imol].nelem; ielem++){
      bool newElem = true;          
      for (int jelem=0; jelem < velem.size(); jelem++){
        if ( velem[jelem] == vmol[imol].elemname[ielem] ) {
          newElem = false;
          break;
        }
      }
      if (newElem)  velem.push_back( vmol[imol].elemname[ielem]  );
    }
  }
  nelem = velem.size(); 
  natfit = nelem - veatom.size();

  // find out total number of atoms, number of energy- and force-equations
  for (int i=0; i< vmol.size(); i++){
    natom += vmol[i].natom;
    if (vmol[i].eweight!=0)  neeq  += 1;
    if (vmol[i].fweight!=0){ 
      if(vmol[i].ncontraineda>0){
        nfeq  += 3 * vmol[i].ncontraineda;
      }else nfeq  += 3 * vmol[i].natom;
    }
  }

  // find out number of additional equations and reaction equations and
  // smoothing equations
     
  naddeq = vadd.size();
  nreaeq = vrea.size();

  // find out number of smoothing equations (one equation per
  // spline-interval (type A) or one equation-neighsmooth (type B))
     
  for (int i=0; i< vpot.size(); i++){
    if (vpot[i].tsmooth != "no" && vpot[i].psmooth != 0) {
      if (vpot[i].tsmooth=="C") nsmootheq += (vpot[i].nknots - vpot[i].neighsmooth) * (vpot[i].derivsmooth + 1);
      else                      nsmootheq +=  vpot[i].nknots - vpot[i].neighsmooth;
    }
  }

  // sum nrows and ncols
  nrows = nspleq + neeq + nfeq + naddeq + nreaeq + nsmootheq;
  ncols = nallequations + natfit;

  if ( ncols > nrows ) {
    nzl    = ncols-nrows;
    nrows += ncols-nrows;
  }

  // size equation system 
  eqmat.resize(nrows,ncols);
  vunknown.resize(ncols);
  vref.resize(nrows);
  vweight.resize(nrows);
  sigma.resize(ncols);
  sigma_e_f.resize(nfree+natfit);
  vres.resize(nrows);
  assembly_cache_ready = false;
  
  // initialize eqmat, vref, vunkown, sigmas to 0, vweight to 1
  for (int i=0; i < nrows; i++){
    vref[i]    = 0;
    vweight[i] = 1;
    vres[i]    = 0;
    for (int j=0; j < ncols; j++){
      eqmat(i,j) = 0;
    }
  }
  for (int j=0; j < ncols; j++) {
    vunknown[j] = 0;
    sigma[j]    = 0;
  }
  for (int i=0; i < nfree+natfit; i++){
    sigma_e_f[i] = 0;
  }
}

int Allequations::include_splineeq(const int ieq_){
// gets number of equation already set up 
// returns total number of equation after this funtion
/////////////////////////////////////////////////////////////////////////////////////////////////////
// f_i(x) = a_i^0 + a_i^1 * (x-x_i) + a_i^2 * (x-x_i)^2 + ...
// example 4th order spline   
// eqmat =
// 1   (x_2 - x_1)    (x_2 - x_1)^2    (x_2 - x_1)^3      (x_2 - x_1)^4   ...  -1    0    0    0   0 ... { f_1(x_2) - f_2(x_2) = 0 }
// 0        1       2*(x_2 - x_1)    3*(x_2 - x_1)^2    4*(x_2 - x_1)^3   ...   0   -1    0    0   0 ... { f'_1(x_2) - f'_2(x_2) = 0 }
// 0        0              2         6*(x_2 - x_1)^2   12*(x_2 - x_1)^2   ...   0    0   -2    0   0 ... { f''_1(x_2) - f''_2(x_2) = 0 }
// 0        0              0                6          24*(x_2 - x_1)     ...   0    0    0   -6   0 ... { f'''_1(x_2) - f'''_2(x_2) = 0 }
// 0        0              0                0                  0          ...   1  (x_3-x_2) ...
//            ...       ...           ...             ...                ...         ...   ........................
//
//            for all inner knots, 
//            than cutoff equation is included { f_last(x_cutoff) = 0 } and derivatives
//            than the equations of the next potential follow...
/////////////////////////////////////////////////////////////////////////////////////////////////////

  int ieq=ieq_; // position variable for row of matrix (current spline eqation)
  int icol=0;   // position variable for column of matrix
  int hlp;

  // for each potential 
  // for all knots: f_i(x_i+1) - f_i+1(x_i+1) = 0 and derivatives
  // for all coefficient within one knot
  // for all derivatives (0 till order-1)
  for (int ipot=0; ipot < vpot.size(); ipot++){
    for (int iknot=0; iknot < vpot[ipot].nknots - 2 ; iknot++){ 
      for (int icoeff=0; icoeff <= vpot[ipot].ordspl ; icoeff++){
        for (int ider=0; ider <= vpot[ipot].ordspl-1; ider++){
          eqmat(ieq+ider,icol) += myfac(icoeff,ider) * pow(vpot[ipot].vr[iknot+1]-vpot[ipot].vr[iknot], icoeff-ider);
        }
        icol++;
      }
      for (int ider=0; ider <= vpot[ipot].ordspl-1; ider++){
        eqmat(ieq+ider,icol+ider) = -myfac(ider,ider);
      }
      ieq += vpot[ipot].ordspl;
    }
    // equation for last knot: f_last(x_cutoff) = 0 and derivatives
    hlp = vpot[ipot].nknots - 1;
    for (int icoeff=0; icoeff <= vpot[ipot].ordspl ; icoeff++){
      for (int ider=0; ider <= vpot[ipot].ordspl-1; ider++){
        eqmat(ieq+ider,icol) += myfac(icoeff,ider) * pow(  vpot[ipot].vr[hlp]-vpot[ipot].vr[hlp-1] , icoeff-ider  );
      }
      icol++;
    }
    ieq += vpot[ipot].ordspl;
  }
  return ieq;
}


int Allequations::include_energyeq(const int ieq_){
////////////////////////////////////////////////////////////////////////////////
//     Ebind       = Etot - sum(Eat) = Eel + Erep - sum(Eat)                  //
// =>  Ebind - Eel = Erep - sum(Eat) = sum_(A!=B)( Erep_AB(R_AB) ) - sum(Eat) //
//     equations added to eqmat and vref respectively                         //
////////////////////////////////////////////////////////////////////////////////
  
  int ieq=ieq_; // position variable for row of matrix (current spline eqation)

  for (int imol=0; imol < vmol.size(); imol++){
    if ( vmol[imol].eweight==0 ) continue;
    // add left hand side of equation system (eqmat)
    add_energy_lhs(ieq,vmol[imol],1);
    // add -1 to coefficient of Eatoms (only those which are getting fitted)
    for (int iat=0; iat < vmol[imol].natom; iat++){
      int j;
      for (j=0; velem[j] != vmol[imol].atomname[iat]; j++){}
      if ( j < veatom.size() ) vref[ieq] += veatom[j];
      else eqmat(ieq, vunknown.size()+j-nelem ) -= 1;
    }
    // add rhs and weight
    vmol[imol].teel=vmol[imol].eel-vref[ieq];
    vref[ieq]    += vmol[imol].ebindmel;
    vweight[ieq]  = vmol[imol].eweight; 
    ieq++;
  }

  return ieq;
}

int Allequations::include_forceeq(const int ieq_){
  if (use_assembly_cache && assembly_cache_ready) {
    // Cached path: pair geometry + pot identity are fixed; only knot
    // segment / basis powers depend on the current vr[].
    for (size_t ic = 0; ic < force_contribs.size(); ic++) {
      const FContrib& c = force_contribs[ic];
      fill_opt_force_pair(c.ieq, c.ipot, c.dist, c.fac);
    }
    for (size_t ia = 0; ia < force_atom_refs.size(); ia++) {
      const ForceAtomRef& a = force_atom_refs[ia];
      vref[a.ieq  ]    += a.vref[0];
      vref[a.ieq+1]    += a.vref[1];
      vref[a.ieq+2]    += a.vref[2];
      vweight[a.ieq]    = a.weight;
      vweight[a.ieq+1]  = a.weight;
      vweight[a.ieq+2]  = a.weight;
    }
    return ieq_ + (int)force_atom_refs.size() * 3;
  }
  return include_forceeq_legacy(ieq_);
}

int Allequations::include_forceeq_legacy(const int ieq_){
///////////////////////////////////////////////////////////////////////////////
//     F_k       = 0   = Fel_k + Frep_k       , k is coord of every atom     //
// =>  Frep_k    = - Fel_k                                                   //
//               = sum_(at1!=at2) f'_at1at2 * (r_k_at1 - r_k_at2) / r_at1at2 //
//     f'_at1at2 = ...                                                       //
//                                                                           //
//     equations added to eqmat and vref respectively                        //
///////////////////////////////////////////////////////////////////////////////
        
  int   ieq=ieq_;
  double dist;
  Colind i;  
  int iconstraineda;
  bool constrained;
  for (int imol=0; imol < vmol.size(); imol++){
    if ( vmol[imol].fweight==0 ) continue;
    for (int at1=0; at1 < vmol[imol].natom; at1++){
      if(vmol[imol].ncontraineda>0){
        constrained=false; 
        for(iconstraineda=0;iconstraineda<vmol[imol].ncontraineda;iconstraineda++){
          if(at1==vmol[imol].contraineda[iconstraineda]) constrained=true;
        }
        if(!constrained) continue;
      }
      for (int at2=0; at2 < vmol[imol].natom; at2++){
        if ( at1==at2 ) continue;
        if ( at1 > at2 )        dist = vmol[imol].dist(at2,at1);
        else if ( at1 < at2 )   dist = vmol[imol].dist(at1,at2); // create a n x n distmatrix??????

        i = findcolind( vmol[imol].atomname[at1] , vmol[imol].atomname[at2] , dist , vmol[imol].name, at1 , at2 );
        if ( i.e == 2 ) continue; // if distance is greater than cutoff
        if ( i.e == 1 ) {     // if potential not selected to optimize
          addextef(ieq,vmol[imol],at1,at2,"force",1);
          continue;
        }
        // add lhs
        for (int xyz=0; xyz<3; xyz++){
          double fac; 
          fac = ( vmol[imol].coord(at1,xyz) - vmol[imol].coord(at2,xyz) ) / dist ;
          for (int j=1; j <= vpot[i.p].ordspl; j++ ) eqmat(ieq+xyz,i.c+j) += -fac * j * pow(  dist - vpot[i.p].vr[i.k-1]  ,  j-1  ); //add force to lhs, not gradient!
        }
      }
      // add rhs and weight
      vref[ieq  ]    += vmol[imol].frefmel(at1,0);
      vref[ieq+1]    += vmol[imol].frefmel(at1,1);
      vref[ieq+2]    += vmol[imol].frefmel(at1,2);
      vweight[ieq]    = vmol[imol].fweight; 
      vweight[ieq+1]  = vmol[imol].fweight; 
      vweight[ieq+2]  = vmol[imol].fweight; 
      ieq += 3;
    }
  }
  return ieq;
}

int Allequations::include_addeq(const int ieq_) {
  int   ieq=ieq_;
  int   iadd,ipot,icol,knot;

  for (iadd=0; iadd < vadd.size(); iadd++){
    // determine icol for the spline coefficient 
    icol = 0;
    for (ipot=0; ipot < vpot.size() ;ipot++) {
      if ( vpot[ipot].potname == vadd[iadd].potname ) break;
      icol += (vpot[ipot].nknots - 1) * (vpot[ipot].ordspl + 1);
    }
    if ( ipot == vpot.size() ) {
      cerr << endl << "ERROR:" << endl << endl
           << "The potential \"" << vadd[iadd].potname << "\" was not asked to optimize, "
           << "thus there cannot be added an additional equation for mentioned potential." 
           << endl << "exit repopt" << endl << endl;
      exit(1);
    }
    for (knot=0; knot < vpot[ipot].vr.size() && vadd[iadd].dist >= vpot[ipot].vr[knot]; knot++ ){ }
    if ( knot == 0 ) {
      cerr << endl << "ERROR: " << endl << endl 
           << vadd[iadd].potname << " distance(bohr,AA):  " << vadd[iadd].dist << "," << vadd[iadd].dist*AA_Bohr
           << " " << vadd[iadd].comment << endl
           << "This distance in this additional equation is smaller than the value "
           << "of the first knot of the respective potential. " 
           << endl << "exit repopt" << endl << endl;
       exit(1);
     }
     if ( knot >= vpot[ipot].vr.size() ){
       cerr << endl << "ERROR: " << endl << endl
            << vadd[iadd].potname << " distance in AA " << vadd[iadd].dist << " " << vadd[iadd].comment << endl
            << "This distance in this additional equation is greater than the cutoff-distance of the respective potential."
            << endl << "exit repopt" << endl << endl;
       exit(1);
    }
    icol += (knot-1) * ( vpot[ipot].ordspl + 1 );
    // icol determined, now add additional equation in eqmat 
    for (int i=0; i <= vpot[ipot].ordspl; i++) {
      if      (vadd[iadd].deriv==0)  eqmat(ieq,icol) +=               pow(  vadd[iadd].dist - vpot[ipot].vr[knot-1]  ,  i    );
      else if (vadd[iadd].deriv==1)  eqmat(ieq,icol) += i*            pow(  vadd[iadd].dist - vpot[ipot].vr[knot-1]  ,  i-1  );
      else if (vadd[iadd].deriv==2)  eqmat(ieq,icol) += i*(i-1)*      pow(  vadd[iadd].dist - vpot[ipot].vr[knot-1]  ,  i-2  );
      else if (vadd[iadd].deriv==3)  eqmat(ieq,icol) += i*(i-1)*(i-2)*pow(  vadd[iadd].dist - vpot[ipot].vr[knot-1]  ,  i-3  );
      icol++;
    }
    vref[ieq]    = vadd[iadd].val;
    vweight[ieq] = vadd[iadd].weight;
    ieq++;
  }

  return ieq;
}

int Allequations::include_reaeq(const int ieq_) {
  if (use_assembly_cache && assembly_cache_ready) {
    int ieq=ieq_;
    for (int irea=0; irea < vrea.size(); irea++ ){
      vrea[irea].treaE=0.0;
      for (int imol=0; imol < vrea[irea].vmol.size(); imol++){
        int ind=vrea[irea].vmol[imol];
        add_energy_lhs_from_contribs(ieq, vreamol_econtrib[ind], vrea[irea].vcoeff[imol]);
        vref[ieq] -= vrea[irea].vcoeff[imol] * vreamol[ind].eel;
        vrea[irea].treaE += vrea[irea].vcoeff[imol] * vreamol[ind].eel;
      }
      vref[ieq]    += vrea[irea].reaE;
      vweight[ieq]  = vrea[irea].reaweight;
      ieq++;
    }
    return ieq;
  }
  return include_reaeq_legacy(ieq_);
}

int Allequations::include_reaeq_legacy(const int ieq_) {
////////////////////////////////////////////////////////////////////////////////////////////////
//  Erea = sum (E_tot^prod) - sum(E_tot^reactants)                                            //
//  e.g.  3 c2h2 --> c6h6                                                                     //
//  Erea = Etot(c6h6) - 3 Etot(c2h2) = Eel(c6h6) + Erep(c6h6) - 3 Eel(c2h2) - 3 Erep(c2h2)    //
//  =>     Erep(c6h6) - 3 Erep(c2h2) = Erea - Eel(c6h6) + 3 Eel(c2h2)                         //
//  now express Erep in Splinecoefficient as described in the energyequation                  //
//                                                                                            //
////////////////////////////////////////////////////////////////////////////////////////////////

  int     ieq=ieq_;
  int ind;
  
  for (int irea=0; irea < vrea.size(); irea++ ){
    vrea[irea].treaE=0.0;
    for (int imol=0; imol < vrea[irea].vmol.size(); imol++){
      ind=vrea[irea].vmol[imol];
      add_energy_lhs_legacy(ieq,vreamol[ind],vrea[irea].vcoeff[imol]); 
      vref[ieq] -= vrea[irea].vcoeff[imol] * vreamol[ind].eel;
      vrea[irea].treaE += vrea[irea].vcoeff[imol] * vreamol[ind].eel;
    }
    vref[ieq]    += vrea[irea].reaE; 
    vweight[ieq]  = vrea[irea].reaweight; 
    ieq++;
  }

  return ieq;
}

int Allequations::include_smootheq(const int ieq_) {
/////////////////////////////////////////////////////////////////////////////////////////////////////
// f_i(x) = a_i^0 + a_i^1 * (x-x_i) + a_i^2 * (x-x_i)^2 + ...                                      //
//                                                                                                 //
// TSMOOTH=A                                                                                       //
// \int_{x_i}^{x_{i+1}} f^{[splord derivative]}_i(x) dx = 0                                        //
// if psmooth is set to infinity the potential turns into the zero-line                            //
// for each knot of each potential one soomthing equation is added as                              //
//                                               (here example 4th order spline):                  //
// eqmat(ieq, ) =  0  0  0  0  24*(x_{i+1}-x_i)*psmooth 0  0  0  0  0 ...  ; vref(ieq) = 0         //
//                                                                                                 //
// TSMOOTH=B                                                                                       //
// d=derivsmooth; j=neighsmooth                                                                    //
// [ f_i^d(x_{i+j}) - f_{i+j}^d(x_{i+j}) ]*psmooth  = 0                                            //
// if psmooth is zero: no smoothing                                                                //
// if psmooth is large: only (j-1) remaining degrees of freedom for the Spline-function            //
//                                                                                                 //
// TSMOOTH=C                                                                                       //
// as B, but conditions for all derivatives 0<d<derivsmooth                                        //
/////////////////////////////////////////////////////////////////////////////////////////////////////
  
  int   ieq=ieq_;
  int   icol=0;
  
  for (int ipot=0; ipot < vpot.size() ; ipot++){
    if (vpot[ipot].tsmooth == "no" || vpot[ipot].psmooth==0) {
      icol += (vpot[ipot].ordspl+1)*(vpot[ipot].nknots-1);
      continue;
    }
    int iord=vpot[ipot].ordspl;

    // TSMOOTH=A
    if (vpot[ipot].tsmooth=="A"){
      for (int iknot=0; iknot < vpot[ipot].nknots - 1 ; iknot++){
        eqmat(ieq,icol+iord) = myfac(iord,iord)*( vpot[ipot].vr[iknot+1] - vpot[ipot].vr[iknot] );
        vweight[ieq] = vpot[ipot].psmooth;
        ieq++;
        icol += vpot[ipot].ordspl+1;
      }
    }

    // TSMOOTH=B and TSMOOTH=C
    else if (vpot[ipot].tsmooth=="B"||vpot[ipot].tsmooth=="C"){
      int iicol   = icol;
      int inknots = vpot[ipot].nknots;
      int inei    = vpot[ipot].neighsmooth;
      int iider;
      if (vpot[ipot].tsmooth=="B")      iider = vpot[ipot].derivsmooth;
      else if (vpot[ipot].tsmooth=="C") iider = 0;
      for (int ider=iider;ider<=vpot[ipot].derivsmooth;ider++){
        iicol=icol;
        for (int iknot=0; iknot < inknots - inei - 1; iknot++){
          for (int ii=ider;ii<=iord;ii++){
            eqmat(ieq,iicol+ii) = myfac(ii,ider)*(pow( vpot[ipot].vr[iknot+inei]-vpot[ipot].vr[iknot] , ii-ider));
          }
          eqmat(ieq,iicol+inei*(iord+1)+ider) = -myfac(ider,ider);
          vweight[ieq] = vpot[ipot].psmooth;
          ieq++;
          iicol += iord+1;
        }
        // last condition
        for (int ii=ider;ii<=iord;ii++){
          eqmat(ieq,iicol+ii) = myfac(ii,ider)*(pow( vpot[ipot].vr[inknots]-vpot[ipot].vr[inknots-inei] , ii-ider));
        }
        vweight[ieq] = vpot[ipot].psmooth;
        ieq++;
        iicol += iord+1;
        // end last condition
      }
      icol += (iord+1)*(inknots-1);
    }

  } // end ipot loop

  return ieq;
}

void Allequations::weighteqsys() {
 
  for (int irow=0;irow<nrows;irow++){
    for (int icol=0;icol<ncols;icol++){
      eqmat(irow,icol) = eqmat(irow,icol)*vweight[irow];
    }
    vref[irow] = vref[irow]*vweight[irow];
  }
}


void Allequations::svd_fulfill_spleq(){
/////////////////////////////////////////////////////////////////////////////////////////////////////
//
// split up eqmat in T, W, E and vunknown in a3 and a4 as follows:
// a3 contains only the unknowns of (order-1) of the unknowns of each spline-segment (vector)
// a4 contains the remaining unknowns (vector)
//
// all spline-equations:
//  T * a3  =  W * a4   ==>  a3 = T^(-1) * W * a4
// remaining equations:
//  E * vunknown = vref'    (' means excluded the spline-equations where vref =0)
//
//  further split-up of E
//  e3 * a3 + e4 * a4 = vref'
//  (e3 * T^(-1) * W + e4) * a4 = vref'
//
// now singular value decomposition for (e3 * T^(-1) * W + e4) = A (as denoted in code)
//      A * a4 = vref' = u * sigma_e_f * vt * a4
//      a4 = v * sigma_e_f^(-1) * ut * vref'
//  
//  calculation of a4 done as follows:
//      a4[i] = v * u_coli * vref' / sigma_e_f(i)
//
// resubsitute to yield a3
//  a3 = T^(-1) * W * a4    
// 
// write a3 and a4 in vunknown
//
// calculate exponential coefficient for every potential (vunknown needed) 
/////////////////////////////////////////////////////////////////////////////////////////////////////
  
  Timer ts;

  //cout << "Ts0 taken: " << ts.elapsed() << " seconds\n"; ts.reset();
  int a_svd = 4;
  // a_svd is which parameter of each ploynomial is the free potential dependent parameter
  // since it has almost no influence we use a fixed coding...

  //cout<<"eqmat\n"<<eqmat<<endl;

  MatrixXd T = MatrixXd::Zero(nspleq,nspleq);
  MatrixXd W = MatrixXd::Zero(nspleq,nfree);
  //cout << "Ts1 taken: " << ts.elapsed() << " seconds\n"; ts.reset();

  for (int ieq=0; ieq < nspleq; ieq++){
    int icoleqmat=0;
    int icolT=0;
    int icolW=0;
    for (int ipot=0; ipot<vpot.size(); ipot++){
      for (int iknot=0; iknot < vpot[ipot].nknots - 1 ; iknot++){
        for (int ia=0; ia < a_svd; ia++){
          T(ieq,icolT) = eqmat(ieq,icoleqmat+ia);
          icolT++;
        }
        W(ieq,icolW) = -eqmat(ieq,icoleqmat+a_svd);
        icolW++;
        for (int ia=a_svd+1; ia <= vpot[ipot].ordspl; ia++){
          T(ieq,icolT) = eqmat(ieq,icoleqmat+ia);
          icolT++;
        }
        icoleqmat += vpot[ipot].ordspl+1 ;
      }
    }
  }
  //cout << "Ts2 taken: " << ts.elapsed() << " seconds\n"; ts.reset();
  //cout<<"T\n"<<T<<endl;
  //cout<<"W\n"<<W<<endl;

  MatrixXd TinvW = MatrixXd::Zero(nspleq,nfree);
  if(idecompose==1){
    TinvW = T.ldlt().solve(W);
  }else if(idecompose==2){
    TinvW = T.partialPivLu().solve(W);
  }else if(idecompose==3){
    TinvW = T.fullPivLu().solve(W);
  }else if(idecompose==4){
    TinvW = T.householderQr().solve(W);
  }else if(idecompose==5){
    TinvW = T.colPivHouseholderQr().solve(W);
  }else if(idecompose==6){
    TinvW = T.fullPivHouseholderQr().solve(W);
  }else if(idecompose==7){
    TinvW = T.completeOrthogonalDecomposition().solve(W);
  }
  //cout<<"TinvW\n"<<TinvW<<endl;
  //double relative_error = (T*TinvW - W).norm() / W.norm(); // norm() is L2 norm
  ////cout << "The relative error is:\n" << relative_error << endl;
  //cout << "Ts3 taken: " << ts.elapsed() << " seconds\n"; ts.reset();

  // now apply to all other eqations but the splineequations
  // first seperate energyandforce - matrix to e3 and e4
  MatrixXd e3 = MatrixXd::Zero( eqmat.rows()-nspleq, nallequations-nfree );  
  MatrixXd e4 = MatrixXd::Zero( eqmat.rows()-nspleq, nfree );  
  //cout << "Ts4 taken: " << ts.elapsed() << " seconds\n"; ts.reset();

  for (int ieq=nspleq; ieq < eqmat.rows(); ieq++){
    int icoleqmat=0;
    int icole3=0;
    int icole4=0;
    for (int ipot=0; ipot<vpot.size(); ipot++){
      for (int iknot=0; iknot < vpot[ipot].nknots - 1 ; iknot++){
        for (int ia=0; ia < a_svd; ia++){
          e3(ieq-nspleq,icole3) = eqmat(ieq,icoleqmat+ia);
          icole3++;
        }
        e4(ieq-nspleq,icole4) = eqmat(ieq,icoleqmat+a_svd);
        icole4++;
        for (int ia=a_svd+1; ia <= vpot[ipot].ordspl; ia++){
          e3(ieq-nspleq,icole3) = eqmat(ieq,icoleqmat+ia);
          icole3++;
        }
        icoleqmat += vpot[ipot].ordspl+1;
      }
    }
  }
  //cout << "Ts5 taken: " << ts.elapsed() << " seconds\n"; ts.reset();

  //cout<<"e3\n"<<e3<<endl;
  //cout<<"e4\n"<<e4<<endl;
  // Construct A 
  MatrixXd A = MatrixXd::Zero(eqmat.rows()-nspleq, nfree); 
  //cout << "Ts6 taken: " << ts.elapsed() << " seconds\n"; ts.reset();

  A = e3*TinvW + e4;
  //cout << "Ts7 taken: " << ts.elapsed() << " seconds\n"; ts.reset();

  // now add e_atom contribution
  //cout<<"A0\n"<<A<<endl;
  //A.resize(A.rows(), A.cols() + natfit);
  A.conservativeResize(A.rows(), A.cols() + natfit);
  for (int iel=-natfit; iel<0; iel++){
    for (int irow=0; irow < A.rows(); irow++){
      A(irow,A.cols()+iel) = eqmat(nspleq+irow,eqmat.cols()+iel);
    }
  }
  //cout<<"A\n"<<A<<endl;
  //cout << "Ts8 taken: " << ts.elapsed() << " seconds\n"; ts.reset();

  // then svd for (e3 * w_in2 + e4)
  //MatrixXd u ( eqmat.rows()-nspleq , nfree+natfit ); // we do not need all u(rowdim,rowdim) 
  //MatrixXd vt( nfree+natfit         , nfree+natfit );
  //BDCSVD<MatrixXd> svdofA(A, ComputeThinU | ComputeThinV);
  //sigma_e_f = svdofA.singularValues();
  //u = svdofA.matrixU();
  //vt = svdofA.matrixV();

  VectorXd vref_(vref.size()-nspleq);
  
  for (int i=0; i < vref_.size(); i++ ){
    vref_[i] = vref[i+nspleq];
  }

  VectorXd a4_el (nfree + natfit);
  VectorXd a4    (nfree);
  VectorXd a3    (nallequations-nfree);
  //cout << "Ts9 taken: " << ts.elapsed() << " seconds\n"; ts.reset();
  //cout<<"a3\n"<<a3<<endl;
  //cout<<"a4\n"<<a4<<endl;
  //cout<<"a4_el\n"<<a4_el<<endl;
  //cout<<"vref\n"<<vref<<endl;
  //cout<<"vref_\n"<<vref_<<endl;

  if(ilmsfit==1){
    a4_el = A.householderQr().solve(vref_);
  }else if(ilmsfit==2){
    a4_el = A.colPivHouseholderQr().solve(vref_);
  }else if(ilmsfit==3){
    a4_el = A.fullPivHouseholderQr().solve(vref_);
  }else if(ilmsfit==4){
    a4_el = A.bdcSvd(ComputeThinU | ComputeThinV).solve(vref_);
  }
  //cout << "Ts10 taken: " << ts.elapsed() << " seconds\n"; ts.reset();

  for (int i=0; i < a4_el.size()-natfit; i++){
    a4[i] = a4_el[i];
  }         
  a3  = TinvW*a4;
  //cout << "Ts11 taken: " << ts.elapsed() << " seconds\n"; ts.reset();

  // tested again beginning here
  // write a3 and a4 in vunknown 
  int icolunknown=0;
  int icola3=0;
  int icola4=0;
  for (int ipot=0; ipot<vpot.size(); ipot++){
    for (int iknot=0; iknot < vpot[ipot].nknots - 1 ; iknot++){
      for (int ia=0; ia < a_svd; ia++){
        vunknown[icolunknown+ia] = a3[icola3];
        icola3++;
      }
      vunknown[icolunknown+a_svd] = a4_el[icola4];
      icola4++;
      //for (int ia=a_svd+1; ia <= vpot[ipot].ordspl; ia++){
      //  vunknown[icolunknown+ia] = a3[icola3];
      //  icola3++;
      //}
      icolunknown += vpot[ipot].ordspl+1 ;
    }
  }
  for (int iel=-natfit; iel<0; iel++){
    vunknown[vunknown.size()+iel] = a4_el[a4_el.size()+iel];
  }
  //cout << "Ts12 taken: " << ts.elapsed() << " seconds\n"; ts.reset();

  // calculate exponential for each potential ( Erep=exp(-Ar+B)+C )
  int icoeff=0;
  for (int i=0; i < vpot.size(); i++){
    vpot[i].expA = -2*vunknown[icoeff+2]/vunknown[icoeff+1];
    vpot[i].expB = log(-vunknown[icoeff+1]/vpot[i].expA) + vpot[i].expA*vpot[i].vr[0];
    vpot[i].expC = vunknown[icoeff] - exp(-vpot[i].expA*vpot[i].vr[0]+vpot[i].expB);
    icoeff += ( vpot[i].vr.size()-1 ) * ( vpot[i].ordspl+1 );
  }
  //cout << "Ts13 taken: " << ts.elapsed() << " seconds\n"; ts.reset();
  //cout <<"vunknown\n"<<vunknown<<endl;
}

void Allequations::get_residual() {

  double tmp;
  vres = eqmat*vunknown - vref;
  
  restotS=0;reseS=0;resfS=0;resaddS=0;resreaS=0;ressmoothS=0;
  restotU=0;reseU=0;resfU=0;resaddU=0;resreaU=0;ressmoothU=0;
  restot2=0;rese2=0;resf2=0;resadd2=0;resrea2=0;ressmooth2=0;
  restot4=0;rese4=0;resf4=0;resadd4=0;resrea4=0;ressmooth4=0;
  restot8=0;rese8=0;resf8=0;resadd8=0;resrea8=0;ressmooth8=0;

  int irow=nspleq;
  for (int i=0;i<neeq;i++){
    reseS += vres[irow];
    reseU += abs(vres[irow]);
    tmp    = vres[irow]*vres[irow];
    rese2 += tmp;//vres[irow]*vres[irow];
    rese4 += tmp*tmp;
    rese8 += tmp*tmp*tmp*tmp;
    irow++;
  }
  for (int i=0;i<nfeq;i++){
    resfS += vres[irow];
    resfU += abs(vres[irow]);
    tmp    = vres[irow]*vres[irow];
    resf2 += tmp;//vres[irow]*vres[irow];
    resf4 += tmp*tmp;
    resf8 += tmp*tmp*tmp*tmp;
    irow++;
  }
  for (int i=0;i<naddeq;i++){
    resaddS += vres[irow];
    resaddU += abs(vres[irow]);
    tmp      = vres[irow]*vres[irow];
    resadd2 += tmp;//vres[irow]*vres[irow];
    resadd4 += tmp*tmp;
    resadd8 += tmp*tmp*tmp*tmp;
    irow++;
  }
  for (int i=0;i<nreaeq;i++){
    resreaS += vres[irow];
    resreaU += abs(vres[irow]);
    tmp      = vres[irow]*vres[irow];
    resrea2 += tmp;//vres[irow]*vres[irow];
    resrea4 += tmp*tmp;
    resrea8 += tmp*tmp*tmp*tmp;
    irow++;
  }
  for (int i=0;i<nsmootheq;i++){
    ressmoothS += vres[irow];
    ressmoothU += abs(vres[irow]);
    tmp         = vres[irow]*vres[irow];
    ressmooth2 += tmp;//vres[irow]*vres[irow];
    ressmooth4 += tmp*tmp;
    ressmooth8 += tmp*tmp*tmp*tmp;
    irow++;
  }
  restotS = reseS +resfS +resaddS +resreaS +ressmoothS ;
  restotU = reseU +resfU +resaddU +resreaU +ressmoothU ;
  restot2 = rese2+resf2+resadd2+resrea2+ressmooth2;
  restot4 = rese4+resf4+resadd4+resrea4+ressmooth4;
  restot8 = rese8+resf8+resadd8+resrea8+ressmooth8;
}

////// OLD FUNCTION NOT UPDATED TO ACTUAL PROGRAMM!!!!!!!!!!!!!!!!
void Allequations::svd_all() {
// singular value decomposition of eqmat and following determination of unknowns
//
// eqmat = A = u * sigma * vt    (note: ut = u^(-1) and vt = v^(-1), t stands for transposed)
// vunknown = x
//
// A * x = vref = u * sigma * vt * x
//     x = v * sigma^(-1) * ut * vref
// 
//  x[i] = v * u_coli * b / sigma(i)
//
// truncation of sigma, in order not to lose information

  int rowdim = eqmat.rows();
  int coldim = eqmat.cols();

  MatrixXd A = eqmat;  

  // singular value decomposition
  //JacobiSVD<MatrixXd> svdofA(A, ComputeThinU | ComputeThinV);
  BDCSVD<MatrixXd> svdofA(A, ComputeThinU | ComputeThinV);
  vunknown = svdofA.solve(vref);
}

double Allequations::myfac(const int a, const int b){
// returns coefficient z=a*(a-1)*(a-2)*... for the b'th derivative
// (x^a)^[b] = a*(x^{a-1})^[b-1] = ... = z*x^{a-b}
  if (b==0)     return 1; 
  else if (b>a) return 0;
  else          return a*myfac(a-1,b-1);
}


void Allequations::add_energy_lhs(const int ieq, const Molecule& mol, const double coeff){
  if (use_assembly_cache && assembly_cache_ready) {
    // Prefer table lookup by molecule address within vmol / vreamol.
    for (size_t im = 0; im < vmol.size(); im++) {
      if (&vmol[im] == &mol) {
        add_energy_lhs_from_contribs(ieq, vmol_econtrib[im], coeff);
        return;
      }
    }
    for (size_t im = 0; im < vreamol.size(); im++) {
      if (&vreamol[im] == &mol) {
        add_energy_lhs_from_contribs(ieq, vreamol_econtrib[im], coeff);
        return;
      }
    }
  }
  add_energy_lhs_legacy(ieq, mol, coeff);
}

void Allequations::add_energy_lhs_legacy(const int ieq, const Molecule& mol, const double coeff){

  Colind i; 
  for (int at1=0; at1 < mol.natom-1; at1++){
    for (int at2=at1+1; at2 < mol.natom; at2++){
      i = findcolind(mol.atomname[at1],mol.atomname[at2],mol.dist(at1,at2),mol.name,at1,at2);
      if ( i.e == 2 ) continue; // if distance is greater than cutoff
      if ( i.e == 1 ){          // if pot is not selected to be optimized
        addextef(ieq,mol,at1,at2,"energy",coeff);
        continue; 
      }
      for (int j=0; j <= vpot[i.p].ordspl; j++) {
        eqmat(ieq,i.c) += coeff * pow(  mol.dist(at1,at2) - vpot[i.p].vr[i.k-1]  ,  j  );
        i.c++;
      } 
    }
  }
}


Colind Allequations::findcolind(const string  atomname1, const string atomname2, const double dist, 
                                const string& molname, const int at1, const int at2 ) {
  Colind i; i.c=0; i.p=0; i.k=0; i.e=0;
  string potname;

  potname = findpotname( atomname1 , atomname2 ); 
  for (i.p=0; i.p < vpot.size() ;i.p++) {
    if ( vpot[i.p].potname == potname ) break;        
    i.c += (vpot[i.p].nknots - 1) * (vpot[i.p].ordspl + 1);
  }
  if ( i.p >= vpot.size() ) {     // if potential is not selected to be optimized
    excludepot( potname , dist , molname);        
    i.e=1; 
  }
  else {
    for (i.k=0; i.k < vpot[i.p].vr.size() && dist >= vpot[i.p].vr[i.k]; i.k++ ){ }  
    if ( i.k == 0 ) {
      cerr  << endl << "ERROR: " << endl << endl
      << atomname1 << at1+1           
      << atomname2 << at2+1 << "    " << dist*AA_Bohr << " AA   " 
      << molname << endl
      << "This distance in this molecule is smaller than the value "
      << "of the first knot of the respective potential. " << endl
      << "Solution: rerun with first knot shifted to a smaller value"
      << endl << "exit repopt" << endl << endl;
      exit(1);
    }
    if ( i.k >= vpot[i.p].vr.size() ) i.e=2;      // if distance is greater than cutoff
    else i.c += (i.k-1) * ( vpot[i.p].ordspl + 1 );
  }
  return i; 
}

string Allequations::findpotname(const string at1, const string at2) const{
// returns potname if it exists otherwise returns a string = at1+at2
  string potname, at1at2, at2at1;
  bool   potEXIST=false;

  at1at2 = at1 + at2;
  at2at1 = at2 + at1;
  for (int i=0; i<vpot.size() ;i++){
    if ( vpot[i].potname == at1at2 ){
      potEXIST = true;
      potname  = at1at2;
      break;
    }
    if ( vpot[i].potname == at2at1 ) {
      potEXIST = true;
      potname  = at2at1;
      break;
    }
  }
  if (potEXIST)
    return potname;
  else
    return at1at2;
}

void Allequations::excludepot(const string potname, const double dist, const string& molname){
// keeps track of excluded potentials due to double occurence e.g. c_h_ and h_c_
  bool newExclPot=true;
  int ipot;
  int size=vexclpot.size();

  for (ipot=0; ipot<size; ipot++){
    if (potname == vexclpot[ipot].potname){
      newExclPot = false;
      break;
    }
  }
  if ( newExclPot ) {
    vexclpot.resize(size+1);
    vexclpot[size].potname = potname;
    vexclpot[size].ndist   = 1;
    vexclpot[size].mindist = dist;
    vexclpot[size].molname = molname;
  }
  else {
    vexclpot[ipot].ndist++;
    if ( dist < vexclpot[ipot].mindist ){
      vexclpot[ipot].mindist = dist;
      vexclpot[ipot].molname = molname;
    }
  }
}

void Allequations::ifnotfile(const string filename) const{
  ifstream file;
  file.open(filename.c_str());
  if ( !file ) {
    cerr << endl << "ERROR: " << endl << endl;
    cerr << filename << ": file not found" 
         << endl << "exit repopt" << endl << endl;
    exit(1);
  }
  file.close();
}

void Allequations::addextef(const int ieq,const Molecule& mol,const int at1, 
      const int at2,const string ef,const double coeff){
  return;
}

int Allequations::find_pot_index(const string& potname) const {
  for (int ipot = 0; ipot < (int)vpot.size(); ipot++) {
    if (vpot[ipot].potname == potname) return ipot;
  }
  return -1;
}

// Same segment scan as findcolind (dist >= vr[k] advances k). Returns k.
// k==0 => below first knot (fatal); k>=nknots => beyond cutoff.
int Allequations::knot_segment_index(const int ipot, const double dist,
                                     const string& atomname1, const string& atomname2,
                                     const string& molname, const int at1, const int at2) const {
  int k = 0;
  for (; k < (int)vpot[ipot].vr.size() && dist >= vpot[ipot].vr[k]; k++) {}
  if (k == 0) {
    cerr << endl << "ERROR: " << endl << endl
         << atomname1 << at1+1
         << atomname2 << at2+1 << "    " << dist*AA_Bohr << " AA   "
         << molname << endl
         << "This distance in this molecule is smaller than the value "
         << "of the first knot of the respective potential. " << endl
         << "Solution: rerun with first knot shifted to a smaller value"
         << endl << "exit repopt" << endl << endl;
    exit(1);
  }
  return k;
}

int Allequations::fill_opt_energy_pair(const int ieq, const int ipot, const double dist, const double coeff) {
  const int k = knot_segment_index(ipot, dist, "?", "?", "cached", 0, 0);
  if (k >= (int)vpot[ipot].vr.size()) return 0; // beyond cutoff
  int ic = pot_col0[ipot] + (k - 1) * (vpot[ipot].ordspl + 1);
  for (int j = 0; j <= vpot[ipot].ordspl; j++) {
    eqmat(ieq, ic) += coeff * pow(dist - vpot[ipot].vr[k - 1], j);
    ic++;
  }
  return 1;
}

int Allequations::fill_opt_force_pair(const int ieq, const int ipot, const double dist, const double fac[3]) {
  // Use dummy names only for the fatal below-first-knot path; cached pairs
  // were valid at prepare time so this should not trigger in normal GA use.
  const int k = knot_segment_index(ipot, dist, "?", "?", "cached", 0, 0);
  if (k >= (int)vpot[ipot].vr.size()) return 0;
  const int ic0 = pot_col0[ipot] + (k - 1) * (vpot[ipot].ordspl + 1);
  for (int xyz = 0; xyz < 3; xyz++) {
    for (int j = 1; j <= vpot[ipot].ordspl; j++) {
      eqmat(ieq + xyz, ic0 + j) += -fac[xyz] * j * pow(dist - vpot[ipot].vr[k - 1], j - 1);
    }
  }
  return 1;
}

void Allequations::add_energy_lhs_from_contribs(const int ieq, const vector<EContrib>& contribs, const double coeff) {
  for (size_t i = 0; i < contribs.size(); i++) {
    fill_opt_energy_pair(ieq, contribs[i].ipot, contribs[i].dist, coeff);
  }
}

void Allequations::collect_energy_contribs(const Molecule& mol, vector<EContrib>& out) {
  out.clear();
  for (int at1 = 0; at1 < mol.natom - 1; at1++) {
    for (int at2 = at1 + 1; at2 < mol.natom; at2++) {
      const double dist = mol.dist(at1, at2);
      const string potname = findpotname(mol.atomname[at1], mol.atomname[at2]);
      const int ipot = find_pot_index(potname);
      if (ipot < 0) {
        // Non-optimized pot: legacy only calls addextef (currently a no-op)
        // and excludepot for bookkeeping. Record exclusion once here.
        excludepot(potname, dist, mol.name);
        continue;
      }
      // Keep pairs even if currently beyond cutoff: cutoff knot can move in GA.
      EContrib c;
      c.ipot = ipot;
      c.dist = dist;
      out.push_back(c);
    }
  }
}

void Allequations::build_assembly_cache() {
  pot_col0.assign(vpot.size(), 0);
  int col = 0;
  for (int ipot = 0; ipot < (int)vpot.size(); ipot++) {
    pot_col0[ipot] = col;
    col += (vpot[ipot].nknots - 1) * (vpot[ipot].ordspl + 1);
  }

  vmol_econtrib.assign(vmol.size(), vector<EContrib>());
  for (size_t im = 0; im < vmol.size(); im++) {
    collect_energy_contribs(vmol[im], vmol_econtrib[im]);
  }

  vreamol_econtrib.assign(vreamol.size(), vector<EContrib>());
  for (size_t im = 0; im < vreamol.size(); im++) {
    collect_energy_contribs(vreamol[im], vreamol_econtrib[im]);
  }

  force_contribs.clear();
  force_atom_refs.clear();
  // Absolute force-row index matches include_forceeq after spline+energy blocks.
  int ieq = nspleq + neeq;
  for (int imol = 0; imol < (int)vmol.size(); imol++) {
    if (vmol[imol].fweight == 0) continue;
    for (int at1 = 0; at1 < vmol[imol].natom; at1++) {
      if (vmol[imol].ncontraineda > 0) {
        bool constrained = false;
        for (int iconstraineda = 0; iconstraineda < vmol[imol].ncontraineda; iconstraineda++) {
          if (at1 == vmol[imol].contraineda[iconstraineda]) constrained = true;
        }
        if (!constrained) continue;
      }
      for (int at2 = 0; at2 < vmol[imol].natom; at2++) {
        if (at1 == at2) continue;
        double dist;
        if (at1 > at2)      dist = vmol[imol].dist(at2, at1);
        else                dist = vmol[imol].dist(at1, at2);
        const string potname = findpotname(vmol[imol].atomname[at1], vmol[imol].atomname[at2]);
        const int ipot = find_pot_index(potname);
        if (ipot < 0) {
          excludepot(potname, dist, vmol[imol].name);
          // addextef is a no-op
          continue;
        }
        FContrib c;
        c.ieq = ieq;
        c.ipot = ipot;
        c.dist = dist;
        for (int xyz = 0; xyz < 3; xyz++) {
          c.fac[xyz] = (vmol[imol].coord(at1, xyz) - vmol[imol].coord(at2, xyz)) / dist;
        }
        force_contribs.push_back(c);
      }
      ForceAtomRef a;
      a.ieq = ieq;
      a.vref[0] = vmol[imol].frefmel(at1, 0);
      a.vref[1] = vmol[imol].frefmel(at1, 1);
      a.vref[2] = vmol[imol].frefmel(at1, 2);
      a.weight = vmol[imol].fweight;
      force_atom_refs.push_back(a);
      ieq += 3;
    }
  }

  size_t n_vmol_e = 0, n_vrea_e = 0;
  for (size_t i = 0; i < vmol_econtrib.size(); i++) n_vmol_e += vmol_econtrib[i].size();
  for (size_t i = 0; i < vreamol_econtrib.size(); i++) n_vrea_e += vreamol_econtrib[i].size();
  assembly_cache_ready = true;
  cout << "#assembly_cache: force_pairs=" << force_contribs.size()
       << " force_atoms=" << force_atom_refs.size()
       << " vmol_energy_pairs=" << n_vmol_e
       << " vreamol_energy_pairs=" << n_vrea_e
       << endl;
}

void Allequations::assemble_equations() {
  reset();
  int ieq = 0;
  ieq = include_splineeq(0);
  ieq = include_energyeq(ieq);
  ieq = include_forceeq(ieq);
  ieq = include_addeq(ieq);
  ieq = include_reaeq(ieq);
  ieq = include_smootheq(ieq);
  (void)ieq;
}

bool Allequations::verify_assembly(int ntrial) {
  if (!assembly_cache_ready) {
    cerr << "verify_assembly: cache not built" << endl;
    return false;
  }

  // Save knot vectors
  vector<vector<double> > vr0(vpot.size());
  for (size_t ip = 0; ip < vpot.size(); ip++) vr0[ip] = vpot[ip].vr;

  auto restore = [&]() {
    for (size_t ip = 0; ip < vpot.size(); ip++) vpot[ip].vr = vr0[ip];
  };

  // Deterministic knot perturbations: shift inner knots by small fractions
  // of the local spacing, staying strictly inside (vr[0], vr[nk-1]).
  for (int it = 0; it < ntrial; it++) {
    for (size_t ip = 0; ip < vpot.size(); ip++) {
      vpot[ip].vr = vr0[ip];
      const int nk = vpot[ip].nknots;
      for (int j = 1; j < nk - 1; j++) {
        const double lo = vpot[ip].vr[j - 1];
        const double hi = (j + 1 < nk) ? vr0[ip][j + 1] : vr0[ip][nk - 1];
        const double mid = 0.5 * (lo + hi);
        const double amp = 0.15 * (hi - lo);
        // alternate signs / magnitudes by trial index
        const double phase = ((it + j) % 2 == 0) ? 1.0 : -1.0;
        double nj = mid + phase * amp * (0.2 + 0.1 * (it % 5));
        if (nj <= lo) nj = lo + 1e-8;
        if (nj >= hi) nj = hi - 1e-8;
        vpot[ip].vr[j] = nj;
      }
      // lightly nudge cutoff within a tiny band (keep ascending)
      if (nk >= 2) {
        const double prev = vpot[ip].vr[nk - 2];
        double cut = vr0[ip][nk - 1] + (it % 3 - 1) * 1e-4;
        if (cut <= prev) cut = prev + 1e-8;
        vpot[ip].vr[nk - 1] = cut;
      }
    }

    // Cached assembly
    use_assembly_cache = true;
    assemble_equations();
    MatrixXd eq_cached = eqmat;
    VectorXd ref_cached = vref;
    VectorXd w_cached = vweight;

    // Legacy assembly
    use_assembly_cache = false;
    assemble_equations();
    MatrixXd eq_legacy = eqmat;
    VectorXd ref_legacy = vref;
    VectorXd w_legacy = vweight;
    use_assembly_cache = true;

    double max_eq = 0.0, max_ref = 0.0, max_w = 0.0;
    for (int r = 0; r < nrows; r++) {
      max_ref = std::max(max_ref, std::abs(ref_cached[r] - ref_legacy[r]));
      max_w   = std::max(max_w,   std::abs(w_cached[r] - w_legacy[r]));
      for (int c = 0; c < ncols; c++) {
        max_eq = std::max(max_eq, std::abs(eq_cached(r, c) - eq_legacy(r, c)));
      }
    }

    // Exact floating-point match expected: same pow() order and arithmetic.
    const double tol = 0.0;
    if (max_eq > tol || max_ref > tol || max_w > tol) {
      cerr << "verify_assembly FAIL trial=" << it
           << " max|eqmat|=" << max_eq
           << " max|vref|=" << max_ref
           << " max|vweight|=" << max_w << endl;
      restore();
      return false;
    }
    cout << "verify_assembly OK trial=" << it
         << " max|eqmat|=" << max_eq
         << " max|vref|=" << max_ref
         << " max|vweight|=" << max_w << endl;
  }

  restore();
  return true;
}

void Allequations::time_assembly(int nrep) {
  if (!assembly_cache_ready) {
    cerr << "time_assembly: cache not built" << endl;
    return;
  }
  if (nrep < 1) nrep = 1;

  Timer t;
  use_assembly_cache = true;
  t.reset();
  for (int i = 0; i < nrep; i++) assemble_equations();
  const double t_cached = t.elapsed();

  use_assembly_cache = false;
  t.reset();
  for (int i = 0; i < nrep; i++) assemble_equations();
  const double t_legacy = t.elapsed();
  use_assembly_cache = true;

  cout << fixed << setprecision(6);
  cout << "#time_assembly nrep=" << nrep
       << " cached_total=" << t_cached
       << " legacy_total=" << t_legacy
       << " cached_per=" << (t_cached / nrep)
       << " legacy_per=" << (t_legacy / nrep)
       << " speedup=" << (t_legacy / t_cached)
       << endl;
}



