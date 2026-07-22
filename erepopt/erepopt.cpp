#include <iostream>
#include <vector>
#include <string>
#include <sstream> 
#include <omp.h> 
#include "erepobj.hpp"
#include <fstream>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <unistd.h>
#include <ga/ga.h>
#include <ga/std_stream.h>
#include "tools.hpp"
#include "ga.hpp"
#include "ddh.hpp"
#include "genes.hpp"
#include "mpi.h"

using namespace std;

double MyObjective(GAGenome& gnome);
double MyComparator(const GAGenome&, const GAGenome&);
void  MyInitializer(GAGenome&);
int   MyMutator(GAGenome&, double);
int   MyOnePointCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);
int   MyTwoPointCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);
int   MyUniformCrossover(const GAGenome&, const GAGenome&, GAGenome*, GAGenome*);

vector<string> addHamiltonian;
string popfinalfile="pop.final.dat";
string popinitialfile="pop.initial.dat";
string scratch="/tmp";
int cpu_number=1;
int power=2;
int ga_popsize=100, ga_ngen=300, ga_scoref=1, ga_flushf=1;
int preserved_num=1, seed=0;
double ga_pmut=0.2, ga_pcross=0.9;
double gtol=0.000001;
bool ga=false,readr=true,restart=false;
bool use_uniform_crossover=false;
// Slice crossover breakpoint stride. Default 3 restores legacy "3*RandomInt" alignment.
// Set to 0 (or negative) for unaligned sites on [0, length].
int crossover_align=3;
bool runtest=false,skfclean=true;
bool endgen=false,initialgen=true,firsteval=true;
int idrefr=3;
double s1=1.17,s2=1.40,s3=1.20,s4=2.0,s5=2.0,s6=2.0,s7=2.0,s8=2.0,s9=2.0,sdenswf=2.0;
fstream infile,outfile,restartfile;
sddh ddh;
Erepobj erepobj;
int main(int argc, char** argv){
  int i,j,k,idx,length,rank;
  MPI_Init(NULL,NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

 if(rank==0){
  if ( argc < 2 ){
    cerr << "usage: erepobj inputfile" << endl;
    exit(1);
  }
 }
  ddh.td3=ddh.tdamph=ddh.thubbardderivs=ddh.thirdorderfull=ddh.tdamphver2=ddh.thubbards=ddh.tvorbes=false;
  
  erepobj.prepare(argv[1]);
 if(rank==0){
  erepobj.fout.open(erepobj.outfilename.c_str());
  erepobj.fout<<std::fixed;
 }
  if(ga && !runtest){
    GARandomSeed(seed);
    length=genome_length(erepobj, ddh);
    if(length==0){
      cerr<<"ERROR: no free parameters in GA chromosome (genome length 0)"<<endl;
      exit(1);
    }

    GA1DArrayGenome<double> genome(length, MyObjective);
    genome.initializer(::MyInitializer);
    genome.mutator(::MyMutator);
    genome.comparator(::MyComparator);
    MyGASimpleGA ga(genome);
    if(use_uniform_crossover) ga.crossover(::MyUniformCrossover);
    else ga.crossover(::MyTwoPointCrossover);
  //ga.maximize();
    ga.minimize();
  //GANoScaling scale;
    GASigmaTruncationScaling scale;
    ga.scaling(scale);
    GARouletteWheelSelector select;
    ga.selector(select);
    ga.populationSize(ga_popsize);
    ga.nGenerations(ga_ngen);
    ga.pMutation(ga_pmut);
    ga.pCrossover(ga_pcross);
    ga.scoreFilename("score.dat");  // name of file for scores
    ga.scoreFrequency(ga_scoref); // keep the scores of every 10th generation
    ga.flushFrequency(ga_flushf); // specify how often to write the score to disk
    ga.selectScores(GAStatistics::Minimum);
  //ga.parameters(argc, argv, gaTrue); // parse commands, complain if bogus args
    
 if(rank==0){
    erepobj.fout << "#initializing...\n"; erepobj.fout.flush();
 }

  //string call="rm -rf skftmp; mkdir -p skftmp;";
  //i=system(call.c_str());	
  
    if(restart==true) restartfile.open(popfinalfile.c_str(), ios::in);
    ga.initialize();
    if(restart==true) restartfile.close();
    initialgen=false;
    ga.initialstep();
   if(rank==0){
    outfile.open(popinitialfile.c_str(), (STD_IOS_OUT | STD_IOS_TRUNC));
    for(i=0; i<ga.population().size(); i++){
      genome = ga.population().individual(i);
      write_free_genes(outfile, genome, erepobj, ddh);
      outfile.precision(9);
      outfile.width(20);
      outfile << genome.score() << "\n";
    }
    outfile.close();

    erepobj.fout << "#evolving...\n"; erepobj.fout.flush();
  }
    while(!ga.done()) {
     if(rank==0){
    //cout<<"generation:  "<<ga.generation()<<" start\n"; cout.flush();
      erepobj.fout.precision(2);
      erepobj.fout<<double(ga.generation())/double(ga_ngen)<<"   "; erepobj.fout.flush();
     }
      ga.step();
     if(rank==0){
      ga.flushScores();
      erepobj.fout.precision(9);
//    erepobj.fout<<ga.population().individual(0).score()<<endl; cout.flush();
      erepobj.fout<<ga.population().best().score()<<endl; cout.flush();
 
      outfile.open(popfinalfile.c_str(), (STD_IOS_OUT | STD_IOS_TRUNC));
      for(i=0; i<ga.population().size(); i++){
        genome = ga.population().individual(i);
        write_free_genes(outfile, genome, erepobj, ddh);
        outfile.precision(9);
        outfile.width(20);
        outfile << genome.score() << "\n";
      }
      outfile.close();
      if(skfclean){
        string call="rm -f "+erepobj.libdir+"/*";
        i=system(call.c_str());	
      }
     }
    }
   if(rank==0){
    ga.flushScores();
    erepobj.fout<<"#ga end!\n";
    genome = ga.population().best();
    endgen=true;
    erepobj.get_residual(genome); 
    apply_genome_to_params(genome, erepobj, ddh);
   }
  }

 if(rank==0){
  erepobj.writeout();

  erepobj.fout << "** erepobj normal termination **" << endl << endl;
  erepobj.fout.close();
 }
  MPI_Finalize();
  return 0;
}



double MyObjective(GAGenome& g) {
  GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
  int ie1,ie2,i,j,k;
  double tscore;
  int id=0;
  char ctline[512];

  clamp_free_genome(g, erepobj, ddh);
  apply_genome_to_params(genome, erepobj, ddh);

  if(!initialgen){
    ostringstream ss;
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    string call,stemp,part1,part2,part3; 
    string skf_dir=scratch+"/slakos.tmp_"+itoa(id,10);

    call="rm -rf "+skf_dir+"; mkdir -p "+skf_dir+"; ";

    for(i=0;i<erepobj.velem.size();i++){
      part1=erepobj.velem[i].name;  
      ss.str(std::string());
      ss.precision(erepobj.velem[i].radius[0].precision);
      ss<<std::fixed<<erepobj.velem[i].radius[0].r;
      part1+="_d"+ss.str();
      for(k=1;k<erepobj.velem[i].lmax+2.;k++){
        ss.str(std::string());
        ss.precision(erepobj.velem[i].radius[k].precision);
        ss<<std::fixed<<erepobj.velem[i].radius[k].r;
        part1+="_w"+ss.str();  
      }
      for(j=0;j<erepobj.velem.size();j++){
        part2=erepobj.velem[j].name;  
        ss.str(std::string());
        ss.precision(erepobj.velem[j].radius[0].precision);
        ss<<std::fixed<<erepobj.velem[j].radius[0].r;
        part2+="_d"+ss.str();  
        for(k=1;k<erepobj.velem[j].lmax+2.;k++){
          ss.str(std::string());
          ss.precision(erepobj.velem[j].radius[k].precision);
          ss<<std::fixed<<erepobj.velem[j].radius[k].r;
          part2+="_w"+ss.str();  
        }
        part3=pair_onsite_suffix(ddh, erepobj.velem[i].name, erepobj.velem[j].name);
        if(erepobj.lc){
          ss.str(std::string());
          ss.precision(2);ss<<std::fixed<<erepobj.omega;
          stemp="omega_"+ss.str()+"_"+part1+"-"+part2+part3+".skf";
        }else{
          stemp=part1+"-"+part2+part3+".skf";
        }
        if(erepobj.addskf){
          string addpair = erepobj.libadddir+"/"+erepobj.velem[i].name+"-"+erepobj.velem[j].name+".skf";
          if(!access(addpair.c_str(), F_OK)){
            continue;
          }
        }
        erepobj.checkskf(stemp,1);
        call+=" cp "+erepobj.libdir+"/"+stemp+" "+ skf_dir+"/"+erepobj.velem[i].name+"-"+erepobj.velem[j].name+".skf; ";
      }
    }
    if(erepobj.addskf) call+=" cp -f "+erepobj.libadddir+"/*.skf " + skf_dir+"/; ";

    i=system(call.c_str());	

    if (erepobj.fit_type==0){
      fstream tinfile, toutfile;
      tinfile.open(erepobj.frep_in.c_str(),ios::in);
      if(!tinfile.is_open()){
        cout<<"unable to open "<<erepobj.frep_in<<" file.\n";
        exit(1);
      }
      call = skf_dir+"/rep.in";
      toutfile.open(call.c_str(), ios::out);
      if(!toutfile.is_open()){
        cout<<"unable to open "<<call<<" file.\n";
        exit(1);
      }
      toutfile<<std::fixed;
      if( ddh.td3 || ddh.tdamph || ddh.thubbardderivs || ddh.tdamphver2){
        if(ddh.td3){
          toutfile<<"\n$d3:\n"; 
          for(i=0;i<ddh.d3.size();i++){ 
            toutfile<<ddh.d3[i].name<<" "<<ddh.d3[i].value<<" "<<ddh.d3[i].value<<" "<<ddh.d3[i].value<<" "<<ddh.d3[i].precision<<endl;
          }
          toutfile<<"$end\n"; 
        }
        if(ddh.tdamph){
          toutfile<<"\n$damph:\n"; 
          toutfile<<ddh.damph.value<<" "<<ddh.damph.value<<" "<<ddh.damph.value<<" "<<ddh.damph.precision<<endl;
          toutfile<<"$end\n"; 
        }
        if(ddh.thubbardderivs){
          if(ddh.thirdorderfull) toutfile<<"\n$hubbardderivsfull:\n"; 
          else toutfile<<"\n$hubbardderivs:\n"; 
          for(i=0;i<ddh.hubbardderivs.size();i++){ 
            toutfile<<ddh.hubbardderivs[i].name<<" "<<ddh.hubbardderivs[i].value<<" "<<ddh.hubbardderivs[i].value<<" "<<ddh.hubbardderivs[i].value<<" "<<ddh.hubbardderivs[i].precision<<endl;
          }
          toutfile<<"$end\n"; 
        }
        if(ddh.tdamphver2){
          toutfile<<"\n$damphver2:\n"; 
          for(i=0;i<ddh.damphver2.size();i++){ 
            toutfile<<ddh.damphver2[i].name<<" "<<ddh.damphver2[i].value<<" "<<ddh.damphver2[i].value<<" "<<ddh.damphver2[i].value<<" "<<ddh.damphver2[i].precision<<endl;
          }
          toutfile<<"$end\n"; 
        }
      }
      while(tinfile.getline(ctline,512)){
        toutfile<<ctline<<endl;
      }
      tinfile.close();
      toutfile.close();
    }

    tscore=erepobj.score(id); 

    return tscore;
  }else{
    return 0.000001;
  }
}

  //for(i=0;i<erepobj.velem.size();i++){
  //  part1=erepobj.velem[i].name;  
  //  ss.precision(erepobj.velem[i].radius[0].precision); ss<<std::fixed<<erepobj.velem[i].radius[0].r;
  //  part1+="_d"+ss.str(); 
  //  for(k=1;k<erepobj.velem[i].lmax+2.;k++){
  //    ss.precision(erepobj.velem[i].radius[k].precision); ss<<std::fixed<<erepobj.velem[i].radius[k].r;
  //    part1+="_w"+ss.str();  
  //  }
  //  for(j=0;j<erepobj.velem.size();j++){
  //    part2=erepobj.velem[j].name;  
  //    ss.precision(erepobj.velem[j].radius[0].precision); ss<<std::fixed<<erepobj.velem[j].radius[0].r;
  //    part2+="_d"+ss.str();  
  //    for(k=1;k<erepobj.velem[j].lmax+2.;k++){
  //      ss.precision(erepobj.velem[j].radius[k].precision); ss<<std::fixed<<erepobj.velem[j].radius[k].r;
  //      part2+="_w"+ss.str();  
  //    }
  //    stemp=erepobj.libdir+"/"+part1+"-"+part2+".skf";
  //    erepobj.checkskf(stemp);
  //    call+=" cp "+stemp+" "+ skf_dir+"; ";
  //  }
  //}

void MyInitializer(GAGenome& g) {
  GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
  int i,j,k,idx;
  double tmp;

  if(readr==true){
    idx=0;
    for(i=0;i<erepobj.velem.size();i++){
      for(k=0;k<erepobj.velem[i].lmax+2.;k++){
        if(gene_is_free(erepobj.velem[i].radius[k].minr, erepobj.velem[i].radius[k].maxr)){
          genome.gene(idx,erepobj.velem[i].radius[k].r);
          idx++;
        }
      }
    }
    if(ddh.td3){
      for(i=0;i<ddh.d3.size();i++){ 
        if(gene_is_free(ddh.d3[i].min, ddh.d3[i].max)){
          genome.gene(idx,ddh.d3[i].value);
          idx++;
        }
      }
    }
    if(ddh.tdamph){
      if(gene_is_free(ddh.damph.min, ddh.damph.max)){
        genome.gene(idx,ddh.damph.value);
        idx++;
      }
    }
    if(ddh.thubbardderivs){
      for(i=0;i<ddh.hubbardderivs.size();i++){ 
        if(gene_is_free(ddh.hubbardderivs[i].min, ddh.hubbardderivs[i].max)){
          genome.gene(idx,ddh.hubbardderivs[i].value);
          idx++;
        }
      }
    }
    if(ddh.tdamphver2){
      for(i=0;i<ddh.damphver2.size();i++){ 
        if(gene_is_free(ddh.damphver2[i].min, ddh.damphver2[i].max)){
          genome.gene(idx,ddh.damphver2[i].value);
          idx++;
        }
      }
    }
    if(ddh.thubbards){
      for(i=0;i<ddh.hubbards.size();i++){ 
        if(gene_is_free(ddh.hubbards[i].min, ddh.hubbards[i].max)){
          genome.gene(idx,ddh.hubbards[i].value);
          idx++;
        }
      }
    }
    if(ddh.tvorbes){
      for(i=0;i<ddh.vorbes.size();i++){ 
        if(gene_is_free(ddh.vorbes[i].min, ddh.vorbes[i].max)){
          genome.gene(idx,ddh.vorbes[i].value);
          idx++;
        }
      }
    }
    readr=false;
  }else{
    idx=0;
    for(i=0;i<erepobj.velem.size();i++){
      for(k=0;k<erepobj.velem[i].lmax+2.;k++){
        if(gene_is_free(erepobj.velem[i].radius[k].minr, erepobj.velem[i].radius[k].maxr)){
          genome.gene(idx,GARandomFloat(erepobj.velem[i].radius[k].minr,erepobj.velem[i].radius[k].maxr));
          idx++;
        }
      }
    }
    if(ddh.td3){
      for(i=0;i<ddh.d3.size();i++){ 
        if(gene_is_free(ddh.d3[i].min, ddh.d3[i].max)){
          genome.gene(idx,GARandomFloat(ddh.d3[i].min,ddh.d3[i].max));
          idx++;
        }
      }
    }
    if(ddh.tdamph){
      if(gene_is_free(ddh.damph.min, ddh.damph.max)){
        genome.gene(idx,GARandomFloat(ddh.damph.min,ddh.damph.max));
        idx++;
      }
    }
    if(ddh.thubbardderivs){
      for(i=0;i<ddh.hubbardderivs.size();i++){ 
        if(gene_is_free(ddh.hubbardderivs[i].min, ddh.hubbardderivs[i].max)){
          genome.gene(idx,GARandomFloat(ddh.hubbardderivs[i].min,ddh.hubbardderivs[i].max));
          idx++;
        }
      }
    }
    if(ddh.tdamphver2){
      for(i=0;i<ddh.damphver2.size();i++){ 
        if(gene_is_free(ddh.damphver2[i].min, ddh.damphver2[i].max)){
          genome.gene(idx,GARandomFloat(ddh.damphver2[i].min,ddh.damphver2[i].max));
          idx++;
        }
      }
    }
    if(ddh.thubbards){
      for(i=0;i<ddh.hubbards.size();i++){ 
        if(gene_is_free(ddh.hubbards[i].min, ddh.hubbards[i].max)){
          genome.gene(idx,GARandomFloat(ddh.hubbards[i].min,ddh.hubbards[i].max));
          idx++;
        }
      }
    }
    if(ddh.tvorbes){
      for(i=0;i<ddh.vorbes.size();i++){ 
        if(gene_is_free(ddh.vorbes[i].min, ddh.vorbes[i].max)){
          genome.gene(idx,GARandomFloat(ddh.vorbes[i].min,ddh.vorbes[i].max));
          idx++;
        }
      }
    }
  }
    
  if(restart==true){
    idx=0;
    for(i=0;i<erepobj.velem.size();i++){
      for(k=0;k<erepobj.velem[i].lmax+2.;k++){
        if(gene_is_free(erepobj.velem[i].radius[k].minr, erepobj.velem[i].radius[k].maxr)){
          tmp=999999999.9;
          restartfile>>tmp;
          if(tmp!=999999999.9){
            genome.gene(idx,tmp);
          }
          idx++;
        }
      }
    }
    if(ddh.td3){
      for(i=0;i<ddh.d3.size();i++){ 
        if(gene_is_free(ddh.d3[i].min, ddh.d3[i].max)){
          tmp=999999999.9;
          restartfile>>tmp;
          if(tmp!=999999999.9){
            genome.gene(idx,tmp);
          }
          idx++;
        }
      }
    }
    if(ddh.tdamph){
      if(gene_is_free(ddh.damph.min, ddh.damph.max)){
        tmp=999999999.9;
        restartfile>>tmp;
        if(tmp!=999999999.9){
          genome.gene(idx,tmp);
        }
        idx++;
      }
    }
    if(ddh.thubbardderivs){
      for(i=0;i<ddh.hubbardderivs.size();i++){ 
        if(gene_is_free(ddh.hubbardderivs[i].min, ddh.hubbardderivs[i].max)){
          tmp=999999999.9;
          restartfile>>tmp;
          if(tmp!=999999999.9){
            genome.gene(idx,tmp);
          }
          idx++;
        }
      }
    }
    if(ddh.tdamphver2){
      for(i=0;i<ddh.damphver2.size();i++){ 
        if(gene_is_free(ddh.damphver2[i].min, ddh.damphver2[i].max)){
          tmp=999999999.9;
          restartfile>>tmp;
          if(tmp!=999999999.9){
            genome.gene(idx,tmp);
          }
          idx++;
        }
      }
    }
    if(ddh.thubbards){
      for(i=0;i<ddh.hubbards.size();i++){ 
        if(gene_is_free(ddh.hubbards[i].min, ddh.hubbards[i].max)){
          tmp=999999999.9;
          restartfile>>tmp;
          if(tmp!=999999999.9){
            genome.gene(idx,tmp);
          }
          idx++;
        }
      }
    }
    if(ddh.tvorbes){
      for(i=0;i<ddh.vorbes.size();i++){ 
        if(gene_is_free(ddh.vorbes[i].min, ddh.vorbes[i].max)){
          tmp=999999999.9;
          restartfile>>tmp;
          if(tmp!=999999999.9){
            genome.gene(idx,tmp);
          }
          idx++;
        }
      }
    }
    restartfile>>tmp;
  }

  clamp_free_genome(g, erepobj, ddh);
}

int MyMutator(GAGenome& g, double pmut){
  GA1DArrayGenome<double>& child = (GA1DArrayGenome<double>&)g;
  int i,j,k,idx,iMut;
  double value;

  if(pmut <= 0.0) return(0);
  idx=0;
  iMut=0;
  for(i=0;i<erepobj.velem.size();i++){ 
    for(j=0; j<erepobj.velem[i].lmax+2; j++){
      if(!gene_is_free(erepobj.velem[i].radius[j].minr, erepobj.velem[i].radius[j].maxr)) continue;
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(erepobj.velem[i].radius[j].delta);
        value += child.gene(idx);
        value  = GAMax(erepobj.velem[i].radius[j].minr, value);
        value  = GAMin(erepobj.velem[i].radius[j].maxr, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  if(ddh.td3){
    for(i=0;i<ddh.d3.size();i++){ 
      if(!gene_is_free(ddh.d3[i].min, ddh.d3[i].max)) continue;
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(ddh.d3[i].delta);
        value += child.gene(idx);
        value  = GAMax(ddh.d3[i].min, value);
        value  = GAMin(ddh.d3[i].max, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  if(ddh.tdamph){
    if(gene_is_free(ddh.damph.min, ddh.damph.max)){
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(ddh.damph.delta);
        value += child.gene(idx);
        value  = GAMax(ddh.damph.min, value);
        value  = GAMin(ddh.damph.max, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  if(ddh.thubbardderivs){
    for(i=0;i<ddh.hubbardderivs.size();i++){ 
      if(!gene_is_free(ddh.hubbardderivs[i].min, ddh.hubbardderivs[i].max)) continue;
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(ddh.hubbardderivs[i].delta);
        value += child.gene(idx);
        value  = GAMax(ddh.hubbardderivs[i].min, value);
        value  = GAMin(ddh.hubbardderivs[i].max, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  if(ddh.tdamphver2){
    for(i=0;i<ddh.damphver2.size();i++){ 
      if(!gene_is_free(ddh.damphver2[i].min, ddh.damphver2[i].max)) continue;
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(ddh.damphver2[i].delta);
        value += child.gene(idx);
        value  = GAMax(ddh.damphver2[i].min, value);
        value  = GAMin(ddh.damphver2[i].max, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  if(ddh.thubbards){
    for(i=0;i<ddh.hubbards.size();i++){ 
      if(!gene_is_free(ddh.hubbards[i].min, ddh.hubbards[i].max)) continue;
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(ddh.hubbards[i].delta);
        value += child.gene(idx);
        value  = GAMax(ddh.hubbards[i].min, value);
        value  = GAMin(ddh.hubbards[i].max, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  if(ddh.tvorbes){
    for(i=0;i<ddh.vorbes.size();i++){ 
      if(!gene_is_free(ddh.vorbes[i].min, ddh.vorbes[i].max)) continue;
      if(GAFlipCoin(pmut)){
        value  = GAGaussianFloat(ddh.vorbes[i].delta);
        value += child.gene(idx);
        value  = GAMax(ddh.vorbes[i].min, value);
        value  = GAMin(ddh.vorbes[i].max, value);
        child.gene(idx,value);
        iMut++;
      }
      idx++;
    }
  }
  return(iMut);
}


double MyComparator(const GAGenome& g1, const GAGenome& g2) {
  GA1DArrayGenome<double>& gnome1 = (GA1DArrayGenome<double>&)g1;
  GA1DArrayGenome<double>& gnome2 = (GA1DArrayGenome<double>&)g2;
  double diff=0.0; 
  int i;
  for(i=0;i<gnome1.length();i++){ 
    diff+=abs(gnome2.gene(i)-gnome1.gene(i));
  }
//return diff;
  return 1; 
}


// Pick a crossover breakpoint in [0, length], optionally snapped to multiples of crossover_align.
static unsigned int aligned_crossover_site(unsigned int length){
  if(length == 0) return 0;
  if(crossover_align <= 0){
    return (unsigned int)GARandomInt(0, (int)length);
  }
  unsigned int stride = (unsigned int)crossover_align;
  unsigned int nblock = length / stride;
  if(nblock == 0){
    return (unsigned int)GARandomInt(0, (int)length);
  }
  return stride * (unsigned int)GARandomInt(0, (int)nblock);
}

int MyOnePointCrossover(const GAGenome& p1, const GAGenome& p2,
  GAGenome* c1, GAGenome* c2){
  const GA1DArrayGenome<double> &mom=DYN_CAST(const GA1DArrayGenome<double> &, p1);
  const GA1DArrayGenome<double> &dad=DYN_CAST(const GA1DArrayGenome<double> &, p2);

  int nc=0;
  unsigned int momsite, momlen;
  unsigned int dadsite, dadlen;

  if(c1 && c2){
    GA1DArrayGenome<double> &sis=DYN_CAST(GA1DArrayGenome<double> &, *c1);
    GA1DArrayGenome<double> &bro=DYN_CAST(GA1DArrayGenome<double> &, *c2);

    if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE &&
    bro.resizeBehaviour() == GAGenome::FIXED_SIZE){
      if(mom.length() != dad.length() ||
      sis.length() != bro.length() ||
      sis.length() != mom.length()){
        GAErr(GA_LOC, mom.className(), "one-point cross", gaErrSameLengthReqd);
        return nc;
      }
      momsite = dadsite = aligned_crossover_site(mom.length());
      momlen = dadlen = mom.length() - momsite;
    }else if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE ||
    bro.resizeBehaviour() == GAGenome::FIXED_SIZE){
      GAErr(GA_LOC, mom.className(), "one-point cross", gaErrSameBehavReqd);
      return nc;
    }else{
      momsite = aligned_crossover_site(mom.length());
      dadsite = aligned_crossover_site(dad.length());
      momlen = mom.length() - momsite;
      dadlen = dad.length() - dadsite;
      sis.resize(momsite+dadlen);
      bro.resize(dadsite+momlen);
    }

    sis.copy(mom, 0, 0, momsite);
    sis.copy(dad, momsite, dadsite, dadlen);
    bro.copy(dad, 0, 0, dadsite);
    bro.copy(mom, dadsite, momsite, momlen);
    nc = 2;
  }else if(c1 || c2){
    GA1DArrayGenome<double> &sis = (c1 ?
    DYN_CAST(GA1DArrayGenome<double> &, *c1) :
    DYN_CAST(GA1DArrayGenome<double> &, *c2));

    if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE){
      if(mom.length() != dad.length() || sis.length() != mom.length()){
        GAErr(GA_LOC, mom.className(), "one-point cross", gaErrSameLengthReqd);
        return nc;
      }
      momsite = dadsite = aligned_crossover_site(mom.length());
      momlen = dadlen = mom.length() - momsite;
    }else{
      momsite = aligned_crossover_site(mom.length());
      dadsite = aligned_crossover_site(dad.length());
      momlen = mom.length() - momsite;
      dadlen = dad.length() - dadsite;
      sis.resize(momsite+dadlen);
    }

    if(GARandomBit()){
      sis.copy(mom, 0, 0, momsite);
      sis.copy(dad, momsite, dadsite, dadlen);
    }else{
      sis.copy(dad, 0, 0, dadsite);
      sis.copy(mom, dadsite, momsite, momlen);
    }
    nc = 1;
  }
  return nc;
}

int MyTwoPointCrossover(const GAGenome& p1, const GAGenome& p2,
  GAGenome* c1, GAGenome* c2){
  const GA1DArrayGenome<double> &mom=DYN_CAST(const GA1DArrayGenome<double> &, p1);
  const GA1DArrayGenome<double> &dad=DYN_CAST(const GA1DArrayGenome<double> &, p2);

  int nc=0;
  unsigned int momsite[2], momlen[2];
  unsigned int dadsite[2], dadlen[2];

  if(c1 && c2){
    GA1DArrayGenome<double> &sis=DYN_CAST(GA1DArrayGenome<double> &, *c1);
    GA1DArrayGenome<double> &bro=DYN_CAST(GA1DArrayGenome<double> &, *c2);

    if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE &&
    bro.resizeBehaviour() == GAGenome::FIXED_SIZE){
      if(mom.length() != dad.length() ||
      sis.length() != bro.length() ||
      sis.length() != mom.length()){
        GAErr(GA_LOC, mom.className(), "two-point cross", gaErrSameLengthReqd);
        return nc;
      }
      momsite[0] = aligned_crossover_site(mom.length());
      momsite[1] = aligned_crossover_site(mom.length());
      if(momsite[0] > momsite[1]) SWAP(momsite[0], momsite[1]);
      momlen[0] = momsite[1] - momsite[0];
      momlen[1] = mom.length() - momsite[1];

      dadsite[0] = momsite[0];
      dadsite[1] = momsite[1];
      dadlen[0] = momlen[0];
      dadlen[1] = momlen[1];
    }else if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE ||
    bro.resizeBehaviour() == GAGenome::FIXED_SIZE){
      return nc;
    }else{
      momsite[0] = aligned_crossover_site(mom.length());
      momsite[1] = aligned_crossover_site(mom.length());
      if(momsite[0] > momsite[1]) SWAP(momsite[0], momsite[1]);
      momlen[0] = momsite[1] - momsite[0];
      momlen[1] = mom.length() - momsite[1];

      dadsite[0] = aligned_crossover_site(dad.length());
      dadsite[1] = aligned_crossover_site(dad.length());
      if(dadsite[0] > dadsite[1]) SWAP(dadsite[0], dadsite[1]);
      dadlen[0] = dadsite[1] - dadsite[0];
      dadlen[1] = dad.length() - dadsite[1];

      sis.resize(momsite[0]+dadlen[0]+momlen[1]);
      bro.resize(dadsite[0]+momlen[0]+dadlen[1]);
    }

    sis.copy(mom, 0, 0, momsite[0]);
    sis.copy(dad, momsite[0], dadsite[0], dadlen[0]);
    sis.copy(mom, momsite[0]+dadlen[0], momsite[1], momlen[1]);
    bro.copy(dad, 0, 0, dadsite[0]);
    bro.copy(mom, dadsite[0], momsite[0], momlen[0]);
    bro.copy(dad, dadsite[0]+momlen[0], dadsite[1], dadlen[1]);

    nc = 2;
   
  }else if(c1 || c2){
    GA1DArrayGenome<double> &sis = (c1 ?
    DYN_CAST(GA1DArrayGenome<double> &, *c1) :
    DYN_CAST(GA1DArrayGenome<double> &, *c2));

    if(sis.resizeBehaviour() == GAGenome::FIXED_SIZE){
      if(mom.length() != dad.length() || sis.length() != mom.length()){
        GAErr(GA_LOC, mom.className(), "two-point cross", gaErrSameLengthReqd);
        return nc;
      }
      momsite[0] = aligned_crossover_site(mom.length());
      momsite[1] = aligned_crossover_site(mom.length());
      if(momsite[0] > momsite[1]) SWAP(momsite[0], momsite[1]);
      momlen[0] = momsite[1] - momsite[0];
      momlen[1] = mom.length() - momsite[1];

      dadsite[0] = momsite[0];
      dadsite[1] = momsite[1];
      dadlen[0] = momlen[0];
      dadlen[1] = momlen[1];
    }else{
      momsite[0] = aligned_crossover_site(mom.length());
      momsite[1] = aligned_crossover_site(mom.length());
      if(momsite[0] > momsite[1]) SWAP(momsite[0], momsite[1]);
      momlen[0] = momsite[1] - momsite[0];
      momlen[1] = mom.length() - momsite[1];

      dadsite[0] = aligned_crossover_site(dad.length());
      dadsite[1] = aligned_crossover_site(dad.length());
      if(dadsite[0] > dadsite[1]) SWAP(dadsite[0], dadsite[1]);
      dadlen[0] = dadsite[1] - dadsite[0];
      dadlen[1] = dad.length() - dadsite[1];

      sis.resize(momsite[0]+dadlen[0]+momlen[1]);
    }

    if(GARandomBit()){
      sis.copy(mom, 0, 0, momsite[0]);
      sis.copy(dad, momsite[0], dadsite[0], dadlen[0]);
      sis.copy(mom, momsite[0]+dadlen[0], momsite[1], momlen[1]);
    }else{
      sis.copy(dad, 0, 0, dadsite[0]);
      sis.copy(mom, dadsite[0], momsite[0], momlen[0]);
      sis.copy(dad, dadsite[0]+momlen[0], dadsite[1], dadlen[1]);
    }

    nc = 1;
  
  }

  return nc;
}

int MyUniformCrossover(const GAGenome& p1, const GAGenome& p2,
  GAGenome* c1, GAGenome* c2){
  const GA1DArrayGenome<double> &mom=DYN_CAST(const GA1DArrayGenome<double> &, p1);
  const GA1DArrayGenome<double> &dad=DYN_CAST(const GA1DArrayGenome<double> &, p2);

  int nc=0;
  unsigned int i;

  if(c1 && c2){
    GA1DArrayGenome<double> &sis=DYN_CAST(GA1DArrayGenome<double> &, *c1);
    GA1DArrayGenome<double> &bro=DYN_CAST(GA1DArrayGenome<double> &, *c2);

    if(sis.length() != mom.length() || bro.length() != dad.length() ||
       mom.length() != dad.length()){
      GAErr(GA_LOC, mom.className(), "uniform cross", gaErrSameLengthReqd);
      return nc;
    }
    for(i=0; i<mom.length(); i++){
      if(GAFlipCoin(0.5)){
        sis.gene(i, mom.gene(i));
        bro.gene(i, dad.gene(i));
      }else{
        sis.gene(i, dad.gene(i));
        bro.gene(i, mom.gene(i));
      }
    }
    nc = 2;
  }else if(c1 || c2){
    GA1DArrayGenome<double> &sis = (c1 ?
    DYN_CAST(GA1DArrayGenome<double> &, *c1) :
    DYN_CAST(GA1DArrayGenome<double> &, *c2));

    if(sis.length() != mom.length() || mom.length() != dad.length()){
      GAErr(GA_LOC, mom.className(), "uniform cross", gaErrSameLengthReqd);
      return nc;
    }
    for(i=0; i<mom.length(); i++){
      if(GAFlipCoin(0.5)){
        sis.gene(i, mom.gene(i));
      }else{
        sis.gene(i, dad.gene(i));
      }
    }
    nc = 1;
  }
  return nc;
}

