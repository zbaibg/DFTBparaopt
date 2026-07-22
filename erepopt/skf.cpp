#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <vector>
#include <cstring>
#include <cmath>
#include <omp.h>
#include <unistd.h>
#include <ga/ga.h>
#include <ga/std_stream.h>
#include "tools.hpp"
#include "ddh.hpp"
#include "erepobj.hpp"
#include "genes.hpp"
#include "mpi.h"

extern int power;
extern sddh ddh;
extern string scratch;
extern int cpu_number;
using namespace std;
const double approximaterezo=0.000001; 

bool Erepobj::checkskf(const string ifilename, int type){
  ifstream ifile;
  ofstream ofile;
  string filename;
  bool passlc,passspline,passgrid,calculating,pass; 
  char cline[5120],ctemp1[512];
  int nline;
  double otmp;
  pass=passlc=passspline=passgrid=calculating=false;
  filename=libdir+"/"+ifilename;
  ifile.open(filename.c_str());
  if ( !ifile ) {
    ifile.close();
  }else{
    nline=0;
    string stemp;
    while(ifile.getline(cline,5120)){
      nline++;
      stemp=cline;
      if(stemp.find("Spline")!=string::npos){
        passspline=true; 
        if(nline>ngrid) passgrid=true;
      }else if(stemp.find("LC")!=string::npos){
        passlc=true;  
        sscanf(cline,"%s %le",ctemp1,&otmp);
        if(abs(otmp-omega)>approximaterezo) passlc=false;
      }
    }
    ifile.close();
  }

  if(lc){
    //if(passspline && passgrid && passlc) pass=true;
    if(passspline && passlc) pass=true;
  }else{
    //if(passspline && passgrid ) pass=true;
    if(passspline) pass=true;
  }

  if(type==1){
    if(!pass){
      cerr << endl << "ERROR: " << endl << endl;
      cerr << libdir<<"/"<<ifilename << ": is not correct" << endl << "exit score" << endl ;
      exit(1);
    }else return true;
  }else{
    if(!pass) return false;  
    else return true;
  }
} 

void Erepobj::writeskgen(const string& tmp_dir, const GAGenome& g) {
  apply_genome_to_params(g, *this, ddh);
  int i,j,k,nvars;
  string sorbital;
  string stemp = tmp_dir+"/"+"skdef.hsd";
  ofstream fout(stemp.c_str());
  bool pass; 

  fout.precision(6);

  fout<<"SkdefVersion = "<<skdefversion<<"\n\n";

  fout<<"Globals {\n";
  if(lc){
    fout<<"  func = xchybrid {omega = "<<omega<<"}\n";
  }else{
    fout<<"  XCFunctional ="<<xcfunctional<<"\n";
  }
  fout<<"  Superposition = "<<superposition<<"\n";
  fout<<"}"<<"\n\n";

  fout<<"AtomParameters {\n\n";
  fout<<"\
  $OccupationsNe {\n\
    1S = 1.0 1.0\n\
    2S = 1.0 1.0\n\
    2P = 3.0 3.0\n\
  }\n\
  \n\
  $OccupationsAr {\n\
    $OccupationsNe\n\
    3S = 1.0 1.0\n\
    3P = 3.0 3.0\n\
  }\n\
  \n\
  $OccupationsKr {\n\
    $OccupationsAr\n\
    3D = 5.0 5.0\n\
    4S = 1.0 1.0\n\
    4P = 3.0 3.0\n\
  }\n\
  \n\
  $OccupationsXe {\n\
    $OccupationsKr\n\
    4D = 5.0 5.0\n\
    5S = 1.0 1.0\n\
    5P = 3.0 3.0\n\
  }\n\
  \n\
  $OccupationsHg {\n\
    $OccupationsXe\n\
    4F = 7.0 7.0\n\
    5D = 5.0 5.0\n\
    6S = 1.0 1.0\n\
  }\n\
  \n\
  $OccupationsRn {\n\
    $OccupationsHg\n\
    6P = 3.0 3.0\n\
  }\n\n";
 
  for(i=0; i<velem.size();i++){
    pass=false;
    for(j=0; j<atomparameters.size();j++){
      if(velem[i].name==atomparameters[j].name){pass=true; break;}
    }
    if(!pass){
      cerr <<"\nERROR:\nCan not find atomparameters for " << velem[i].name << endl ;
      exit(1);
    }
    fout<<"  "<<atomparameters[j].name<<" {\n";
    for(k=0; k<atomparameters[j].atomconfig.size();k++){
      fout<<atomparameters[j].atomconfig[k]<<"\n";
    }
    fout<<"    DftbAtom {\n";
    fout<<"      ShellResolved = "<<atomparameters[j].shellresolved<<"\n";
    fout.precision(velem[i].radius[0].precision);
    fout<<"      DensityCompression = PowerCompression { Power = "<<atomparameters[j].power<<"; Radius = "<<std::fixed<<velem[i].radius[0].r<<"}\n";
    fout<<"      WaveCompressions = SingleAtomCompressions {\n";
    nvars=velem[i].lmax+2;
    for(k=1; k<nvars; k++){
      fout.precision(velem[i].radius[k].precision);
      if(k==1) sorbital="S"; else if(k==2) sorbital="P"; else if(k==3) sorbital="D";
      fout<<"        "<<sorbital<<" = PowerCompression { Power = "<<atomparameters[j].power<<"; Radius = "<<std::fixed<<velem[i].radius[k].r<<"}\n";
    }
    fout<<"      }\n";
    if(ddh.thubbards){
      fout<<"      CustomizedHubbards {\n";
      for(k=0;k<ddh.hubbards.size();k++){ 
        if(velem[i].name==ddh.hubbards[k].name){
          fout.precision(ddh.hubbards[k].precision);
          fout<<"        "<<ddh.hubbards[k].type<<" = "<< ddh.hubbards[k].value<<endl;
        }
      }
      fout<<"      }\n";
    }
    if(ddh.tvorbes){
      fout<<"      CustomizedOnsites {\n";
      for(k=0;k<ddh.vorbes.size();k++){ 
        if(velem[i].name==ddh.vorbes[k].name){
          fout.precision(ddh.vorbes[k].precision);
          fout<<"        "<<ddh.vorbes[k].type<<" = "<< ddh.vorbes[k].value<<endl;
        }
      }
      fout<<"      }\n";
    }
    fout<<"    }\n";
    fout<<"  }\n";
  } 
  fout<<"}\n\n";
       
  fout<<"OnecenterParameters {\n\n";
  for(k=0;k<onecenterparameters.size();k++){
    fout<<onecenterparameters[k]<<"\n";
  }
  fout<<"}\n\n";

  fout<<"TwoCenterParameters {\n\n";
  for(k=0;k<twocenterparameters.size();k++){
    fout<<twocenterparameters[k]<<"\n";
  }
  fout<<"}\n\n";
} 

void Erepobj::makeskf(const GAGenome& g) {
  apply_genome_to_params(g, *this, ddh);
  int ie1,ie2,i,j,k,it;
  bool pass,pass2;
  ostringstream ss;
  int id=0;
  string sfilenametmp;

  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  string call,stemp,stemp2,part1,part2,part3; 
  string skf_dir=scratch+"/buildsks.tmp_"+itoa(id,10);
  call="rm -rf "+skf_dir+"; mkdir -p "+skf_dir+";";
	it=system(call.c_str());	
  writeskgen(skf_dir,g);

  ie1=0;
  for(i=0;i<velem.size();i++){
    part1=velem[i].name;  
    ss.str(std::string());
    ss.precision(velem[i].radius[0].precision); ss<<std::fixed<<velem[i].radius[0].r;
    part1+="_d"+ss.str();
    for(k=1;k<velem[i].lmax+2.;k++){
      ss.str(std::string());
      ss.precision(velem[i].radius[k].precision); ss<<std::fixed<<velem[i].radius[k].r;
      part1+="_w"+ss.str();  
    }
    for(j=0;j<velem.size();j++){
      part2=velem[j].name;  
      ss.str(std::string());
      ss.precision(velem[j].radius[0].precision); ss<<std::fixed<<velem[j].radius[0].r;
      part2+="_d"+ss.str();  
      for(k=1;k<velem[j].lmax+2.;k++){
        ss.str(std::string());
        ss.precision(velem[j].radius[k].precision); ss<<std::fixed<<velem[j].radius[k].r;
        part2+="_w"+ss.str();  
      }
      part3=pair_onsite_suffix(ddh, velem[i].name, velem[j].name);
      if(lc){
        ss.str(std::string());
        ss.precision(2);ss<<std::fixed<<omega;
        stemp="omega_"+ss.str()+"_"+part1+"-"+part2+part3+".skf";
        stemp2="omega_"+ss.str()+"_"+part2+"-"+part1+part3+".skf";
      }else{
        stemp=part1+"-"+part2+part3+".skf";
        stemp2=part2+"-"+part1+part3+".skf";
      }
      if( i!=j){
//      #pragma omp critical (checkskf02)
//      if(skfexist[i][j]){
        if(addskf) sfilenametmp=libadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
        if(addskf && !access(sfilenametmp.c_str(), F_OK )){
          pass=true;
          pass2=true;
        }else{
          pass=checkskf(stemp,0);
          pass2=checkskf(stemp2,0);
        }
        if((!pass)&&(!pass2)){
          call=" cd "+skf_dir+"; export OMP_NUM_THREADS="+itoa(cpu_number,10)+";";
          call+= skgen +" -l error -o "+onecent+" -t "+twocent+" sktable "+velem[i].name+" "+velem[j].name+" >log;";
          if(lc){
//          if(repexist[i][j]){
            if(addrep) sfilenametmp=repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
            if(addrep && !access(sfilenametmp.c_str(), F_OK )){
              call+="omega=`grep \"LC\" "+velem[i].name+"-"+velem[j].name+".skf`;\
                  sed -i '/RangeSep/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                  sed -i '/LC/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                  cat "+ repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf >> "+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"RangeSep\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"$omega\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  sed -i '/RangeSep/d' "+velem[j].name+"-"+velem[i].name+".skf;\
                  sed -i '/LC/d' "+velem[j].name+"-"+velem[i].name+".skf;\
                  cat "+ repadddir+"/"+velem[j].name+"-"+velem[i].name+".skf >> "+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"RangeSep\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"$omega\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  cd ../;\
                  mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";\
                  mv -n "+skf_dir+"/"+velem[j].name+"-"+velem[i].name+".skf "+ libdir+"/"+stemp2+";";
            }else{
              call+="omega=`grep \"LC\" "+velem[i].name+"-"+velem[j].name+".skf`;\
                  sed -i '/RangeSep/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                  sed -i '/LC/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"Spline\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"1 9.99\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"0.00 9.99 0.0 0.0 0.0 0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"RangeSep\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"$omega\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  sed -i '/RangeSep/d' "+velem[j].name+"-"+velem[i].name+".skf;\
                  sed -i '/LC/d' "+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"Spline\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"1 9.99\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"0.0 0.0 0.0\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"0.00 9.99 0.0 0.0 0.0 0.0 0.0 0.0\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"RangeSep\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"$omega\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  cd ../;\
                  mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";\
                  mv -n "+skf_dir+"/"+velem[j].name+"-"+velem[i].name+".skf "+ libdir+"/"+stemp2+";";
            }
          }else{
//          if(repexist[i][j]){
            if(addrep) sfilenametmp=repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
            if(addrep && !access(sfilenametmp.c_str(), F_OK )){
              call+="cat "+ repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf >> "+velem[i].name+"-"+velem[j].name+".skf;\
                  cat "+ repadddir+"/"+velem[j].name+"-"+velem[i].name+".skf >> "+velem[j].name+"-"+velem[i].name+".skf;\
                  cd ../;\
                  mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";\
                  mv -n "+skf_dir+"/"+velem[j].name+"-"+velem[i].name+".skf "+ libdir+"/"+stemp2+";";
            }else{
              call+="echo \"Spline\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"1 9.99\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"0.00 9.99 0.0 0.0 0.0 0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                  echo \"Spline\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"1 9.99\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"0.0 0.0 0.0\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  echo \"0.00 9.99 0.0 0.0 0.0 0.0 0.0 0.0\" >>"+velem[j].name+"-"+velem[i].name+".skf;\
                  cd ../;\
                  mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";\
                  mv -n "+skf_dir+"/"+velem[j].name+"-"+velem[i].name+".skf "+ libdir+"/"+stemp2+";";
            }
          }
          it=system(call.c_str());	
        }
      }else{
 //     #pragma omp critical (checkskf01)
//      if(skfexist[i][j]){
        if(addskf) sfilenametmp=libadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
        if(addskf && !access(sfilenametmp.c_str(), F_OK )){
          pass=true;
        }else{
          pass=checkskf(stemp,0);
        }
        if(!pass){
          call=" cd "+skf_dir+"; export OMP_NUM_THREADS="+itoa(cpu_number,10)+";";
          call+= skgen +" -l error -o "+onecent+" -t "+twocent+" sktable "+velem[i].name+" "+velem[j].name+" >log;";
          if(lc){
//          if(repexist[i][j]){
            if(addrep) sfilenametmp=repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
            if(addrep && !access(sfilenametmp.c_str(), F_OK )){
              call+="omega=`grep \"LC\" "+velem[i].name+"-"+velem[j].name+".skf`;\
                   sed -i '/RangeSep/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                   sed -i '/LC/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                   cat "+ repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf >> "+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"RangeSep\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"$omega\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   cd ../;\
                   mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";";
            }else{
              call+="omega=`grep \"LC\" "+velem[i].name+"-"+velem[j].name+".skf`;\
                   sed -i '/RangeSep/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                   sed -i '/LC/d' "+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"Spline\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"1 9.99\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"0.00 9.99 0.0 0.0 0.0 0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"RangeSep\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"$omega\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   cd ../;\
                   mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";";
            }
          }else{
//          if(repexist[i][j]){
            if(addrep) sfilenametmp=repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
            if(addrep && !access(sfilenametmp.c_str(), F_OK )){
              call+="cat "+ repadddir+"/"+velem[i].name+"-"+velem[j].name+".skf >> "+velem[i].name+"-"+velem[j].name+".skf;\
                   cd ../;\
                   mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";";
            }else{
              call+="echo \"Spline\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"1 9.99\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   echo \"0.00 9.99 0.0 0.0 0.0 0.0 0.0 0.0\" >>"+velem[i].name+"-"+velem[j].name+".skf;\
                   cd ../;\
                   mv -n "+skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf "+ libdir+"/"+stemp+";";
            }
          }
          //cout<<"####################\n"<<call<<endl;
          it=system(call.c_str());	
        }
      }
    }
  }

  call="rm -rf "+skf_dir+";";
	it=system(call.c_str());	

} 

