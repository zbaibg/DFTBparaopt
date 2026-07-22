#include <iostream>
#include <string>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include <unistd.h>
#include "ga/ga.h"
#include <ga/std_stream.h>
#include "erepobj.hpp"
#include "tools.hpp"
#include "ddh.hpp"
#include "genes.hpp"

extern sddh ddh;
extern string scratch;
extern int cpu_number;

using namespace std;
extern bool runtest;

Erepobj::Erepobj(){
}

Erepobj::Erepobj(const string inputfile) {
  readinp(inputfile); 
}

void Erepobj::prepare(const string inputfile) {
  readinp(inputfile); 
}

double Erepobj::score(int id){
  int i,j,iq,ia;
  double a,score,tscore,tmscore;
  string stemp,skf_dir,call;
  char cline[512];
  skf_dir=scratch+"/slakos.tmp_"+itoa(id,10);


  if (fit_type==1){
    int index,length;
    double  * ev;
    length=0;
    for(i=0;i<vmol.size();i++){
      length=length+ vmol[i].nelectron/2;
    }
    ev = new double [length];
    index=0;
    for(i=0;i<vmol.size();i++){
      vmol[i].getev(id,index,ev,dftbversion);
    }
    a=0.1;
    score=999999999.9;
    for(ia=0;ia<9000;ia++){
      index=0;
      tscore=0.0;
      a=a+0.0001;
      a=1.0;
      for(i=0;i<vmol.size();i++){
        tmscore=0.0;
        for(j=0;j<vmol[i].nelectron/2;j++){
          if(score_type==0) tmscore+=vmol[i].weight*(ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0);
          else if(score_type==1) tmscore+=vmol[i].weight*abs(ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0);
          else if(score_type==2) tmscore+=vmol[i].weight*pow((ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0),2.0);
          else if(score_type==4) tmscore+=vmol[i].weight*pow((ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0),4.0);
          index++;
        }
        tscore=tscore+abs(tmscore);
      }
      if(tscore < score) score=tscore;
    }
    delete  [] ev;
  }else if (fit_type==2){
    score=0.0;
    for(i=0;i<vmol.size();i++){
      score+=vmol[i].score(id,score_type,dftbversion,false);
    } 
  }else{
    score=999999999.9;
    call  = "cp -rf "+fgrid+" "+skf_dir+";";
//  call += "cp "+frep_in+" "+skf_dir+"/rep.in;";
//  call += "cd "+skf_dir+"/; "+ repopt +" "+ frep_in +">rep.out; " + "cd ../ ";
    call += "cd "+skf_dir+"/; export OMP_NUM_THREADS="+itoa(cpu_number,10)+";"+ repopt +" rep.in >rep.out; " + "cd ../ ";

    iq = system(call.c_str());	

    stemp=skf_dir+"/rep.out";
    ifstream fin(stemp.c_str());
    if ( !fin ) {
      cerr << endl << "ERROR: " << endl << endl;
      cerr << "The system call: \"" << call << "\" did not produce the file: \""<<skf_dir<<"/rep.out\""
           << endl << "calculation of repopt was failed."
           << endl << "exit erepobj" << endl << endl;
      exit(1);
    }
    while(fin.getline(cline,512)){
      stemp=cline;
      if(stemp.find("final_score:")!=string::npos){
        fin >> score;
        break;
      }
    }
    fin.close();
  }

  return score;
}

void Erepobj::get_residual(const GAGenome& g) {
  restotU=0; restotS=0; restot2=0; restot4=0;
  GA1DArrayGenome<double>& genome = (GA1DArrayGenome<double>&)g;
  apply_genome_to_params(genome, *this, ddh);
  int ie1,ie2,i,j,k,ia;
  double tscore,tmscore;
  int id=0;

    ostringstream ss;
    id=0;
    string call,stemp,part1,part2,part3; 
    string skf_dir=scratch+"/slakos.tmp_"+itoa(id,10);
    call=" mkdir -p "+skf_dir+"; ";

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
        }else{
          stemp=part1+"-"+part2+part3+".skf";
        }
        if (fit_type!=0){
          // If an additional skf for this pair exists in libadddir, skip
          // checking the parameterized skf in libdir to keep behavior
          // consistent with makeskf. The actual copy from libadddir will
          // be done later by the bulk "cp libadddir/*.skf" command.
          if(addskf){
            string addpair = libadddir+"/"+velem[i].name+"-"+velem[j].name+".skf";
            if(!access(addpair.c_str(), F_OK)){
              continue;
            }
          }
          checkskf(stemp,1);
          call+=" cp "+libdir+"/"+stemp+" "+ skf_dir+"/"+velem[i].name+"-"+velem[j].name+".skf; ";
        }
      }
    }
    if(addskf) call+=" cp -f "+libadddir+"/*.skf " + skf_dir+"/; ";

    i=system(call.c_str());	
//}
  
  //tscore=score(id); 

  //tscore=0.0;
  //for(i=0;i<vmol.size();i++){
  //  tscore+=vmol[i].score(id,score_type,dftbversion,true);
  //}

  if (fit_type==1){
    int index,length;
    double  * ev,a,af,fscore;
    length=0;
    for(i=0;i<vmol.size();i++){
      length=length+ vmol[i].nelectron/2;
    }
    ev = new double [length];
    index=0;
    for(i=0;i<vmol.size();i++){
      vmol[i].getev(id,index,ev,dftbversion);
    }
    //index=0;
    //for(i=0;i<vmol.size();i++){
    //  for(j=0;j<vmol[i].nelectron/2;j++){
    //    cout<<ev[index]<<endl;
    //    index++;
    //  }
    //}
    af=a=0.1;
    fscore=999999999.9;
    for(ia=0;ia<9000;ia++){
      index=0;
      tscore=0.0;
      a=a+0.0001;
      a=1.0;
      for(i=0;i<vmol.size();i++){
        tmscore=0.0;
        for(j=0;j<vmol[i].nelectron/2;j++){
          if(score_type==0) tmscore+=vmol[i].weight*(ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0);
          else if(score_type==1) tmscore+=vmol[i].weight*abs(ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0);
          else if(score_type==2) tmscore+=vmol[i].weight*pow((ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0),2.0);
          else if(score_type==4) tmscore+=vmol[i].weight*pow((ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0),4.0);
          index++;
        }
        tscore=tscore+abs(tmscore);
      }
      if(tscore < fscore) {fscore=tscore; af=a;}
    }
    a=af;
    index=0;
    for(i=0;i<vmol.size();i++){
      tscore=0.0;
      for(j=0;j<vmol[i].nelectron/2;j++){
        if(score_type==0) tscore+=vmol[i].weight*(ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0);
        else if(score_type==1) tscore+=vmol[i].weight*abs(ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0);
        else if(score_type==2) tscore+=vmol[i].weight*pow((ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0),2.0);
        else if(score_type==4) tscore+=vmol[i].weight*pow((ev[index]-a*vmol[i].evref[j]-(1.0-a)*vmol[i].b0),4.0);
        index++;
      }
      vmol[i].error=tscore;
    }
    evscaling=af;
    delete  [] ev;
  }else if (fit_type==2){
    tscore=0.0;
    for(i=0;i<vmol.size();i++){
      tscore+=vmol[i].score(id,score_type,dftbversion,true);
    }
  }
}

void Erepobj::get_residual() {
  double tmp;
  restotU=0; restotS=0; restot2=0; restot4=0;

}

