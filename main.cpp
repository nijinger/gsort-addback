#include <iostream>
#include <TSystem.h>

#include "gsort.h"

using namespace std;

int main(int argc,char **argv){
  if(argc!=3) {
//    cout<<argv[0]<<" # runnum"<<endl;
    cout<<argv[0]<<" INPUTFILE OUTPUTFILE"<<endl;
    return 1;
  }
//  int runnum = atoi(argv[1]);
//  TFile *ipf=new TFile(Form("./rootfile/run%04d.root",runnum));
  if(gSystem->AccessPathName(argv[1])) return 0;
  TFile *ipf=new TFile(argv[1]);
  if(!ipf->IsOpen()) return 0;
  TTree *ipt = (TTree *)ipf->Get("t");
  gsort *sr=new gsort(ipt);
  
//  TFile *opf=new TFile(Form("./gsortfile/run%04d.root",runnum),"RECREATE");
  TFile *opf=new TFile(argv[2],"RECREATE");
  TTree *opt=new TTree("t","t");
  sr->Loop(opt);
  opt->Write();
  opf->Close();
  ipf->Close();
  
  return 0;
}
