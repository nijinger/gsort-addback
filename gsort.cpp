#include "gsort.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <TVector3.h>
#include <Math/Interpolator.h>
#include <ctime>

typedef ROOT::Math::Interpolator interp;

using namespace std;

interp *zner;
interp *znre;
interp *caer;
interp *care;

Double_t cluster_angle = 25/180.*TMath::Pi();

#include "masslcut.cpp"
#include "massrcut.cpp"

#define P0 3886.37 // updated 5/10/2018
//#define P0 3406.37 // updated 5/10/2018
//#define P0 3841.2
#define Mtot 118.0
#define beamE 168.8 // updated 5/10/2018 beamE=175,target thickness 1mg/cm2
//#define beamE 165.0
#define beamA 48.0
#define targetA 70.0
#define QNum 7
#define GTPosOffset 19.4
#define atomic_mass 931.4940954 // unit MeV/c2
#define v_c 299792458000 // mm/s

Double_t angle(Double_t x1,Double_t y1, Double_t z1,Double_t x2,Double_t y2,Double_t z2){
    return TVector3(x1,y1,z1).Angle(TVector3(x2,y2,z2));
}

bool BuildTable(const char*filename,interp* &er,interp* &re){
  ifstream in(filename);
  if(!in.is_open()) {
    cout<<"Cannot open file: "<<filename<<endl;
    return false;
  }
  Double_t energy[200],range[200];
  int cnt = 0;
  Double_t x,y,z;
  while(1){
    in>>energy[cnt]>>range[cnt];
    if(!in.good()){
      energy[cnt]=0;
      range[cnt]=0;
      break;
    }
    cnt++;
    if(cnt>=200) break;
  }
  
  er = new interp(cnt,ROOT::Math::Interpolation::kCSPLINE);
  re = new interp(cnt,ROOT::Math::Interpolation::kCSPLINE);
  er->SetData(cnt,energy,range);
  re->SetData(cnt,range,energy);
  return true;
}

Double_t updatep(Double_t pin,Double_t thick, int index=0){
  if(pin<=0) return 0;
  interp* er,*re;
  Double_t mass;
  if(index==0) {
    er = zner;
    re = znre;
    mass = 70;
  }else{
    er = caer;
    re = care;
    mass = 48;
  }
  Double_t ein = pin*pin*0.5/(mass*atomic_mass);
  Double_t range = er->Eval(ein)-thick;
  ein = re->Eval(range);
  return sqrt(2*mass*atomic_mass*ein);
}

Float_t offGT[QNum*4],gainGT[QNum*4];

template <class T,int m>
void readp(const char* fn, T (& off)[m], T (& gain)[m]){
  ifstream in(fn);
  if(!in.is_open()){
    cout<<"cannot open file: "<<fn<<endl;
    return;
  }
  int cnt=0;
  float x,y;
  while(1){
    in>>x>>y;
    if(!in.good()) break;
    off[cnt]=x;gain[cnt]=y;
    cnt++;
    if(cnt>=m) break;
  }
}

void gsort::Loop(TTree *opt)
{
//   In a ROOT session, you can do:
//      Root > .L gsort.C
//      Root > gsort t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;
   if (opt==0) return;
   BranchOpt(opt);

   Long64_t nentries = fChain->GetEntriesFast();

   std::clock_t start = std::clock();
   Double_t duration;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
       Clear();
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      if(!GetMass()) continue;
      gammadc();
      gamma_addback_dc();
//      gamma_bddback_dc();
      tkgammadc();
      opt->Fill();
      // if (Cut(ientry) < 0) continue;
      if(jentry%10000==0){
          duration = (std::clock() - start) / (double) CLOCKS_PER_SEC;
          cout<<"Process Rate: "<<jentry/duration<<"\t total:"<<jentry<<"\r";
          cout<<flush;
      }
   }
   cout<<"\n";
}

gsort::gsort(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("rootfile/run0009.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("rootfile/run0009.root");
      }
      f->GetObject("t",tree);

   }
   Init(tree);
   masslcut=load_masslcut();
   massrcut=load_massrcut();
   readp("GT_112014_2.cal",offGT,gainGT);
   BuildTable("calcium_range.txt",caer,care);
   BuildTable("zinc_range.txt",zner,znre);
}

gsort::~gsort()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t gsort::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t gsort::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void gsort::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("np", &np, &b_np);
   fChain->SetBranchAddress("cts", cts, &b_cts);
   fChain->SetBranchAddress("id", id, &b_id);
   fChain->SetBranchAddress("dT", dT, &b_dT);
   fChain->SetBranchAddress("dL", dL, &b_dL);
   fChain->SetBranchAddress("dR", dR, &b_dR);
   fChain->SetBranchAddress("phiL", phiL, &b_phiL);
   fChain->SetBranchAddress("phiR", phiR, &b_phiR);
   fChain->SetBranchAddress("fphiL", fphiL, &b_fphiL);
   fChain->SetBranchAddress("fphiR", fphiR, &b_fphiR);
   fChain->SetBranchAddress("thetaL", thetaL, &b_thetaL);
   fChain->SetBranchAddress("thetaR", thetaR, &b_thetaR);
   fChain->SetBranchAddress("fthetaL", fthetaL, &b_fthetaL);
   fChain->SetBranchAddress("fthetaR", fthetaR, &b_fthetaR);
   fChain->SetBranchAddress("ng", &ng, &b_ng);
   fChain->SetBranchAddress("nint", nint, &b_nint);
   fChain->SetBranchAddress("x", x, &b_x);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("z", z, &b_z);
   fChain->SetBranchAddress("e", e, &b_e);
   fChain->SetBranchAddress("xm", xm, &b_xm);
   fChain->SetBranchAddress("ym", ym, &b_ym);
   fChain->SetBranchAddress("zm", zm, &b_zm);
   fChain->SetBranchAddress("xs", xs, &b_xs);
   fChain->SetBranchAddress("ys", ys, &b_ys);
   fChain->SetBranchAddress("zs", zs, &b_zs);
   fChain->SetBranchAddress("ee1", ee1, &b_ee1);
   fChain->SetBranchAddress("ee2", ee2, &b_ee2);
   fChain->SetBranchAddress("eem", eem, &b_eem);
   fChain->SetBranchAddress("dx1", dx1, &b_dx1);
   fChain->SetBranchAddress("dy1", dy1, &b_dy1);
   fChain->SetBranchAddress("dz1", dz1, &b_dz1);
   fChain->SetBranchAddress("dx2", dx2, &b_dx2);
   fChain->SetBranchAddress("dy2", dy2, &b_dy2);
   fChain->SetBranchAddress("dz2", dz2, &b_dz2);
   fChain->SetBranchAddress("dxm", dxm, &b_dxm);
   fChain->SetBranchAddress("dym", dym, &b_dym);
   fChain->SetBranchAddress("dzm", dzm, &b_dzm);
   fChain->SetBranchAddress("dis1", dis1, &b_dis1);
   fChain->SetBranchAddress("dis2", dis2, &b_dis2);
   fChain->SetBranchAddress("dism", dism, &b_dism);
   fChain->SetBranchAddress("de1", de1, &b_de1);
   fChain->SetBranchAddress("de2", de2, &b_de2);
   fChain->SetBranchAddress("dem", dem, &b_dem);
   fChain->SetBranchAddress("theta", theta, &b_theta);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("t", t, &b_t);
   fChain->SetBranchAddress("t0", t0, &b_t0);
   fChain->SetBranchAddress("cc_id", cc_id, &b_cc_id);
   fChain->SetBranchAddress("cry_id", cry_id, &b_cry_id);
   fChain->SetBranchAddress("ntg", &ntg, &b_ntg);
   fChain->SetBranchAddress("tsTK", &tsTK, &b_tsTK);
   fChain->SetBranchAddress("pad", pad, &b_pad);
   fChain->SetBranchAddress("tracked", tracked, &b_tracked);
   fChain->SetBranchAddress("esum", esum, &b_esum);
   fChain->SetBranchAddress("ndet", ndet, &b_ndet);
   fChain->SetBranchAddress("gtkts", gtkts, &b_gtkts);
   fChain->SetBranchAddress("fom", fom, &b_fom);
   fChain->SetBranchAddress("x0", x0, &b_x0);
   fChain->SetBranchAddress("y0", y0, &b_y0);
   fChain->SetBranchAddress("z0", z0, &b_z0);
   fChain->SetBranchAddress("e0", e0, &b_e0);
   fChain->SetBranchAddress("x1", x1, &b_x1);
   fChain->SetBranchAddress("y1", y1, &b_y1);
   fChain->SetBranchAddress("z1", z1, &b_z1);
   fChain->SetBranchAddress("e1", e1, &b_e1);
   Notify();
}

Bool_t gsort::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void gsort::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t gsort::Cut(Long64_t entry)
{
// This function may1 be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}

void gsort::BranchOpt(TTree *opt){
  assert(opt!=0);

  opt->Branch("MassT",&MassT,"MassT/D");
  opt->Branch("MassP",&MassP,"MassP/D");
  opt->Branch("bt",&bt,"bt/D");
  opt->Branch("bp",&bp,"bp/D");
  opt->Branch("theT",&theT,"theT/D");
  opt->Branch("phiT",&phiT,"phiT/D");
  opt->Branch("theP",&theP,"theP/D");
  opt->Branch("phiP",&phiP,"phiP/D");
  opt->Branch("qval",&qval,"qval/D");
  opt->Branch("qval2",&qval2,"qval2/D");
  opt->Branch("np",&np,"np/I");
  opt->Branch("cts",&cts,"cts[np]/L");

  opt->Branch("ng",&ng,"ng/I");
  opt->Branch("t",&t,"t[ng]/D");
  opt->Branch("t0",&t0,"t0[ng]/D");
  opt->Branch("x",&x,"x[ng]/F");
  opt->Branch("y",&y,"y[ng]/F");
  opt->Branch("z",&z,"z[ng]/F");
  opt->Branch("ge",&ge,"ge[ng]/D");
  opt->Branch("get",&get,"get[ng]/D");
  opt->Branch("gep",&gep,"gep[ng]/D");
  opt->Branch("gthe",&gthe,"gthe[ng]/D");
  opt->Branch("gphi",&gphi,"gphi[ng]/D");
  opt->Branch("gthet",&gthet,"gthet[ng]/D");
  opt->Branch("gthep",&gthep,"gthep[ng]/D");

  opt->Branch("ntg",&ntg,"ntg/I");
  opt->Branch("x0",&x0,"x0[ntg]/F");
  opt->Branch("y0",&y0,"y0[ntg]/F");
  opt->Branch("z0",&z0,"z0[ntg]/F");
  opt->Branch("e0",&e0,"e0[ntg]/F");
  opt->Branch("x1",&x1,"x1[ntg]/F");
  opt->Branch("y1",&y1,"y1[ntg]/F");
  opt->Branch("z1",&z1,"z1[ntg]/F");
  opt->Branch("e1",&e1,"e1[ntg]/F");
  opt->Branch("ndet",&ndet,"ndet[ntg]/I");
  opt->Branch("fom",&fom,"fom[ntg]/F");
  opt->Branch("gtkts",&gtkts,"gtkts[ntg]/L");
  opt->Branch("tge",&tge,"tge[ntg]/D");
  opt->Branch("tget",&tget,"tget[ntg]/D");
  opt->Branch("tgep",&tgep,"tgep[ntg]/D");
  opt->Branch("tgthe",&tgthe,"tgthe[ntg]/D");
  opt->Branch("tgphi",&tgphi,"tgphi[ntg]/D");
  opt->Branch("tgthet",&tgthet,"tgthet[ntg]/D");
  opt->Branch("tgthep",&tgthep,"tgthep[ntg]/D");

  opt->Branch("nag",&nag,"nag/I");
  opt->Branch("anint",&anint,"anint[nag]/I");
  opt->Branch("ax1",&ax1,"ax1[nag]/D");
  opt->Branch("ay1",&ay1,"ay1[nag]/D");
  opt->Branch("az1",&az1,"az1[nag]/D");
  opt->Branch("ax2",&ax2,"ax2[nag]/D");
  opt->Branch("ay2",&ay2,"ay2[nag]/D");
  opt->Branch("az2",&az2,"az2[nag]/D");
  opt->Branch("age",&age,"age[nag]/D");
  opt->Branch("aget",&aget,"aget[nag]/D");
  opt->Branch("agep",&agep,"agep[nag]/D");
  opt->Branch("agthe",&agthe,"agthe[nag]/D");
  opt->Branch("agphi",&agphi,"agphi[nag]/D");
  opt->Branch("agthet",&agthet,"agthet[nag]/D");
  opt->Branch("agthep",&agthep,"agthep[nag]/D");
  opt->Branch("openangle",&openangle,"openangle[nag]/D");
  opt->Branch("distance",&distance,"distance[nag]/D");
  opt->Branch("aggthe",&aggthe,"aggthe[nag]/D");
  opt->Branch("aggphi",&aggphi,"aggphi[nag]/D");

  opt->Branch("nbg",&nbg,"nbg/I");
  opt->Branch("bnint",&bnint,"bnint[nbg]/I");
  opt->Branch("bge",&bge,"bge[nbg]/D");
  opt->Branch("bget",&bget,"bget[nbg]/D");
  opt->Branch("bgep",&bgep,"bgep[nbg]/D");
  opt->Branch("bgthe",&bgthe,"bgthe[nbg]/D");
  opt->Branch("bgphi",&bgphi,"bgphi[nbg]/D");
  opt->Branch("bgthet",&bgthet,"bgthet[nbg]/D");
  opt->Branch("bgthep",&bgthep,"bgthep[nbg]/D");
  opt->Branch("bopenangle",&bopenangle,"bopenangle[nbg]/D");
  opt->Branch("bdistance",&bdistance,"bdistance[nbg]/D");
  opt->Branch("bggthe",&bggthe,"bggthe[nbg]/D");
  opt->Branch("bggphi",&bggphi,"bggphi[nbg]/D");

}

bool gsort::GetMass(){
  float theta1=0.0,theta2=0.0;
  float phi1=0,phi2=0;
  float L1=0.0,L2=0.0;
  float dt=0;
  float f1,f2,p1,p2;
  if(masslcut->IsInside(fthetaL[0],dT[0])){
    theta1 = fthetaL[0];
    theta2 = fthetaR[0];
    phi1=fphiL[0];
    phi2=fphiR[0];
    L1 = dL[0];
    L2 = dR[0];
//    dt = (dT[0]-2.5)*0.1; // Shaofei's code
    dt = (dT[0]-0.0)*0.1;
  }else if(massrcut->IsInside(fthetaR[0],dT[0])){
    theta1 = fthetaR[0];
    theta2 = fthetaL[0];
    phi1 = fphiR[0];
    phi2 = fphiL[0];
    L1 = dR[0];
    L2 = dL[0];
//    dt = -1.*(dT[0]-2.5)*0.1;
    dt = -1.*(dT[0]-0.0)*0.1;
  }else return false;
//  f1 = 1.1094 - 0.12008/TMath::Cos(theta1) - 0.0096272/TMath::Power(TMath::Cos(theta1),2);
//  f1*=1.06;
//  f2 = 1.96556 - 4.3451*TMath::Cos(theta2) + 4.2348*TMath::Power(TMath::Cos(theta2),2);
  p1 = P0 * TMath::Sin(theta2)/TMath::Sin(theta1+theta2);//*f1;
  p2 = P0 * TMath::Sin(theta1)/TMath::Sin(theta1+theta2);//*f2;
  
//  MassT = (-0.032206*dt + Mtot*L1/p1)/(L1/p1+L2/p2);
//  MassP = (0.032206*dt + Mtot*L2/p2)/(L1/p1+L2/p2);

  L1 = L1 * 10;
  L2 = L2 * 10;
//  p1 = updatep(p1,0.70028011/cos(theta1),1);
//  p2 = updatep(p2,0.70028011/cos(theta2),0);
  p1 = updatep(p1,0.45028011/cos(theta1),1);
  p2 = updatep(p2,0.45028011/cos(theta2),0);
  double temp = (L1/p2/v_c+L2/p1/v_c)*atomic_mass;
  MassP = (dt*1e-9+L2/p2*Mtot*atomic_mass/v_c)/temp;
  MassT = (-dt*1e-9+L1/p1*Mtot*atomic_mass/v_c)/temp;

  bt = p2/MassT/atomic_mass;
  bp = p1/MassP/atomic_mass;
//  bt = p2/targetA/atomic_mass;
//  bp = p1/beamA/atomic_mass;

  Double_t et,et2,ep;
  ep = TMath::Power(TMath::Sin(theta2)/TMath::Sin(theta1+theta2),2)*beamE;
//  bp = TMath::Sqrt(ep/beamA)*0.046;
//  bp *= f1;
  et = TMath::Power(TMath::Sin(theta1)/TMath::Sin(theta1+theta2),2)*beamE*beamA/targetA;
  et2 = TMath::Power(TMath::Sin(theta1)/TMath::Sin(theta1+theta2),2)*beamE*48./70.;
//  bt = sqrt(et/targetA)*0.046;
//  bt *=f2;

  theP = theta1;
  phiP = phi1;
  theT = theta2;
  phiT = phi2;

  qval = ep+et-beamE;
  qval2 = ep+et2-beamE;
  return true;
}

void gsort::Clear(){
  MassT=0;
  MassP=0;
  bt=0;
  bp=0;
  qval=0;
  memset(tget,0,sizeof(tget));
  nag=0;
  nbg=0;
}

void gsort::gammadc(){
  if(ng<=0) return;
  for(int i=0;i<ng;i++){
    int gid = (cc_id[i]-1)*4+cry_id[i];
    assert(gid<28);
    ge[i]=e[i]*gainGT[gid]+offGT[gid];
    TVector3 gp(x[i],y[i],z[i]);
    TVector3 pp,pt;
    pt.SetMagThetaPhi(1,theT,phiT);
    pp.SetMagThetaPhi(1,theP,phiP);
    gthe[i]=gp.Angle(TVector3(0,0,1));
    gphi[i]=gp.Phi();
    gthet[i]=gp.Angle(pt);
    gthep[i]=gp.Angle(pp);
    get[i]=ge[i]*(1-bt*TMath::Cos(gthet[i]))/TMath::Sqrt(1-bt*bt);
    gep[i]=ge[i]*(1-bp*TMath::Cos(gthep[i]))/TMath::Sqrt(1-bp*bp);
  }
}

void gsort::gamma_addback_dc(){
    if(ng<=0) return;
    Int_t *index = new Int_t[ng];
    TMath::Sort(ng,ee1,index); // sorted by the 1st interaction point

    if(ng<2) {
        nag = 1;
        anint[0] = nint[0];
        age[0]=ge[0];
        ax1[0] = x[0];
        ay1[0] = y[0];
        az1[0] = z[0];
        ax2[0] = xs[0];
        ay2[0] = ys[0];
        az2[0] = zs[0];
        openangle[0] = 0;
        distance[0] = 0;
    } else
    for(int i =0; i<ng-1; i++){
        TVector3 dr1(xm[index[i]],ym[index[i]],zm[index[i]]);
        TVector3 dr2(x[index[i+1]],y[index[i+1]],z[index[i+1]]);
        openangle[nag] = dr1.Angle(dr2);
        distance[nag] = (dr1-dr2).Mag();
        
//        if(openangle[nag]>cluster_angle || eem[index[i+1]] > ee1[index[i]]){
        if(openangle[nag]>cluster_angle){
            anint[nag] = nint[index[i]];
            age[nag] = ge[index[i]];
            ax1[nag] = x[index[i]];
            ay1[nag] = y[index[i]];
            az1[nag] = z[index[i]];
            ax2[nag] = xs[index[i]];
            ay2[nag] = ys[index[i]];
            az2[nag] = zs[index[i]];
            nag++;
            if(i == ng - 2){
               i++;
               anint[nag] = nint[index[i]];
               age[nag] = ge[index[i]];
               ax1[nag] = x[index[i]];
               ay1[nag] = y[index[i]];
               az1[nag] = z[index[i]];
               ax2[nag] = xs[index[i]];
               ay2[nag] = ys[index[i]];
               az2[nag] = zs[index[i]];
               nag++;
            }
        }else{
            if( i <= ng - 3 && angle(xm[index[i]],ym[index[i]],zm[index[i]],x[index[i+2]],y[index[i+2]],z[index[i+2]]) < cluster_angle && angle(xm[index[i+1]],ym[index[i+1]],zm[index[i+1]],x[index[i+2]],y[index[i+2]],z[index[i+2]]) < cluster_angle){ // n3 events
                age[nag] = ge[index[i]]+ge[index[i+1]]+ge[index[i+2]];
                anint[nag] = nint[index[i]]+nint[index[i+1]]+nint[index[i+2]];
                ax1[nag] = x[index[i]];
                ay1[nag] = y[index[i]];
                az1[nag] = z[index[i]];
                if(nint[index[i]]<2){
                  ax2[nag] = x[index[i+1]];
                  ay2[nag] = y[index[i+1]];
                  az2[nag] = z[index[i+1]];
                }else{
                  ax2[nag] = xs[index[i]];
                  ay2[nag] = ys[index[i]];
                  az2[nag] = zs[index[i]];
                }
                nag++;
                i+=2;
            }else{
    
                age[nag]=ge[index[i]]+ge[index[i+1]];
                anint[nag] = nint[index[i]]+nint[index[i+1]];
                ax1[nag] = x[index[i]];
                ay1[nag] = y[index[i]];
                az1[nag] = z[index[i]];
                if(nint[index[i]]<2){
                  ax2[nag] = x[index[i+1]];
                  ay2[nag] = y[index[i+1]];
                  az2[nag] = z[index[i+1]];
                }else{
                  ax2[nag] = xs[index[i]];
                  ay2[nag] = ys[index[i]];
                  az2[nag] = zs[index[i]];
                }
                nag++;
                i++;
            }
        }
    }

    delete [] index;

    TVector3 pp,pt;
    pp.SetMagThetaPhi(1,theP,phiP);
    pt.SetMagThetaPhi(1,theT,phiT);
    for(int i = 0;i<nag;i++){
        TVector3 gg(ax1[i],ay1[i],az1[i]);
        agthe[i] = gg.Angle(TVector3(0,0,1));
        agphi[i] = gg.Phi();
        agthet[i] = gg.Angle(pt);
        agthep[i] = gg.Angle(pp);
        aget[i]=age[i]*(1-bt*TMath::Cos(agthet[i]))/TMath::Sqrt(1-bt*bt);
        agep[i]=age[i]*(1-bp*TMath::Cos(agthep[i]))/TMath::Sqrt(1-bp*bp);

        if(anint[i]>1){
            TVector3 gg2(ax2[i]-ax1[i],ay2[i]-ay1[i],az2[i]-az1[i]);
            aggthe[i] = gg.Angle(gg2);
            TVector3 reaction = pt.Cross(gg);
            TVector3 compton = gg.Cross(gg2);
            aggphi[i] = reaction.Angle(compton);
        }else{
            aggthe[i] = -1;
            aggphi[i] = -1;
        }
    }
}

void gsort::gamma_bddback_dc(){
    if(ng<=0) return;
    Int_t *index = new Int_t[ng];
    TMath::Sort(ng,dis1,index,false);

    if(ng<2) {
        nbg = 1;
        bnint[0] = nint[0];
        bge[0]=ge[0];
        bx1[0] = dx1[0];
        by1[0] = dy1[0];
        bz1[0] = dz1[0];
        bx2[0] = dx2[0];
        by2[0] = dy2[0];
        bz2[0] = dz2[0];
        bopenangle[0] = 0;
        bdistance[0] = 0;
    } else
    for(int i =0; i<ng-1; i++){
        TVector3 dr1(dxm[index[i]],dym[index[i]],dzm[index[i]]);
        TVector3 dr2(dx1[index[i+1]],dy1[index[i+1]],dz1[index[i+1]]);
        bopenangle[nbg] = dr1.Angle(dr2);
        bdistance[nbg] = (dr1-dr2).Mag();
        
        if(bopenangle[nbg]>cluster_angle){
            bnint[nbg] = nint[index[i]];
            bge[nbg] = ge[index[i]];
            bx1[nbg] = dx1[index[i]];
            by1[nbg] = dy1[index[i]];
            bz1[nbg] = dz1[index[i]];
            bx2[nbg] = dx2[index[i]];
            by2[nbg] = dy2[index[i]];
            bz2[nbg] = dz2[index[i]];
            nbg++;
            if(i == ng - 2){
               i++;
               bnint[nbg] = nint[index[i]];
               bge[nbg] = ge[index[i]];
               bx1[nbg] = dx1[index[i]];
               by1[nbg] = dy1[index[i]];
               bz1[nbg] = dz1[index[i]];
               bx2[nbg] = dx2[index[i]];
               by2[nbg] = dy2[index[i]];
               bz2[nbg] = dz2[index[i]];
               nbg++;
            }
        }else{
            if( i <= ng-3 && angle(dxm[index[i]],dym[index[i]],dzm[index[i]],dx1[index[i+2]],dy1[index[i+2]],dz1[index[i+2]]) < cluster_angle && angle(dxm[index[i+1]],dym[index[i+1]],dzm[index[i+1]],dx1[index[i+2]],dy1[index[i+2]],dz1[index[i+2]]) < cluster_angle){
                bge[nbg] = ge[index[i]]+ge[index[i+1]]+ge[index[i+2]];
                bnint[nbg] = nint[index[i]]+nint[index[i+1]]+nint[index[i+2]];
                bx1[nbg] = x[index[i]];
                by1[nbg] = y[index[i]];
                bz1[nbg] = z[index[i]];
                if(nint[index[i]]<2){
                  bx2[nbg] = dx1[index[i+1]];
                  by2[nbg] = dy1[index[i+1]];
                  bz2[nbg] = dz1[index[i+1]];
                }else{
                  bx2[nbg] = dx2[index[i]];
                  by2[nbg] = dy2[index[i]];
                  bz2[nbg] = dz2[index[i]];
                }
                nbg++;
                i+=2;
            } else {
                bge[nbg]=ge[index[i]]+ge[index[i+1]];
                bnint[nbg] = nint[index[i]]+nint[index[i+1]];
                bx1[nbg] = dx1[index[i]];
                by1[nbg] = dy1[index[i]];
                bz1[nbg] = dz1[index[i]];
                if(nint[index[i]]<2){
                  bx2[nbg] = dx1[index[i+1]];
                  by2[nbg] = dy1[index[i+1]];
                  bz2[nbg] = dz1[index[i+1]];
                }else{
                  bx2[nbg] = dx2[index[i]];
                  by2[nbg] = dy2[index[i]];
                  bz2[nbg] = dz2[index[i]];
                }
                nbg++;
                i++;
            }
        }
    }

    delete [] index;

    TVector3 pp,pt;
    pp.SetMagThetaPhi(1,theP,phiP);
    pt.SetMagThetaPhi(1,theT,phiT);
    for(int i = 0;i<nbg;i++){
        TVector3 gg(bx1[i],by1[i],bz1[i]);
        bgthe[i] = gg.Angle(TVector3(0,0,1));
        bgphi[i] = gg.Phi();
        bgthet[i] = gg.Angle(pt);
        bgthep[i] = gg.Angle(pp);
        bget[i]=bge[i]*(1-bt*TMath::Cos(bgthet[i]))/TMath::Sqrt(1-bt*bt);
        bgep[i]=bge[i]*(1-bp*TMath::Cos(bgthep[i]))/TMath::Sqrt(1-bp*bp);

        if(bnint[i]>1){
            TVector3 gg2(bx2[i]-bx1[i],by2[i]-by1[i],bz2[i]-bz1[i]);
            bggthe[i] = gg.Angle(gg2);
            TVector3 reaction = pt.Cross(gg);
            TVector3 compton = gg.Cross(gg2);
            bggphi[i] = reaction.Angle(compton);
        }else{
            bggthe[i] = -1;
            bggphi[i] = -1;
        }
    }
}

void gsort::tkgammadc(){
  if(ntg<=0) return;
  TVector3 pp,pt;
  pt.SetMagThetaPhi(1,theT,phiT);
  pp.SetMagThetaPhi(1,theP,phiP);
  for(int i=0;i<ntg;i++){
    tge[i]=esum[i];
    TVector3 gp1(x0[i],y0[i]-GTPosOffset,z0[i]);
    TVector3 gp2(x1[i]-x0[i],y1[i]-y0[i],z1[i]-z0[i]);
    tgthe[i]=gp1.Angle(TVector3(0,0,1));
    tgphi[i]=gp1.Phi();
    tgthet[i]=gp1.Angle(pt);
    tgthep[i]=gp1.Angle(pp);
    tget[i]=tge[i]*(1-bt*TMath::Cos(tgthet[i]))/TMath::Sqrt(1-bt*bt);
    tgep[i]=tge[i]*(1-bp*TMath::Cos(tgthep[i]))/TMath::Sqrt(1-bp*bp);
  }
}
