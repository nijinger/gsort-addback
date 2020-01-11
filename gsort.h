//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Sun Dec 10 10:08:20 2017 by ROOT version 5.34/36
// from TTree t/opt
// found on file: rootfile/run0009.root
//////////////////////////////////////////////////////////

#ifndef gsort_h
#define gsort_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TCutG.h>
#include <TMath.h>
#include <assert.h>

#define MAX_G_N 42 // maiimum gamma number

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class gsort {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           np;
   Long64_t        cts[1];   //[np]
   Int_t           id[1];   //[np]
   Int_t           dT[1];   //[np]
   Float_t         dL[1];   //[np]
   Float_t         dR[1];   //[np]
   Int_t           phiL[1];   //[np]
   Int_t           phiR[1];   //[np]
   Float_t         fphiL[1];   //[np]
   Float_t         fphiR[1];   //[np]
   Int_t           thetaL[1];   //[np]
   Int_t           thetaR[1];   //[np]
   Float_t         fthetaL[1];   //[np]
   Float_t         fthetaR[1];   //[np]

   Int_t           ng;
   Int_t           nint[MAX_G_N];
   Float_t         x[MAX_G_N];   //[ng]
   Float_t         y[MAX_G_N];   //[ng]
   Float_t         z[MAX_G_N];   //[ng]
   Float_t         e[MAX_G_N];   //[ng]
   Float_t         xm[MAX_G_N];   //[ng]
   Float_t         ym[MAX_G_N];   //[ng]
   Float_t         zm[MAX_G_N];   //[ng]
   Float_t         xs[MAX_G_N];   //[ng]
   Float_t         ys[MAX_G_N];   //[ng]
   Float_t         zs[MAX_G_N];   //[ng]
   Float_t         ee1[MAX_G_N];
   Float_t         ee2[MAX_G_N];
   Float_t         eem[MAX_G_N];
   Float_t         dx1[MAX_G_N];   //[ng]
   Float_t         dy1[MAX_G_N];   //[ng]
   Float_t         dz1[MAX_G_N];   //[ng]
   Float_t         dx2[MAX_G_N];   //[ng]
   Float_t         dy2[MAX_G_N];   //[ng]
   Float_t         dz2[MAX_G_N];   //[ng]
   Float_t         dxm[MAX_G_N];   //[ng]
   Float_t         dym[MAX_G_N];   //[ng]
   Float_t         dzm[MAX_G_N];   //[ng]
   Float_t         dis1[MAX_G_N];
   Float_t         dis2[MAX_G_N];
   Float_t         dism[MAX_G_N];
   Float_t         de1[MAX_G_N];
   Float_t         de2[MAX_G_N];
   Float_t         dem[MAX_G_N];
   Float_t         theta[MAX_G_N];   //[ng]
   Float_t         phi[MAX_G_N];   //[ng]
   Double_t        t[MAX_G_N];   //[ng]
   Double_t        t0[MAX_G_N];   //[ng]
   Int_t           cc_id[MAX_G_N];   //[ng]
   Int_t           cry_id[MAX_G_N];   //[ng]

   Int_t           ntg;
   Long64_t        tsTK;
   Int_t           pad[MAX_G_N];   //[ntg]
   Int_t           tracked[MAX_G_N];   //[ntg]
   Float_t         esum[MAX_G_N];   //[ntg]
   Int_t           ndet[MAX_G_N];   //[ntg]
   Long64_t        gtkts[MAX_G_N];   //[ntg]
   Float_t         fom[MAX_G_N];   //[ntg]
   Float_t         x0[MAX_G_N];   //[ntg]
   Float_t         y0[MAX_G_N];   //[ntg]
   Float_t         z0[MAX_G_N];   //[ntg]
   Float_t         e0[MAX_G_N];   //[ntg]
   Float_t         x1[MAX_G_N];   //[ntg]
   Float_t         y1[MAX_G_N];   //[ntg]
   Float_t         z1[MAX_G_N];   //[ntg]
   Float_t         e1[MAX_G_N];   //[ntg]

   // List of branches
   TBranch        *b_np;   //!
   TBranch        *b_cts;   //!
   TBranch        *b_id;   //!
   TBranch        *b_dT;   //!
   TBranch        *b_dL;   //!
   TBranch        *b_dR;   //!
   TBranch        *b_phiL;   //!
   TBranch        *b_phiR;   //!
   TBranch        *b_fphiL;   //!
   TBranch        *b_fphiR;   //!
   TBranch        *b_thetaL;   //!
   TBranch        *b_thetaR;   //!
   TBranch        *b_fthetaL;   //!
   TBranch        *b_fthetaR;   //!
   TBranch        *b_ng;   //!
   TBranch        *b_nint;
   TBranch        *b_x;   //!
   TBranch        *b_y;   //!
   TBranch        *b_z;   //!
   TBranch        *b_e;   //!
   TBranch        *b_xm;   //!
   TBranch        *b_ym;   //!
   TBranch        *b_zm;   //!
   TBranch        *b_xs;   //!
   TBranch        *b_ys;   //!
   TBranch        *b_zs;   //!
   TBranch        *b_ee1;   //!
   TBranch        *b_ee2;   //!
   TBranch        *b_eem;   //!
   TBranch        *b_dx1;   //!
   TBranch        *b_dy1;   //!
   TBranch        *b_dz1;   //!
   TBranch        *b_dx2;   //!
   TBranch        *b_dy2;   //!
   TBranch        *b_dz2;   //!
   TBranch        *b_dxm;   //!
   TBranch        *b_dym;   //!
   TBranch        *b_dzm;   //!
   TBranch        *b_dis1;   //!
   TBranch        *b_dis2;   //!
   TBranch        *b_dism;   //!
   TBranch        *b_de1;
   TBranch        *b_de2;
   TBranch        *b_dem;
   TBranch        *b_theta;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_t;   //!
   TBranch        *b_t0;   //!
   TBranch        *b_cc_id;   //!
   TBranch        *b_cry_id;   //!
   TBranch        *b_ntg;   //!
   TBranch        *b_tsTK;   //!
   TBranch        *b_pad;   //!
   TBranch        *b_tracked;   //!
   TBranch        *b_esum;   //!
   TBranch        *b_ndet;   //!
   TBranch        *b_gtkts;   //!
   TBranch        *b_fom;   //!
   TBranch        *b_x0;   //!
   TBranch        *b_y0;   //!
   TBranch        *b_z0;   //!
   TBranch        *b_e0;   //!
   TBranch        *b_x1;   //!
   TBranch        *b_y1;   //!
   TBranch        *b_z1;   //!
   TBranch        *b_e1;   //!

   TCutG   *masslcut;
   TCutG   *massrcut;
   // Output variables 
   Double_t   MassT;
   Double_t   MassP;
   Double_t   bt; // beta valuw assuming target
   Double_t   bp; // beta value assuming projectile
   Double_t   theT; // 
   Double_t   phiT;
   Double_t   theP;
   Double_t   phiP;
   Double_t   qval;
   Double_t   qval2;

   Double_t   ge[MAX_G_N];
   Double_t   get[MAX_G_N];
   Double_t   gep[MAX_G_N];
   Double_t   gthe[MAX_G_N];
   Double_t   gphi[MAX_G_N];
   Double_t   gthet[MAX_G_N];
   Double_t   gthep[MAX_G_N];

   Double_t   tge[MAX_G_N];
   Double_t   tget[MAX_G_N];
   Double_t   tgep[MAX_G_N];
   Double_t   tgthe[MAX_G_N];
   Double_t   tgphi[MAX_G_N];
   Double_t   tgthet[MAX_G_N];
   Double_t   tgthep[MAX_G_N];

   Int_t      nag;
   Int_t      anint[MAX_G_N];
   Double_t   age[MAX_G_N];
   Double_t   aget[MAX_G_N];
   Double_t   agep[MAX_G_N];
   Double_t   agthe[MAX_G_N];
   Double_t   agphi[MAX_G_N];
   Double_t   agthet[MAX_G_N];
   Double_t   agthep[MAX_G_N];
   Double_t   ax1[MAX_G_N]; // based on energy order
   Double_t   ay1[MAX_G_N];
   Double_t   az1[MAX_G_N];
   Double_t   ax2[MAX_G_N]; 
   Double_t   ay2[MAX_G_N];
   Double_t   az2[MAX_G_N];
   Double_t   openangle[MAX_G_N];
   Double_t   distance[MAX_G_N];
   Double_t   aggthe[MAX_G_N]; // gamma gamma angle
   Double_t   aggphi[MAX_G_N]; // angle between comton scattering plane and reaction plane

   Int_t      nbg;
   Int_t      bnint[MAX_G_N];
   Double_t   bge[MAX_G_N];
   Double_t   bget[MAX_G_N];
   Double_t   bgep[MAX_G_N];
   Double_t   bgthe[MAX_G_N];
   Double_t   bgphi[MAX_G_N];
   Double_t   bgthet[MAX_G_N];
   Double_t   bgthep[MAX_G_N];
   Double_t   bx1[MAX_G_N]; // based on energy order
   Double_t   by1[MAX_G_N];
   Double_t   bz1[MAX_G_N];
   Double_t   bx2[MAX_G_N]; 
   Double_t   by2[MAX_G_N];
   Double_t   bz2[MAX_G_N];
   Double_t   bopenangle[MAX_G_N];
   Double_t   bdistance[MAX_G_N];
   Double_t   bggthe[MAX_G_N]; // gamma gamma angle
   Double_t   bggphi[MAX_G_N]; // angle between comton scattering plane and reaction plane

   gsort(TTree *tree=0);
   virtual ~gsort();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TTree *opt);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);

   void     BranchOpt(TTree *opt);
   bool     GetMass();
   void     Clear();
   void     gammadc(); // gamma dop. correc.
   void     tkgammadc(); // tracked gamma dop. correc.
   void     gamma_addback_dc();
   void     gamma_bddback_dc();
};

#endif

