//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Feb  8 13:19:07 2016 by ROOT version 5.34/09
// from TTree tt/AMPT tree
// found on file: user.phuo.7638948._000002.ampt.root
//////////////////////////////////////////////////////////

#ifndef UUana_h
#define UUana_h

#include <TROOT.h>
#include <TChain.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <cmath>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TPad.h>
#include <iostream>
#include <fstream>
#include <string>

const int MaxTrack = 100000;

class UUana {
    public :

        //User define
        char name[200];
        enum{
            NHAR = 6,
            PBIN=200,
            EBIN=240,
            NAUM=238,
        };

        UUana(string filelist, int fr, int tr);

        void InitHist();
        bool calcEcc();
        void eventInfo();
        void looptrack2();

//        int getphibin(float mphi);
//        int getetabin(float meta);
//        void makeCan();

//        Short_t phibin[PBIN];
//        Short_t etabin[EBIN];
        int evtno;
        TTree* outtree;
        TCanvas* can;

        bool fill_tree;
        bool flag;

        Int_t mevts;

        int from;
        int to;
    
    
        int npartP;
        int npartT;
        int Np_neutron;
        int Nt_neutron;
        int N_neutron;

        float qnx[3][NHAR];
        float qny[3][NHAR];
        float qn[3][NHAR];
        float qn0[3];
        float ecc[NHAR];
        float eccP[NHAR];
        float eccT[NHAR];
    
        float eccangP[NHAR];
        float eccangT[NHAR];
        float eccang[NHAR];
        float psiang[NHAR];


        float mrawvnc[NHAR];
        float mrawvns[NHAR];


        int eff_trk;
        int eff_npart;


        TFile* fout;
//        TFile* fout2;

        TH2* hecc_psi[NHAR];
        TH2* hecc_vn[NHAR];

        TH1* h_ecc[NHAR];
        TH1* h_qn[3][NHAR];
        TH1* h_psin[3][NHAR];
        TH2* h_qnxy[3][NHAR];
        TH2* h_qndxy[NHAR];
    
        TH2* h_spectator;
        TH1* h_gausP;
        TH1* h_gausT;
        TH2* h_ZDC;
    
        TH1* h_vnc[NHAR];
        TH1* h_vns[NHAR];

        TH1* h_pt;
        TH1* h_phi;
        TH1* h_eta;
        TH1* h_bimp;
        TH1* h_mul;
        TH1* h_mulraw;
        TH1* h_multotal;
        TH1* h_npart;
        TH2* h_phiPT;
        TH2* h_thetaPT;
        TH2* h_etaphi;


        //end of user define

        TTree          *fChain;   //!pointer to the analyzed TTree or TChain
        Int_t           fCurrent; //!current Tree number in a TChain

        int badtrk;
        // Declaration of leaf types
        Int_t           ieve;
        Int_t           jeve;
        Float_t         bimp;
        Int_t           eid[6];
        Int_t           ntrk;
        Int_t           id[MaxTrack];   //[ntrk]
        Float_t         phi[MaxTrack];   //[ntrk]
        Float_t         pt[MaxTrack];   //[ntrk]
        Float_t         eta[MaxTrack];   //[ntrk]
        Float_t         ms[MaxTrack];   //[ntrk]
        Int_t           npart;
        Float_t         thetaP;
        Float_t         phiP;
        Float_t         thetaT;
        Float_t         phiT;
        Float_t         x_part[476];   //[npart]
        Float_t         y_part[476];   //[npart]
        Int_t           id_part[476];   //[npart]
        Int_t           st_part[476];   //[npart]
        Int_t           pid_part[476];   //[npart]

        // List of branches
        TBranch        *b_ieve;   //!
        TBranch        *b_jeve;   //!
        TBranch        *b_bimp;   //!
        TBranch        *b_eid;   //!
        TBranch        *b_ntrk;   //!
        TBranch        *b_id;   //!
        TBranch        *b_phi;   //!
        TBranch        *b_pt;   //!
        TBranch        *b_eta;   //!
        TBranch        *b_ms;   //!
        TBranch        *b_npart;   //!
        TBranch        *b_thetaP;   //!
        TBranch        *b_phiP;   //!
        TBranch        *b_thetaT;   //!
        TBranch        *b_phiT;   //!
        TBranch        *b_x_part;   //!
        TBranch        *b_y_part;   //!
        TBranch        *b_id_part;   //!
        TBranch        *b_st_part;   //!
        TBranch        *b_pid_part;   //!


        UUana(TTree *tree=0);
        virtual ~UUana();
        virtual Int_t    Cut(Long64_t entry);
        virtual Int_t    GetEntry(Long64_t entry);
        virtual Long64_t LoadTree(Long64_t entry);
        virtual void     Init(TTree *tree);
        virtual void     Loop();
        virtual Bool_t   Notify();
        virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef UUana_cxx
UUana::UUana(TTree *tree) : fChain(0)
{
    // if parameter tree is not specified (or zero), connect the file
    // used to generate this class and read the Tree.
    if (tree == 0) {
        TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("user.phuo.7638948._000002.ampt.root");
        if (!f || !f->IsOpen()) {
            f = new TFile("user.phuo.7638948._000002.ampt.root");
        }
        f->GetObject("tt",tree);
        
    }
    Init(tree);
}

UUana::~UUana()
{
    if (!fChain) return;
    delete fChain->GetCurrentFile();
}

Int_t UUana::GetEntry(Long64_t entry)
{
    // Read contents of entry.
    if (!fChain) return 0;
    return fChain->GetEntry(entry);
}
Long64_t UUana::LoadTree(Long64_t entry)
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

void UUana::Init(TTree *tree)
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

    fChain->SetBranchAddress("ieve", &ieve, &b_ieve);
    fChain->SetBranchAddress("jeve", &jeve, &b_jeve);
    fChain->SetBranchAddress("bimp", &bimp, &b_bimp);
    fChain->SetBranchAddress("eid", eid, &b_eid);
    fChain->SetBranchAddress("ntrk", &ntrk, &b_ntrk);
    fChain->SetBranchAddress("id", id, &b_id);
    fChain->SetBranchAddress("phi", phi, &b_phi);
    fChain->SetBranchAddress("pt", pt, &b_pt);
    fChain->SetBranchAddress("eta", eta, &b_eta);
    fChain->SetBranchAddress("ms", ms, &b_ms);
    fChain->SetBranchAddress("npart", &npart, &b_npart);
    fChain->SetBranchAddress("thetaP", &thetaP, &b_thetaP);
    fChain->SetBranchAddress("phiP", &phiP, &b_phiP);
    fChain->SetBranchAddress("thetaT", &thetaT, &b_thetaT);
    fChain->SetBranchAddress("phiT", &phiT, &b_phiT);
    fChain->SetBranchAddress("x_part", x_part, &b_x_part);
    fChain->SetBranchAddress("y_part", y_part, &b_y_part);
    fChain->SetBranchAddress("id_part", id_part, &b_id_part);
    fChain->SetBranchAddress("st_part", st_part, &b_st_part);
    fChain->SetBranchAddress("pid_part", pid_part, &b_pid_part);
    Notify();
}

Bool_t UUana::Notify()
{
    // The Notify() function is called when a new file is opened. This
    // can be either for a new TTree in a TChain or when when a new TTree
    // is started when using PROOF. It is normally not necessary to make changes
    // to the generated code, but the routine can be extended by the
    // user if needed. The return value is currently not used.

    return kTRUE;
}

void UUana::Show(Long64_t entry)
{
    // Print contents of entry.
    // If entry is not specified, print current entry
    if (!fChain) return;
    fChain->Show(entry);
}
Int_t UUana::Cut(Long64_t entry)
{
    // This function may be called from Loop.
    // returns  1 if entry is accepted.
    // returns -1 otherwise.
    return 1;
}
#endif // #ifdef UUana_cxx
