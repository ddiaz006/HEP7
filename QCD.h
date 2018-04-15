//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Apr 12 22:38:24 2018 by ROOT version 6.08/00
// from TTree Events/created: Wed Apr 11 23:58:23 2018 HepMC 2.06.09
// found on file: qcd_bu.root
//////////////////////////////////////////////////////////

#ifndef QCD_h
#define QCD_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class QCD {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain
   //static const float R; //= 0.8;

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           Event_number;
   Int_t           Event_numberMP;
   Double_t        Event_scale;
   Double_t        Event_alphaQCD;
   Double_t        Event_alphaQED;
   Int_t           Event_barcodeSPV;
   Int_t           Event_numberV;
   Int_t           Event_barcodeBP1;
   Int_t           Event_barcodeBP2;
   Int_t           Event_numberP;
   Double_t        Xsection_value;
   Double_t        Xsection_error;
   Int_t           PDF_parton1;
   Int_t           PDF_parton2;
   Double_t        PDF_x1;
   Double_t        PDF_x2;
   Double_t        PDF_Q2;
   Double_t        PDF_x1f;
   Double_t        PDF_x2f;
   Int_t           PDF_id1;
   Int_t           PDF_id2;
   Double_t        Particle_x[1524];   //[Event_numberP]
   Double_t        Particle_y[1524];   //[Event_numberP]
   Double_t        Particle_z[1524];   //[Event_numberP]
   Double_t        Particle_ctau[1524];   //[Event_numberP]
   Double_t        Particle_barcode[1524];   //[Event_numberP]
   Int_t           Particle_pid[1524];   //[Event_numberP]
   Double_t        Particle_px[1524];   //[Event_numberP]
   Double_t        Particle_py[1524];   //[Event_numberP]
   Double_t        Particle_pz[1524];   //[Event_numberP]
   Double_t        Particle_energy[1524];   //[Event_numberP]
   Double_t        Particle_mass[1524];   //[Event_numberP]
   Int_t           Particle_status[1524];   //[Event_numberP]
   Int_t           Particle_d1[1524];   //[Event_numberP]
   Int_t           Particle_d2[1524];   //[Event_numberP]
   float R = 0.8;

   // List of branches
   TBranch        *b_Event_number;   //!
   TBranch        *b_Event_numberMP;   //!
   TBranch        *b_Event_scale;   //!
   TBranch        *b_Event_alphaQCD;   //!
   TBranch        *b_Event_alphaQED;   //!
   TBranch        *b_Event_barcodeSPV;   //!
   TBranch        *b_Event_numberV;   //!
   TBranch        *b_Event_barcodeBP1;   //!
   TBranch        *b_Event_barcodeBP2;   //!
   TBranch        *b_Event_numberP;   //!
   TBranch        *b_Xsection_value;   //!
   TBranch        *b_Xsection_error;   //!
   TBranch        *b_PDF_parton1;   //!
   TBranch        *b_PDF_parton2;   //!
   TBranch        *b_PDF_x1;   //!
   TBranch        *b_PDF_x2;   //!
   TBranch        *b_PDF_Q2;   //!
   TBranch        *b_PDF_x1f;   //!
   TBranch        *b_PDF_x2f;   //!
   TBranch        *b_PDF_id1;   //!
   TBranch        *b_PDF_id2;   //!
   TBranch        *b_Particle_x;   //!
   TBranch        *b_Particle_y;   //!
   TBranch        *b_Particle_z;   //!
   TBranch        *b_Particle_ctau;   //!
   TBranch        *b_Particle_barcode;   //!
   TBranch        *b_Particle_pid;   //!
   TBranch        *b_Particle_px;   //!
   TBranch        *b_Particle_py;   //!
   TBranch        *b_Particle_pz;   //!
   TBranch        *b_Particle_energy;   //!
   TBranch        *b_Particle_mass;   //!
   TBranch        *b_Particle_status;   //!
   TBranch        *b_Particle_d1;   //!
   TBranch        *b_Particle_d2;   //!

   QCD(TTree *tree=0);
   virtual ~QCD();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   double           delta_R(double Eta, double Phi);
   void             Clustering(vector<int> &P_list, vector<int> &J_list,vector<int> &R_list, vector<double> v_PPt, vector<double> v_PEta, vector<double> v_PPhi);
};

#endif

#ifdef QCD_cxx
QCD::QCD(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("qcd.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("qcd.root");
      }
      f->GetObject("Events",tree);

   }
   Init(tree);
}

QCD::~QCD()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}



Int_t QCD::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t QCD::LoadTree(Long64_t entry)
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

void QCD::Init(TTree *tree)
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

   fChain->SetBranchAddress("Event_number", &Event_number, &b_Event_number);
   fChain->SetBranchAddress("Event_numberMP", &Event_numberMP, &b_Event_numberMP);
   fChain->SetBranchAddress("Event_scale", &Event_scale, &b_Event_scale);
   fChain->SetBranchAddress("Event_alphaQCD", &Event_alphaQCD, &b_Event_alphaQCD);
   fChain->SetBranchAddress("Event_alphaQED", &Event_alphaQED, &b_Event_alphaQED);
   fChain->SetBranchAddress("Event_barcodeSPV", &Event_barcodeSPV, &b_Event_barcodeSPV);
   fChain->SetBranchAddress("Event_numberV", &Event_numberV, &b_Event_numberV);
   fChain->SetBranchAddress("Event_barcodeBP1", &Event_barcodeBP1, &b_Event_barcodeBP1);
   fChain->SetBranchAddress("Event_barcodeBP2", &Event_barcodeBP2, &b_Event_barcodeBP2);
   fChain->SetBranchAddress("Event_numberP", &Event_numberP, &b_Event_numberP);
   fChain->SetBranchAddress("Xsection_value", &Xsection_value, &b_Xsection_value);
   fChain->SetBranchAddress("Xsection_error", &Xsection_error, &b_Xsection_error);
   fChain->SetBranchAddress("PDF_parton1", &PDF_parton1, &b_PDF_parton1);
   fChain->SetBranchAddress("PDF_parton2", &PDF_parton2, &b_PDF_parton2);
   fChain->SetBranchAddress("PDF_x1", &PDF_x1, &b_PDF_x1);
   fChain->SetBranchAddress("PDF_x2", &PDF_x2, &b_PDF_x2);
   fChain->SetBranchAddress("PDF_Q2", &PDF_Q2, &b_PDF_Q2);
   fChain->SetBranchAddress("PDF_x1f", &PDF_x1f, &b_PDF_x1f);
   fChain->SetBranchAddress("PDF_x2f", &PDF_x2f, &b_PDF_x2f);
   fChain->SetBranchAddress("PDF_id1", &PDF_id1, &b_PDF_id1);
   fChain->SetBranchAddress("PDF_id2", &PDF_id2, &b_PDF_id2);
   fChain->SetBranchAddress("Particle_x", Particle_x, &b_Particle_x);
   fChain->SetBranchAddress("Particle_y", Particle_y, &b_Particle_y);
   fChain->SetBranchAddress("Particle_z", Particle_z, &b_Particle_z);
   fChain->SetBranchAddress("Particle_ctau", Particle_ctau, &b_Particle_ctau);
   fChain->SetBranchAddress("Particle_barcode", Particle_barcode, &b_Particle_barcode);
   fChain->SetBranchAddress("Particle_pid", Particle_pid, &b_Particle_pid);
   fChain->SetBranchAddress("Particle_px", Particle_px, &b_Particle_px);
   fChain->SetBranchAddress("Particle_py", Particle_py, &b_Particle_py);
   fChain->SetBranchAddress("Particle_pz", Particle_pz, &b_Particle_pz);
   fChain->SetBranchAddress("Particle_energy", Particle_energy, &b_Particle_energy);
   fChain->SetBranchAddress("Particle_mass", Particle_mass, &b_Particle_mass);
   fChain->SetBranchAddress("Particle_status", Particle_status, &b_Particle_status);
   fChain->SetBranchAddress("Particle_d1", Particle_d1, &b_Particle_d1);
   fChain->SetBranchAddress("Particle_d2", Particle_d2, &b_Particle_d2);
   Notify();
}

Bool_t QCD::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void QCD::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t QCD::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef QCD_cxx
