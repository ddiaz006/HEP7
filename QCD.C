#define QCD_cxx
#include "QCD.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void QCD::Loop()
{
//      root> .L QCD.C
//      root> QCD t
//      root> t.Loop();       // Loop on all entries
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;
   TH1F *h_Pt  = new TH1F("P_t","P_t", 100, 0,500);
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
   
    for (int i = 0; i<NP; i++){

      //cout <<"PT = "<<sqrt(Particle_px[i]*Particle_px[i] + Particle_py[i]*Particle_py[i])<<endl;
      h_Pt->Fill(sqrt(Particle_px[i]*Particle_px[i] + Particle_py[i]*Particle_py[i]));
    }
   }

TFile *outfile = new TFile("histos_QCD.root","RECREATE");
outfile->cd();
h_Pt->Write();
outfile->Close();
}//end Loop()
