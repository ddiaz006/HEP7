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
   const float R = 0.8;
   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;

   //Defining histograms to save
   TH1F *h_PPt   = new TH1F("PP_t","PP_t", 100, 0,700);
   TH1F *h_PPhi  = new TH1F("PPhi","PPhi", 30, -1.65,1.65);
   TH1F *h_PEta  = new TH1F("PEta","PEta", 30, -3.2,3.2);


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    //make arrays to store individ. values for later
    float a_PPhi[Event_numberP];
    float a_PEta[Event_numberP];
    float a_PPt[Event_numberP];

    for (int z = 0; z<Event_numberP; z++){
      //Particle Pt,phi
      float PPt = sqrt(Particle_px[z]*Particle_px[z] + Particle_py[z]*Particle_py[z]);
      float PPhi = atan(Particle_px[z]/Particle_py[z]);
      float PEta = 0.5*log( (Particle_energy[z] + Particle_pz[z])/(Particle_energy[z] - Particle_pz[z]) );

      //fill particle histograms and arrays
      h_PPt->Fill(PPt);   a_PPt[z] =PPt;
      h_PPhi->Fill(PPhi); a_PPhi[z]=PPhi;
      h_PEta->Fill(PEta); a_PEta[z]=PEta; 

    }//getting individual data
    
    
    float d_ij[Event_numberP][Event_numberP];
    for(int i = 0; i<Event_numberP; i++){
      float PtInv_i = 1.0/(a_PPt[i]*a_PPt[i]);
      float Pt_min = 9999999;
      float dR_ij[Event_numberP][Event_numberP];
      for(int j = 0; j <Event_numberP; j++){
      float PtInv_j = 1.0/(a_PPt[j]*a_PPt[j]);
         if(i!=j){
           //delta R
           dR_ij[i][j] = sqrt( (a_PPhi[i]-a_PPhi[j])*(a_PPhi[i]-a_PPhi[j]) + (a_PEta[i]-a_PEta[j])*(a_PEta[i]-a_PEta[j]));
           //min inverse pt-squared
           if( PtInv_i< PtInv_j ) Pt_min = PtInv_i;
           else Pt_min = PtInv_j;
         }
         else {dR_ij[i][j] = 9999999;}   
         d_ij[i][j] = Pt_min*(1.0/(R*R))*(dR_ij[i][j]*dR_ij[i][j]);
      }
     }
   }//loop over events

TFile *outfile = new TFile("histos_QCD.root","RECREATE");
outfile->cd();
h_PPt->Write();
h_PPhi->Write();
h_PEta->Write();
outfile->Close();
}//end Loop()
