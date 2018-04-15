#define QCD_cxx
#include "QCD.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

   /* 
    float d_ij[Event_numberP][Event_numberP];
    for(int i = 0; i<Event_numberP; i++){
      float dR_ij[Event_numberP][Event_numberP];
      for(int j = 0; j <Event_numberP; j++){
         if(i!=j){
           //delta R
           dR_ij[i][j] = sqrt( (a_PPhi[i]-a_PPhi[j])*(a_PPhi[i]-a_PPhi[j]) + (a_PEta[i]-a_PEta[j])*(a_PEta[i]-a_PEta[j]));
         }
         else {dR_ij[i][j] = 9999999;}   
         d_ij[i][j] = Pt_min*(1.0/(R*R))*(dR_ij[i][j]*dR_ij[i][j]);
      }
     }*/
float QCD::delta_R(float dEta, float dPhi)
{
     return sqrt( (dEta*dEta) + (dPhi*dPhi));
}

void QCD::Clustering(vector<int> &P_list, vector<int> &J_list, vector<float> v_PPt, vector<float> v_PEta, vector<float> v_PPhi)
{
   vector<float> d_ij;
   for(int i = 0; i<P_list.size()-1; i++)
   {  
      d_ij.clear();
      float dmin = 99999999;
      float id = -1;
      float PtInv2_i = 1.0/(v_PPt[i]*v_PPt[i]);
      float PtInv2_min = 9999999;
      for(int j = i+1; j<P_list.size(); j++)
      {
         float PtInv2_j = 1.0/(v_PPt[j]*v_PPt[j]);
         //min inverse pt-squared
         if( PtInv2_i< PtInv2_j ) PtInv2_min = PtInv2_i;
         else PtInv2_min = PtInv2_j;
         //cout <<"DeltaR_"<<i<<","<<j<<": "<<delta_R( (v_PEta[i]-v_PEta[j]) , (v_PPhi[i]-v_PPhi[j]) )<<" PEtai: "<<v_PEta[i]<<" PEtaj: "<<v_PEta[j]<<endl;
         d_ij.push_back( PtInv2_min*delta_R( (v_PEta[i]-v_PEta[j]) , (v_PPhi[i]-v_PPhi[j]) )*delta_R( (v_PEta[i]-v_PEta[j]) , (v_PPhi[i]-v_PPhi[j]) )*(1.0/(R*R)));  
      }
    /*  for(int j = i+1; j<P_list.size(); j++)
      {
         if (d_ij[j] < dmin){ dmin = d_ij[j]; id = j;}
      }
      cout << i<<":   dmin: "<<dmin<<"   id: "<<id<<endl; */
   }

}

void QCD::Loop()
{
//      root> .L QCD.C
//      root> QCD t
//      root> t.Loop();       // Loop on all entries
   if (fChain == 0) return;
   Long64_t nentries = fChain->GetEntriesFast();
   
   Long64_t nbytes = 0, nb = 0;

   //Defining histograms to save
   TH1F *h_PPt   = new TH1F("PP_t","PP_t", 100, 0,700);
   TH1F *h_PPhi  = new TH1F("PPhi","PPhi", 30, -1.65,1.65);
   TH1F *h_PEta  = new TH1F("PEta","PEta", 30, -3.2,3.2);
   

   for (Long64_t jentry=0; jentry<1/*nentries*/;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
    
    //make arrays to store individ. values for later
    vector<float> v_PPhi;
    vector<float> v_PEta;
    vector<float> v_PPt;
    vector<int>   P_list;
    vector<int>   J_list;
    vector<TLorentzVector> P;
    
    for (int z = 0; z<Event_numberP; z++){
      //Particle Pt,phi
      float PPt = sqrt(Particle_px[z]*Particle_px[z] + Particle_py[z]*Particle_py[z]);
      float PPhi = atan2(Particle_px[z],Particle_py[z]);
      float PP = sqrt(Particle_px[z]*Particle_px[z] + Particle_py[z]*Particle_py[z] + Particle_pz[z]*Particle_pz[z]);
      double denom = (Particle_energy[z]*1.00000001) - Particle_pz[z];
      double num  = (Particle_energy[z]) + Particle_pz[z];
      float PEta = 0.5*log( num/denom );//fix eta range = -inf/nan cases
      TLorentzVector dummy;
      dummy.SetPxPyPzE(Particle_px[z], Particle_py[z], Particle_pz[z], Particle_energy[z]);
      P.push_back(dummy);
      P[z].Print();
      cout<<z<<":     Eta: "<<PEta <<" Phi: "<<PPhi<<" Px: "<<Particle_px[z]<<" Py: "<<Particle_py[z]<<" Pz: "<<Particle_pz[z]<< " E: "<<Particle_energy[z] <<" x: "<<Particle_x[z]<<" y: "<<Particle_y[z]<<" z: "<<Particle_z[z]<<" Num: "<<num<<" Denom: "<<denom<<endl;

      //fill particle histograms and arrays
      h_PPt->Fill(PPt);   v_PPt.push_back(PPt);
      h_PPhi->Fill(PPhi); v_PPhi.push_back(PPhi);
      h_PEta->Fill(PEta); v_PEta.push_back(PEta); 
      //make list of particles
      P_list.push_back(z);
    }//getting individual data
    
      //if(jentry==0)//cout<<P_list.size()<<"    "<<v_PPt.size()<<"   "<<v_PPhi.size()<<endl;
      Clustering( P_list, J_list, v_PPt, v_PEta, v_PPhi);

   /* 
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
     }*/
   }//loop over events

TFile *outfile = new TFile("histos_QCD.root","RECREATE");
outfile->cd();
h_PPt->Write();
h_PPhi->Write();
h_PEta->Write();
outfile->Close();
}//end Loop()
