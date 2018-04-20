#define QCD_cxx
#include "QCD.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TROOT.h"
#include <TChain.h>
#include <vector>
#include <cmath>
#include "TGraph.h"
#include "TStyle.h"
#include "TColor.h"
#include "TF2.h"
#include "TExec.h"

double QCD::delta_R(double dEta, double dPhi)
{
     return sqrt( (dEta*dEta) + (dPhi*dPhi));
}

void QCD::Erase(std::vector<double> &v_PPt, std::vector<double> &v_PEta, std::vector<double> &v_PPhi, std::vector<double> &v_PPx, std::vector<double> &v_PPy, std::vector<double> &v_PPz, std::vector<double> &v_PEn){
   bool clean = false;
   int count = 0;
   while(clean !=true){
   //for(int l  = 0; l<v_PPx.size(); l++){std::cout<<count <<"............."<<v_PPx[l]<<std::endl;}
   for (int i = 0; i< v_PPx.size(); i++){
      bool found = false;
      if(v_PPx[i] >= Null_d){
         v_PPt .erase(v_PPt .begin()  + i); 
         v_PEta.erase(v_PEta.begin()  + i); 
         v_PPhi.erase(v_PPhi.begin()  + i);  
         v_PPx .erase(v_PPx .begin()  + i); 
         v_PPy .erase(v_PPy .begin()  + i); 
         v_PPz .erase(v_PPz .begin()  + i); 
         v_PEn .erase(v_PEn .begin()  + i);
         found = true;
         count = count + 1; 
      }
      if (found == true) break;
      if (i == (v_PPx.size() - 1)) clean =true;
      }
   }
}

void QCD::Clustering( int &dP_listSize, std::vector<int> &dJ_list, std::vector<int> &dR_list, std::vector<double> &dv_PPt, std::vector<double> &dv_PEta, std::vector<double> &dv_PPhi, std::vector<double> &dv_PPx, std::vector<double> &dv_PPy, std::vector<double> &dv_PPz, std::vector<double> &dv_PEn, std::vector<double> &dv_JPt, std::vector<double> &dv_JEta, std::vector<double> &dv_JPhi, std::vector<double> &dv_JPx, std::vector<double> &dv_JPy, std::vector<double> &dv_JPz, std::vector<double> &dv_JEn, int &dcounter)
{
   std::vector<double> d_ij;
   dcounter = 0;
   double dimin = 999999;
   //find di_min
   for(int ii=0; ii<dP_listSize; ii++){if ((1.0/(dv_PPt[ii]*dv_PPt[ii])) <dimin  ) dimin = 1.0/(dv_PPt[ii]*dv_PPt[ii]);}
   for(int i = 0; i<dP_listSize-1; i++)
   {  
      d_ij.clear();
      double dmin = 99999999;
      double id = -1;
      double PtInv2_i = 1.0/(dv_PPt[i]*dv_PPt[i]);
      double PtInv2_min = 9999999;
      for(int j = i+1; j<dP_listSize; j++)
      {
         double PtInv2_j = 1.0/(dv_PPt[j]*dv_PPt[j]);
         //min inverse pt-squared
         if( PtInv2_i< PtInv2_j ) PtInv2_min = PtInv2_i;
         else PtInv2_min = PtInv2_j;
         d_ij.push_back( PtInv2_min*delta_R( (dv_PEta[i]-dv_PEta[j]) , (dv_PPhi[i]-dv_PPhi[j]) )*delta_R( (dv_PEta[i]-dv_PEta[j]) , (dv_PPhi[i]-dv_PPhi[j]) )*(1.0/(R*R)));  
      }
      
      for(int j = i+1; j<dP_listSize; j++)
      {  
         if (d_ij[j] < dmin){ dmin = d_ij[j]; id = j;}
      }
      if(dmin > dimin && dv_PPt[i]>1000.0 && fabs(dv_PEta[i])<1.0) {
         dJ_list.push_back(i);
         //std::cout<<"JetEta:  "<<dv_PEta [i]<<std::endl;
         //put jet into jet std::vectors
         dv_JPhi.push_back(dv_PPhi[i]); dv_PPhi[i] = Null_d;
         dv_JEta.push_back(dv_PEta[i]); dv_PEta[i] = Null_d;
         dv_JPt .push_back(dv_PPt [i]); dv_PPt [i] = Null_d;
         dv_JPx .push_back(dv_PPz [i]); dv_PPx [i] = Null_d;
         dv_JPy .push_back(dv_PPy [i]); dv_PPy [i] = Null_d;
         dv_JPz .push_back(dv_PPz [i]); dv_PPz [i] = Null_d;
         dv_JEn .push_back(dv_PEn [i]); dv_PEn [i] = Null_d;
         Erase(dv_PPt, dv_PEta, dv_PPhi, dv_PPx, dv_PPy, dv_PPz, dv_PEn);
         dP_listSize = dv_PPt.size();
      }    
         else {
         if(i<dP_listSize-1){
         dR_list.push_back(i); dR_list.push_back(id);
         //so add as new entries the merged objects, set the old quantites to a value marked for erasure.
         dv_PPt .push_back(dv_PPt [i]  + dv_PPt [id]); dv_PPt [i] = Null_d; dv_PPt [id] = Null_d;
      //   dv_PPhi.push_back(dv_PPhi[i]  + dv_PPhi[id]); dv_PPhi[i] = Null_d; dv_PPhi[id] = Null_d;  
         dv_PPx .push_back(dv_PPx [i]  + dv_PPx [id]); dv_PPx [i] = Null_d; dv_PPx [id] = Null_d;
         dv_PPy .push_back(dv_PPy [i]  + dv_PPy [id]); dv_PPy [i] = Null_d; dv_PPy [id] = Null_d;
         dv_PPz .push_back(dv_PPz [i]  + dv_PPz [id]); dv_PPz [i] = Null_d; dv_PPz [id] = Null_d;
         dv_PEn .push_back(dv_PEn [i]  + dv_PEn [id]); dv_PEn [i] = Null_d; dv_PEn [id] = Null_d;
      //   dv_PEta.push_back(dv_PEta[i]  + dv_PEta[id]); dv_PEta[i] = Null_d; dv_PEta[id] = Null_d;
         double ddenom = (dv_PEn[dv_PEn.size() -1]*1.00000001) - dv_PPz[dv_PEn.size() -1];
         double dnum  = dv_PEn[dv_PEn.size()-1] + dv_PPz[dv_PEn.size()-1];
         double dPEta = 0.5*log( fabs(dnum/ddenom) );
         double dPPhi = atan2(dv_PPx[dv_PEn.size()-1],dv_PPy[dv_PEn.size()-1]);
         dv_PPhi.push_back(dPPhi); dv_PPhi[i] = Null_d; dv_PPhi[id] = Null_d;  
         dv_PEta.push_back(dPEta); dv_PEta[i] = Null_d; dv_PEta[id] = Null_d;
         dcounter= dcounter + 1;
         Erase(dv_PPt, dv_PEta, dv_PPhi, dv_PPx, dv_PPy, dv_PPz, dv_PEn);
         dP_listSize = dv_PPt.size();
         }
         else{
         //so add as new entries the merged objects, set the old quantites to a value marked for erasure.
         dv_PPt [i] = Null_d;
      //   dv_PPhi[i] = Null_d;  
         dv_PPx [i] = Null_d;
         dv_PPy [i] = Null_d;
         dv_PPz [i] = Null_d;
         dv_PEn [i] = Null_d;
         //dv_PEta[i] = Null_d;
         double ddenom = (dv_PEn[dv_PEn.size() -1]*1.00000001) - dv_PPz[dv_PEn.size() -1];
         double dnum  = dv_PEn[dv_PEn.size()-1] + dv_PPz[dv_PEn.size()-1];
         double dPEta = 0.5*log( fabs(dnum/ddenom) );
         double dPPhi = atan2(dv_PPx[dv_PEn.size()-1],dv_PPy[dv_PEn.size()-1]);
         dv_PPhi.push_back(dPPhi); dv_PPhi[i] = Null_d; dv_PPhi[id] = Null_d;  
         dv_PEta.push_back(dPEta); dv_PEta[i] = Null_d; dv_PEta[id] = Null_d;
         Erase(dv_PPt, dv_PEta, dv_PPhi, dv_PPx, dv_PPy, dv_PPz, dv_PEn);
         dP_listSize = dv_PPt.size();
         }
      }
   }
}

void QCD::Loop()
{
    if (fChain == 0) return;
    Long64_t nentries = fChain->GetEntriesFast();
   
    Long64_t nbytes = 0, nb = 0;

    //Defining histograms to save
    TH1F *h_PPt   = new TH1F("PP_t","PP_t", 100, 0,700);
    TH1F *h_PPhi  = new TH1F("PPhi","PPhi", 30, -1.65,1.65);
    TH1F *h_PEta  = new TH1F("PEta","PEta", 30, -3.2,3.2);
    TH1F *h_JPt = new TH1F("JPt","JPt",30,999.5,2999.5);
    TH1F *h_NJets =new TH1F("NJets","NJets",4,0.5,3.5);
   
    for (Long64_t jentry=0; jentry<nentries/*nentries*/;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
       
     //make arrays to store individ. values for later
    std::vector<double> v_PPhi; std::vector<double> v_JPhi;
    std::vector<double> v_PEta; std::vector<double> v_JEta;
    std::vector<double> v_PPt;  std::vector<double> v_JPt;
    std::vector<double> v_PPx;  std::vector<double> v_JPx;
    std::vector<double> v_PPy;  std::vector<double> v_JPy;
    std::vector<double> v_PPz;  std::vector<double> v_JPz;
    std::vector<double> v_PEn;  std::vector<double> v_JEn;
    std::vector<int>   P_list;  int P_listSize = 0; 
    std::vector<int>   J_list;
    std::vector<int>   R_list;
    std::vector<int>   v_Particle_d1;
    //std::vector<TLorentzVector> P;
    int counter;
    
    for (int z = 0; z<Event_numberP; z++){
      //Particle Pt,phi
      double PPt = sqrt(Particle_px[z]*Particle_px[z] + Particle_py[z]*Particle_py[z]);
      double PPhi = atan2(Particle_px[z],Particle_py[z]);
      double PP = sqrt(Particle_px[z]*Particle_px[z] + Particle_py[z]*Particle_py[z] + Particle_pz[z]*Particle_pz[z]);
      double denom = (Particle_energy[z]*1.00000001) - Particle_pz[z];
      double num  = Particle_energy[z] + Particle_pz[z];
      double PEta = 0.5*log( fabs(num/denom) );
      // double JPt=v_JPt[z]
      //TLorentzVector dummy;
      //dummy.SetPxPyPzE(Particle_px[z], Particle_py[z], Particle_pz[z], Particle_energy[z]);
      //P.push_back(dummy);
      //P[z].Print();
//    std::cout<<z<<":     Eta: "<<PEta <<" Phi: "<<PPhi<<" Px: "<<Particle_px[z]<<" Py: "<<Particle_py[z]<<" Pz: "<<Particle_pz[z]<< " E: "<<Particle_energy[z] <<" x: "<<Particle_x[z]<<" y: "<<Particle_y[z]<<" z: "<<Particle_z[z]<<" Num: "<<num<<" Denom: "<<denom<<std::endl;

      //fill particle histograms and arrays
      if(fabs(PEta) <5 && Particle_d1[z] <0){
        h_PPt->Fill(PPt);   v_PPt.push_back(PPt);
        h_PPhi->Fill(PPhi); v_PPhi.push_back(PPhi);
        h_PEta->Fill(PEta); v_PEta.push_back(PEta);
        v_PPx.push_back(Particle_px[z]); v_PPy.push_back(Particle_py[z]); v_PPz.push_back(Particle_pz[z]); v_PEn.push_back(Particle_energy[z]); 
        //make list of particles
        P_list.push_back(z); P_listSize= P_listSize + 1; 
      }

    }//getting individual data
 //trying a dummy particle to avoid missing one in clustering algortithm
 v_PPx.push_back(0.0); v_PPy.push_back(0.0); v_PPz.push_back(0.0); v_PEn.push_back(0.0); 
 v_PPt.push_back(0.0); v_PEta.push_back(0.0); v_PPhi.push_back(0.0);
 
//    for(int k=0; k<v_PEta.size(); k++){std::cout<<v_PEta[k]<<"   "<<v_PPhi[k]<<"   "<<v_PPt[k]<<std::endl;} //for initial particle dist. 
    std::cout<<"Before Clustering: NParticles: "<<v_PPt.size()<<"        NJets: "<<v_JPt.size()<<"       PList size:   "<<P_listSize<<std::endl;
    while(v_PPt.size()>1){
      Clustering( P_listSize, J_list, R_list, v_PPt, v_PEta, v_PPhi, v_PPx, v_PPy, v_PPz, v_PEn, v_JPt, v_JEta, v_JPhi, v_JPx, v_JPy, v_JPz, v_JEn, counter);
      //for(int k=0; k<v_JEta.size(); k++){std::cout<<v_JEta[k]<<"   "<<v_JPhi[k]<<std::endl;} //for plotting subsequent dist. 
      //std::cout<<"NParticles: "<<v_PPt.size()<<"        NJets: "<<v_JPt.size()<<"       PList size:   "<<P_listSize<<std::endl;
     }
    std::cout<<"NParticles: "<<v_PPt.size()<<"       -NJets: "<<v_JPt.size()<<"       PList size:   "<<P_listSize<<std::endl;
//    for(int k=0; k<v_JEta.size(); k++){std::cout<<v_JEta[k]<<"   "<<v_JPhi[k]<<"   "<<v_JPt[k]<<std::endl;} //for plotting final Jet dist. 
 
    for (int jt=0;jt<v_JPt.size();jt++)
      {
	h_JPt->Fill(v_JPt[jt]);
	h_NJets->Fill(v_JPt.size());
      }
    
    std::cout<<"**************************  Event: "<<jentry<<"   ******************************" <<std::endl;
    }//loop over events
    TFile *outfile = new TFile("histos_QCD.root","RECREATE");
    outfile->cd();
    h_JPt->SetFillColor(40);
    h_JPt->SetLineColor(35);
    h_JPt->GetXaxis()->SetTitle("Events");
    h_JPt->GetYaxis()->SetTitle("Jet Transverse Momentum(GeV)");
    h_JPt->Write();
    h_NJets->SetFillColor(45);
    h_NJets->GetXaxis()->SetTitle("Events");
    h_NJets->GetYaxis()->SetTitle("Number of Jets");
    h_NJets->SetLineColor(36);
    h_NJets->Write();
    h_PEta->Write();
    h_PPhi->Write();
    h_PPt->Write();
    outfile->Close();
}//end Loop()
