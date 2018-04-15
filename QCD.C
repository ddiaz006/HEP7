#define QCD_cxx
#include "QCD.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

double QCD::delta_R(double dEta, double dPhi)
{
     return sqrt( (dEta*dEta) + (dPhi*dPhi));
}

void QCD::Clustering(vector<int> &P_list, vector<int> &J_list, vector<int> &R_list, vector<double> v_PPt, vector<double> v_PEta, vector<double> v_PPhi)
{
   vector<double> d_ij;
   for(int i = 0; i<P_list.size()-1; i++)
   {  
      d_ij.clear();
      double dmin = 99999999;
      double id = -1;
      double PtInv2_i = 1.0/(v_PPt[i]*v_PPt[i]);
      double PtInv2_min = 9999999;
      for(int j = i+1; j<P_list.size(); j++)
      {
         double PtInv2_j = 1.0/(v_PPt[j]*v_PPt[j]);
         //min inverse pt-squared
         if( PtInv2_i< PtInv2_j ) PtInv2_min = PtInv2_i;
         else PtInv2_min = PtInv2_j;
         //cout <<"DeltaR_"<<i<<","<<j<<": "<<delta_R( (v_PEta[i]-v_PEta[j]) , (v_PPhi[i]-v_PPhi[j]) )<<" PEtai: "<<v_PEta[i]<<" PEtaj: "<<v_PEta[j]<<endl;
         d_ij.push_back( PtInv2_min*delta_R( (v_PEta[i]-v_PEta[j]) , (v_PPhi[i]-v_PPhi[j]) )*delta_R( (v_PEta[i]-v_PEta[j]) , (v_PPhi[i]-v_PPhi[j]) )*(1.0/(R*R)));  
      }
      for(int j = i+1; j<P_list.size(); j++)
      {  //cout <<"****j=:"<<j<<"*** d_ij[j]: "<<d_ij[j]<<endl;
         if (d_ij[j] < dmin){ dmin = d_ij[j]; id = j;}
      }
      //cout << i<<":   dmin: "<<dmin<<"   id: "<<id<<endl;

      if(dmin < PtInv2_i) {R_list.push_back(i); R_list.push_back(id);}
      else J_list.push_back(i);
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
    vector<double> v_PPhi; vector<double> v_JPhi;
    vector<double> v_PEta; vector<double> v_JEta;
    vector<double> v_PPt;  vector<double> v_JPt;
    vector<double> v_PPx;  vector<double> v_JPx;
    vector<double> v_PPy;  vector<double> v_JPy;
    vector<double> v_PPz;  vector<double> v_JPz;
    vector<double> v_PEn;  vector<double> v_JEn;
    vector<int>   P_list;
    vector<int>   J_list;
    vector<int>   R_list;
    vector<TLorentzVector> P;
    
    for (int z = 0; z<Event_numberP; z++){
      //Particle Pt,phi
      double PPt = sqrt(Particle_px[z]*Particle_px[z] + Particle_py[z]*Particle_py[z]);
      double PPhi = atan2(Particle_px[z],Particle_py[z]);
      double PP = sqrt(Particle_px[z]*Particle_px[z] + Particle_py[z]*Particle_py[z] + Particle_pz[z]*Particle_pz[z]);
      double denom = (Particle_energy[z]*1.00000001) - Particle_pz[z];
      double num  = Particle_energy[z] + Particle_pz[z];
      double PEta = 0.5*log( abs(num/denom) );
      TLorentzVector dummy;
      dummy.SetPxPyPzE(Particle_px[z], Particle_py[z], Particle_pz[z], Particle_energy[z]);
      P.push_back(dummy);
      //P[z].Print();
  //    cout<<z<<":     Eta: "<<PEta <<" Phi: "<<PPhi<<" Px: "<<Particle_px[z]<<" Py: "<<Particle_py[z]<<" Pz: "<<Particle_pz[z]<< " E: "<<Particle_energy[z] <<" x: "<<Particle_x[z]<<" y: "<<Particle_y[z]<<" z: "<<Particle_z[z]<<" Num: "<<num<<" Denom: "<<denom<<endl;

      //fill particle histograms and arrays
      h_PPt->Fill(PPt);   v_PPt.push_back(PPt);
      h_PPhi->Fill(PPhi); v_PPhi.push_back(PPhi);
      h_PEta->Fill(PEta); v_PEta.push_back(PEta);
      v_PPx.push_back(Particle_px[z]); v_PPy.push_back(Particle_py[z]); v_PPz.push_back(Particle_pz[z]); v_PEn.push_back(Particle_energy[z]); 
      //make list of particles
      P_list.push_back(z);
    }//getting individual data
    
      //if(jentry==0)//cout<<P_list.size()<<"    "<<v_PPt.size()<<"   "<<v_PPhi.size()<<endl;
      Clustering( P_list, J_list, R_list, v_PPt, v_PEta, v_PPhi);
      for (int k = 0; k < J_list.size(); k ++)
      {//addin jet info to jet vectors
        v_JPhi.push_back(v_PPhi[J_list[k]] );
        v_JEta.push_back(v_PEta[J_list[k]] );
        v_JPt .push_back(v_PPt [J_list[k]] );
        v_JPx .push_back(v_PPz [J_list[k]] );
        v_JPy .push_back(v_PPy [J_list[k]] );
        v_JPz .push_back(v_PPz [J_list[k]] );
        v_JEn .push_back(v_PEn [J_list[k]] );
      }
      for(int k = 0; k < R_list.size()-1; k= k+2){
        cout<<k<<"         "<<R_list[k]<<endl;
       //merge i,j pairs
       v_PPhi.push_back(v_PPhi[k]+v_PPhi[k+1]);
       v_PEta.push_back(v_PEta[k]+v_PEta[k+1]);
       v_PPt .push_back(v_PPt [k]+v_PPt [k+1]);
       v_PPx .push_back(v_PPx [k]+v_PPx [k+1]);
       v_PPy .push_back(v_PPy [k]+v_PPy [k+1]);
       v_PPz .push_back(v_PPz [k]+v_PPz [k+1]);
       v_PEn .push_back(v_PEn [k]+v_PEn [k+1]);
       //remove merged pairs
       //P_list.erase(P_list.begin());
       //P_list.erase(P_list.begin()+1);
       //v_PPhi  .erase(v_PPhi  .begin()+k);
       //v_PPhi  .erase(v_PPhi  .begin()+(k+1));
       //v_PEta  .erase(v_PEta  .begin()+k);
       //v_PEta  .erase(v_PEta  .begin()+(k+1));
       //v_PPt   .erase(v_PPt   .begin()+k);
       //v_PPt   .erase(v_PPt   .begin()+(k+1));
       //v_PPx   .erase(v_PPx   .begin()+k);
       //v_PPx   .erase(v_PPx   .begin()+(k+1));
       //v_PPy   .erase(v_PPy   .begin()+k);
       //v_PPy   .erase(v_PPy   .begin()+(k+1));
       //v_PPz   .erase(v_PPz   .begin()+k);
       //v_PPz   .erase(v_PPz   .begin()+(k+1));
       //v_PEn   .erase(v_PEn   .begin()+k);
       //v_PEn   .erase(v_PEn   .begin()+(k+1));
      }
      //Clustering( P_list, J_list, R_list, v_PPt, v_PEta, v_PPhi);
  /*    for (int k = 0; k < J_list.size(); k ++)
      {//addin jet info to jet vectors
        v_JPhi.push_back(v_PPhi[J_list[k]] );
        v_JEta.push_back(v_PEta[J_list[k]] );
        v_JPt .push_back(v_PPt [J_list[k]] );
        v_JPx .push_back(v_PPz [J_list[k]] );
        v_JPy .push_back(v_PPy [J_list[k]] );
        v_JPz .push_back(v_PPz [J_list[k]] );
        v_JEn .push_back(v_PEn [J_list[k]] );
      }
      for(int k = 0; k < R_list.size(); k= k+2){
        cout<<k<<"         "<<R_list[k]<<endl;
       //merge i,j pairs
       v_PPhi.push_back(v_PPhi[k]+v_PPhi[k+1]);
       v_PEta.push_back(v_PEta[k]+v_PEta[k+1]);
       v_PPt .push_back(v_PPt [k]+v_PPt [k+1]);
       v_PPx .push_back(v_PPx [k]+v_PPx [k+1]);
       v_PPy .push_back(v_PPy [k]+v_PPy [k+1]);
       v_PPz .push_back(v_PPz [k]+v_PPz [k+1]);
       v_PEn .push_back(v_PEn [k]+v_PEn [k+1]);
       //remove merged pairs
       P_list.erase(P_list.begin()+k);
       P_list.erase(P_list.begin()+(k+1));
       v_PPhi  .erase(v_PPhi  .begin()+k);
       v_PPhi  .erase(v_PPhi  .begin()+(k+1));
       v_PEta  .erase(v_PEta  .begin()+k);
       v_PEta  .erase(v_PEta  .begin()+(k+1));
       v_PPt   .erase(v_PPt   .begin()+k);
       v_PPt   .erase(v_PPt   .begin()+(k+1));
       v_PPx   .erase(v_PPx   .begin()+k);
       v_PPx   .erase(v_PPx   .begin()+(k+1));
       v_PPy   .erase(v_PPy   .begin()+k);
       v_PPy   .erase(v_PPy   .begin()+(k+1));
       v_PPz   .erase(v_PPz   .begin()+k);
       v_PPz   .erase(v_PPz   .begin()+(k+1));
       v_PEn   .erase(v_PEn   .begin()+k);
       v_PEn   .erase(v_PEn   .begin()+(k+1));
      }
*/
   }//loop over events

TFile *outfile = new TFile("histos_QCD.root","RECREATE");
outfile->cd();
h_PPt->Write();
h_PPhi->Write();
h_PEta->Write();
outfile->Close();
}//end Loop()
