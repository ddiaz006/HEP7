#define QCD_cxx
#include "QCD.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

double QCD::delta_R(double dEta, double dPhi)
{
     return sqrt( (dEta*dEta) + (dPhi*dPhi));
}

void QCD::Erase(vector<double> &v_PPt, vector<double> &v_PEta, vector<double> &v_PPhi, vector<double> &v_PPx, vector<double> &v_PPy, vector<double> &v_PPz, vector<double> &v_PEn){
   bool clean = false;
   while(clean !=true){
   for (int i = 0; i< v_PPx.size(); i++){
      bool found = false;
      if(v_PPx[i] == Null_d){
         v_PPt .erase(v_PPt .begin()  + i); 
         v_PEta.erase(v_PEta.begin()  + i); 
         v_PPhi.erase(v_PPhi.begin()  + i);  
         v_PPx .erase(v_PPx .begin()  + i); 
         v_PPy .erase(v_PPy .begin()  + i); 
         v_PPz .erase(v_PPz .begin()  + i); 
         v_PEn .erase(v_PEn .begin()  + i);
         found = true; 
      }
      if (found == true) break;
      if (i == (v_PPx.size() - 1)) clean =true;
      }
   }
}

void QCD::Clustering( int &dP_listSize, vector<int> &dJ_list, vector<int> &dR_list, vector<double> &dv_PPt, vector<double> &dv_PEta, vector<double> &dv_PPhi, vector<double> &dv_PPx, vector<double> &dv_PPy, vector<double> &dv_PPz, vector<double> &dv_PEn, vector<double> &dv_JPt, vector<double> &dv_JEta, vector<double> &dv_JPhi, vector<double> &dv_JPx, vector<double> &dv_JPy, vector<double> &dv_JPz, vector<double> &dv_JEn, int &dcounter)
{
   vector<double> d_ij;
   dcounter = 0;
   double dimin = 999999;
   //find di_min
   for(int ii=0; ii<dP_listSize-1; ii++){if ((1.0/(dv_PPt[ii]*dv_PPt[ii])) <dimin  ) dimin = 1.0/(dv_PPt[ii]*dv_PPt[ii]);}
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
         //cout <<"DeltaR_"<<i<<","<<j<<": "<<delta_R( (v_PEta[i]-v_PEta[j]) , (v_PPhi[i]-v_PPhi[j]) )<<" PEtai: "<<v_PEta[i]<<" PEtaj: "<<v_PEta[j]<<endl;
         d_ij.push_back( PtInv2_min*delta_R( (dv_PEta[i]-dv_PEta[j]) , (dv_PPhi[i]-dv_PPhi[j]) )*delta_R( (dv_PEta[i]-dv_PEta[j]) , (dv_PPhi[i]-dv_PPhi[j]) )*(1.0/(R*R)));  
      }
      for(int j = i+1; j<dP_listSize; j++)
      {  //cout <<"****j=:"<<j<<"*** d_ij[j]: "<<d_ij[j]<<endl;
         if (d_ij[j] < dmin){ dmin = d_ij[j]; id = j;}
      }
      cout << i<<":   dmin: "<<dmin<<"   dimin: "<<dimin<<endl;
      //plan to add stuff not marked for removal to dummy vecotrs then at end delete old vectors and save dummy vectors in new ones 
      if(dmin < dimin) {
         dR_list.push_back(i); dR_list.push_back(id);
         //so add as new entries the merged objects, set the old quantites to a value marked for erasure.
         dv_PPt .push_back(dv_PPt [i]  + dv_PPt [id]); dv_PPt [i] = Null_d; dv_PPt [id] = Null_d;
         dv_PEta.push_back(dv_PEta[i]  + dv_PEta[id]); dv_PEta[i] = Null_d; dv_PEta[id] = Null_d;
         dv_PPhi.push_back(dv_PPhi[i]  + dv_PPhi[id]); dv_PPhi[i] = Null_d; dv_PPhi[id] = Null_d;  
         dv_PPx .push_back(dv_PPx [i]  + dv_PPx [id]); dv_PPx [i] = Null_d; dv_PPx [id] = Null_d;
         dv_PPy .push_back(dv_PPy [i]  + dv_PPy [id]); dv_PPy [i] = Null_d; dv_PPy [id] = Null_d;
         dv_PPz .push_back(dv_PPz [i]  + dv_PPz [id]); dv_PPz [i] = Null_d; dv_PPz [id] = Null_d;
         dv_PEn .push_back(dv_PEn [i]  + dv_PEn [id]); dv_PEn [i] = Null_d; dv_PEn [id] = Null_d;
         dcounter= dcounter + 1;
      }
      else {
         dJ_list.push_back(i);
         //put jet into jet vectors
         dv_JPhi.push_back(dv_PPhi[i]); dv_PPhi[i] = Null_d;
         dv_JEta.push_back(dv_PEta[i]); dv_PEta[i] = Null_d;
         dv_JPt .push_back(dv_PPt [i]); dv_PPt [i] = Null_d;
         dv_JPx .push_back(dv_PPz [i]); dv_PPx [i] = Null_d;
         dv_JPy .push_back(dv_PPy [i]); dv_PPy [i] = Null_d;
         dv_JPz .push_back(dv_PPz [i]); dv_PPz [i] = Null_d;
         dv_JEn .push_back(dv_PEn [i]); dv_PEn [i] = Null_d;
       
      }    
   }
//      for(int zz = 0; zz <dcounter; zz++){dP_list.pop_back();} 
//      for(int k= 0; k < dv_PPt.size(); k++){cout <<"before erase("<<k<<"): "<<dv_PPt[k]<<endl;}
      Erase(dv_PPt, dv_PEta, dv_PPhi, dv_PPx, dv_PPy, dv_PPz, dv_PEn);
      dP_listSize = dv_PPt.size();
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
    vector<int>   P_list;  int P_listSize = 0; 
    vector<int>   J_list;
    vector<int>   R_list;
    vector<int>   v_Particle_d1;
    //vector<TLorentzVector> P;
    int counter;
    
    for (int z = 0; z<Event_numberP; z++){
      //Particle Pt,phi
      double PPt = sqrt(Particle_px[z]*Particle_px[z] + Particle_py[z]*Particle_py[z]);
      double PPhi = atan2(Particle_px[z],Particle_py[z]);
      double PP = sqrt(Particle_px[z]*Particle_px[z] + Particle_py[z]*Particle_py[z] + Particle_pz[z]*Particle_pz[z]);
      double denom = (Particle_energy[z]*1.00000001) - Particle_pz[z];
      double num  = Particle_energy[z] + Particle_pz[z];
      double PEta = 0.5*log( abs(num/denom) );
      //TLorentzVector dummy;
      //dummy.SetPxPyPzE(Particle_px[z], Particle_py[z], Particle_pz[z], Particle_energy[z]);
      //P.push_back(dummy);
      //P[z].Print();
  //    cout<<z<<":     Eta: "<<PEta <<" Phi: "<<PPhi<<" Px: "<<Particle_px[z]<<" Py: "<<Particle_py[z]<<" Pz: "<<Particle_pz[z]<< " E: "<<Particle_energy[z] <<" x: "<<Particle_x[z]<<" y: "<<Particle_y[z]<<" z: "<<Particle_z[z]<<" Num: "<<num<<" Denom: "<<denom<<endl;

      //h_PPt->Fill(PPt);
      ////h_PPhi->Fill(PPhi);
      ///h_PEta->Fill(PEta);
      //fill particle histograms and arrays
      if(abs(PEta) <5 && Particle_d1[z] <0){
        h_PPt->Fill(PPt);   v_PPt.push_back(PPt);
        h_PPhi->Fill(PPhi); v_PPhi.push_back(PPhi);
        h_PEta->Fill(PEta); v_PEta.push_back(PEta);
        v_PPx.push_back(Particle_px[z]); v_PPy.push_back(Particle_py[z]); v_PPz.push_back(Particle_pz[z]); v_PEn.push_back(Particle_energy[z]); 
        //make list of particles
        P_list.push_back(z); P_listSize= P_listSize + 1;
      }

    }//getting individual data
    
      //if(jentry==0)//cout<<P_listSize<<"    "<<v_PPt.size()<<"   "<<v_PPhi.size()<<endl;
      cout<<"NParticles: "<<v_PPt.size()<<"        NJets: "<<v_JPt.size()<<"       PList size:   "<<P_listSize<<endl;
//      for(int k= 0; k < v_PPt.size(); k++){cout <<"intial("<<k<<"): "<<v_PPt[k]<<endl;}
      Clustering( P_listSize, J_list, R_list, v_PPt, v_PEta, v_PPhi, v_PPx, v_PPy, v_PPz, v_PEn, v_JPt, v_JEta, v_JPhi, v_JPx, v_JPy, v_JPz, v_JEn, counter);
      cout<<"NParticles: "<<v_PPt.size()<<"        NJets: "<<v_JPt.size()<<"       PList size:   "<<P_listSize<<endl;
      Clustering( P_listSize, J_list, R_list, v_PPt, v_PEta, v_PPhi, v_PPx, v_PPy, v_PPz, v_PEn, v_JPt, v_JEta, v_JPhi, v_JPx, v_JPy, v_JPz, v_JEn, counter);
      cout<<"NParticles: "<<v_PPt.size()<<"        NJets: "<<v_JPt.size()<<"       PList size:   "<<P_listSize<<endl;
//      for(int k= 0; k < v_PPt.size(); k++){cout <<"final("<<k<<"): "<<v_PPt[k]<<endl;}
   }//loop over events

TFile *outfile = new TFile("histos_QCD.root","RECREATE");
outfile->cd();
h_PPt->Write();
h_PPhi->Write();
h_PEta->Write();
outfile->Close();
}//end Loop()
