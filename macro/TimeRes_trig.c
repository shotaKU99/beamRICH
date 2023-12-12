#include <math.h>

#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TTree.h"
#include "TCanvas.h"

void TimeRes_trig(int run){

  const Int_t NRGBs = 5;
  const Int_t NCont = 255;
  Double_t stops[NRGBs] = { 0.00, 0.34, 0.61, 0.84, 1.00 };
  Double_t Red[NRGBs] = { 0.00, 0.00, 0.87, 1.00, 0.51 };
  Double_t Green[NRGBs] = { 0.00, 0.81, 1.00, 0.20, 0.00 };
  Double_t Blue[NRGBs] = { 0.51, 1.00, 0.12, 0.00, 0.00 };
  TColor::CreateGradientColorTable(NRGBs, stops, Red, Green, Blue, NCont);
  gStyle->SetNumberContours(NCont);
  //  gStyle->SetOptStat(0);
  //gStyle->SetPadGridX(true);
  //gStyle->SetPadGridY(true);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.15);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetTextFont(132);
  gStyle->SetTitleFont(132,"XYZ");
  gStyle->SetLabelFont(132,"XYZ");
  gStyle->SetStatH(6);
  gStyle->SetStatW(0.30);
  gStyle->SetStatFontSize(0.08);
  gStyle->SetLabelSize(0.05,"XY");
  gStyle->SetTitleSize(0.05,"XY");
  gStyle->SetTitleOffset(0.7,"X");
  gStyle->SetTitleOffset(1.3,"X");
  gStyle->SetTitleSize(0.05);


  // Ch map
  //std::vector<int> ChMapL(12);
  std::vector<int> ChMapX(16);
  ChMapX[0]=0; //X1
  ChMapX[1]=1; //X2
  ChMapX[2]=2; //X3
  ChMapX[3]=3; //X4
  ChMapX[4]=4; //X5
  ChMapX[5]=5; //X6
  ChMapX[6]=6; //X7
  ChMapX[7]=7; //X8
  ChMapX[8]=8; //X9
  ChMapX[9]=9; //X10
  ChMapX[10]=10; //X11
  ChMapX[11]=11; //X12
  ChMapX[12]=12; //X13
  ChMapX[13]=13; //X14
  ChMapX[14]=14; //X15
  ChMapX[15]=15; //X16

  //std::vector<int> ChMapR(12)
  std::vector<int> ChMapY(16);
  ChMapY[0]=22; //Y1
  ChMapY[1]=23; //Y2
  ChMapY[2]=24; //Y3
  ChMapY[3]=25; //Y4
  ChMapY[4]=28; //Y5
  ChMapY[5]=29; //Y6
  ChMapY[6]=30; //Y7
  ChMapY[7]=31; //Y8
  ChMapY[8]=40; //Y9
  ChMapY[9]=41; //Y10
  ChMapY[10]=42; //Y11
  ChMapY[11]=43; //Y12
  ChMapY[12]=44; //Y13
  ChMapY[13]=45; //Y14
  ChMapY[14]=46; //Y15
  ChMapY[15]=47; //Y16
  


//tigger ch
  int trig1ch = 33;
  int trig2ch = 37;
  int trig3ch = 20;
  //int trig4ch = 41;

  //tdc cut
  double tdcrangelowL = 630;
  double tdcrangehighL = 690;
  double tdcrangelowR = 630;
  double tdcrangehighR = 690;

  
  TFile *fin = new TFile(Form("/mnt/c/Users/shota/ELPHroot/data/run%04d.root",run));
  cout << Form("root/run00%2d",run) << " was opened" << endl;
  TTree *tree = (TTree*)fin->Get("tree");

  std::vector<std::vector<double>> *ltdc = 0;
  std::vector<std::vector<double>> *ttdc = 0;


    
  TBranch *Bltdc = 0;
  TBranch *Bttdc = 0;



  tree->SetBranchAddress("ltdc", &ltdc, &Bltdc );
  tree->SetBranchAddress("ttdc", &ttdc, &Bttdc );


  //definitions of histograms
  TH1F *htof_21 = new TH1F("tof between 1 and 2", "tof between 1 and 2", 100,-4,4);
  TH1F *htof_31 = new TH1F("tof between 1 and 3", "tof between 1 and 3", 100,27,35);
  TH1F *htof_32 = new TH1F("tof between 2 and 3", "tof between 2 and 3", 100,27,35);

  TH2F *hThrough_23 = new TH2F("scatterplot for throughing correction of tof_32", "; TOT1; TOF_32",100,-5,20,100,20,40);


  TH2F *hThrough_12_1 = new TH2F("scatterplot for throughing correction of tof_21", "; TOT1; TOF_21",100,0,20,100,-4,4);
  TH2F *hThrough_12_2 = new TH2F("throughing correction of tof_21","; TOT2; TOF_21",100,0,20,100,-4,4);
  TH2F *hThrough_12_1comp_1 = new TH2F("12 after a compensation of TOT1","after calib of TOT1;TOT1;TOF_12_comp",100,0,20,100,-4,4);
  TH2F *hThrough_12_1comp_2 = new TH2F("mpensation of TOT1","after calib of TOT1;TOT2;TOF_12_comp",100,0,20,100,-4,4);
  TH2F *hThrough_12_2comp_1 = new TH2F(" TOT1 & TOT2",";TOT1;TOF_12_comps",100,0,20,100,-4,4);
  TH2F *hThrough_12_2comp_2 = new TH2F("12 after compensations of TOT1 & TOT2",";TOT2;TOF_12_comps",100,0,20,100,-4,4);


  TH2F *hThrough_13_1 = new TH2F("hThrough131", "; TOT1; TOF_31",100,0,25,100,26,34);
  TH2F *hThrough_13_2 = new TH2F("hThrough132", "; TOT3; TOF_31",100,20,60,100,26,34);
  TH2F *hThrough_13_1comp_1 = new TH2F("hThrough131comp1","after calib of TOT1;TOT1;TOF_13_comp",100,0,20,100,-4,4);
  TH2F *hThrough_13_1comp_2 = new TH2F("hThrough131comp2","after calib of TOT1;TOT3;TOF_13_comp",100,20,60,100,-4,4);
  TH2F *hThrough_13_2comp_1 = new TH2F("hThrough132comp1",";TOT1;TOF_13_comps",100,0,20,100,-4,4);
  TH2F *hThrough_13_2comp_2 = new TH2F("hThrough132comp2",";TOT3;TOF_13_comps",100,20,60,100,-4,4);

  TH2F *hThrough_23_1 = new TH2F("hThrough231", "; TOT2; TOF_32",100,0,35,100,26,34);
  TH2F *hThrough_23_2 = new TH2F("hThrough232", "; TOT3; TOF_32",100,20,60,100,26,34);
  TH2F *hThrough_23_1comp_1 = new TH2F("hThrough231comp1","after caib of TOT2;TOT1;TOF_23_comp",100,0,35,100,-4,4);
  TH2F *hThrough_23_1comp_2 = new TH2F("hThrough231comp2","after calb of TOT3;TOT3;TOF_23_comp",100,20,60,100,-4,4);
  TH2F *hThrough_23_2comp_1 = new TH2F("hThrough232comp1","after cali of TOT2;TOT2;TOF_23_comps",100,0,35,100,-4,4);
  TH2F *hThrough_23_2comp_2 = new TH2F("hThrough232comp2","after compensations of TOT3 &2;TOT3;TOF_23_comps",100,20,60,100,-4,4);
 
  //definitions of variables
  double ttdc1,ttdc2,ttdc3;
  double tot1,tot2,tot3;
  double tof_21,tof_31,tof_32;

  double ltdc1, ltdc2, ltdc3;
  
  //definitions of fitting functions
  TF1 *f1 = new TF1("f1","[0] + [1]/sqrt([2]*(x-[3]))",0.000001,15);
  f1->SetParameters(2.12, -4, 0.46, -1.32);
  TF1 *f2 = new TF1("f2","[0] + [1]/sqrt([2]*(x-[3]))",0.000001,15);
  f2->SetParameters(-1.4, 3.6, 0.5, -0.61);

  TF1 *f3 = new TF1("f3","[0] + [1]/sqrt([2]*(x-[3]))",0.000001,15);
  f3->SetParameters(34, -11, 3.5, -1.29);
  //TF1 *f4 = new TF1("f4","[0] + [1]/sqrt([2]*(x-[3]))",20,60);
  //f4->SetParameters(-73,200,0,-1000);
  TF1 *f4 = new TF1("f4","[0] + [1]*x",20,60);
  f4->SetParameters(1.5, -0.039);
  TF1 *f5 = new TF1("f5","[0] + [1]/sqrt([2]*(x-[3]))",0.000001,15);
  f5->SetParameters(33, -7, 2.7, -0.32);
  //TF1 *f6 = new TF1("f6","[0] + [1]/sqrt([2]*(x-[3]))",20,60);
  //f6->SetParameters(-19,41,0,-200);
  TF1 *f6 = new TF1("f6","[0] + [1]*x",20,60);
  f6->SetParameters(1.5, -0.036);


  TF1 *gaus12 = new TF1("gaus12","gaus(0)",-2,2);
  TF1 *gaus13 = new TF1("gaus13","gaus(0)",-2,2);
  TF1 *gaus23 = new TF1("gaus23","gaus(0)",-2,2);

  const int n = tree->GetEntries();





  // Fill
  for (int i=0; i<n; i++){
    if(i%10000 == 0) cout << "event num = " << i << endl;
    tree->GetEntry(i);
    
    std::vector<int> HitChX;
    std::vector<int> HitChY;

    if( ltdc->at(trig1ch).size()<1 || ttdc->at(trig1ch).size()<1 ) continue;
    if( ltdc->at(trig2ch).size()<1 || ttdc->at(trig2ch).size()<1 ) continue;
    if( ltdc->at(trig3ch).size()<1 || ttdc->at(trig3ch).size()<1 ) continue;


    for(int j=0; j<ChMapX.size(); j++){
      if( ltdc->at(ChMapX[j]).size() < 1 || ttdc->at(ChMapX[j]).size() < 1 ) continue;

      HitChX.push_back(j);

    }

    for(int j=0; j<ChMapY.size(); j++){
      if( ltdc->at(ChMapY[j]).size() < 1 || ttdc->at(ChMapY[j]).size() < 1 ) continue;
     
      HitChY.push_back(j);
    }
    
    
    if (HitChX.size() == 1 && HitChY.size() == 1){
      ltdc1 = ltdc->at(trig1ch).at(0);
      ltdc2 = ltdc->at(trig2ch).at(0);
      ltdc3 = ltdc->at(trig3ch).at(0);

      ttdc1 = ttdc->at(trig1ch).at(0);
      ttdc2 = ttdc->at(trig2ch).at(0);
      ttdc3 = ttdc->at(trig3ch).at(0);
      
      tof_21 = ltdc1-ltdc2;
      tof_31 = ltdc1-ltdc3;
      tof_32 = ltdc2-ltdc3;

      tot1 = ltdc1-ttdc1;
      tot2 = ltdc2-ttdc2;
      tot3 = ltdc3-ttdc3;

      htof_21->Fill(tof_21);
      htof_31->Fill(tof_31);
      htof_32->Fill(tof_32);

      hThrough_12_1->Fill(tot1,tof_21);
      hThrough_12_2->Fill(tot2,tof_21);
      hThrough_12_1comp_1->Fill(tot1,tof_21 - f1->Eval(tot1));
      hThrough_12_1comp_2->Fill(tot2,tof_21 - f1->Eval(tot1));
      hThrough_12_2comp_1->Fill(tot1,tof_21 - f1->Eval(tot1) - f2->Eval(tot2));
      hThrough_12_2comp_2->Fill(tot2,tof_21 - f1->Eval(tot1) - f2->Eval(tot2));

      hThrough_13_1->Fill(tot1,tof_31);
      hThrough_13_2->Fill(tot3,tof_31);
      hThrough_13_1comp_1->Fill(tot1,tof_31 - f3->Eval(tot1));
      hThrough_13_1comp_2->Fill(tot3,tof_31 - f3->Eval(tot1));
      hThrough_13_2comp_1->Fill(tot1,tof_31 - f3->Eval(tot1) - f4->Eval(tot3));
      hThrough_13_2comp_2->Fill(tot3,tof_31 - f3->Eval(tot1) - f4->Eval(tot3));

      hThrough_23_1->Fill(tot2,tof_32);
      hThrough_23_2->Fill(tot3,tof_32);
      hThrough_23_1comp_1->Fill(tot2,tof_32 - f5->Eval(tot2));
      hThrough_23_1comp_2->Fill(tot3,tof_32 - f5->Eval(tot2));
      hThrough_23_2comp_1->Fill(tot2,tof_32 - f5->Eval(tot2) - f6->Eval(tot3));
      hThrough_23_2comp_2->Fill(tot3,tof_32 - f5->Eval(tot2) - f6->Eval(tot3));

    }

  }



  // Fitting
  hThrough_12_1->Fit(f1,"Q");
  hThrough_12_1comp_2->Fit(f2,"Q");

  hThrough_13_1->Fit(f3,"Q");
  hThrough_13_1comp_2->Fit(f4,"Q");

  hThrough_23_1->Fit(f5,"Q");
  hThrough_23_1comp_2->Fit(f6,"Q");





  // Draw
  TCanvas *c0 = new TCanvas("c0", "tof", 3,2,1200,900);
  TCanvas *c1 = new TCanvas("c1", "throughing compensations for tof_21", 3,2,1200,900);
  TCanvas *c2 = new TCanvas("c2", "throughing compensations for tof_31", 3,2,1200,900);
  TCanvas *c3 = new TCanvas("c3", "throughing compensations for tof_32", 3,2,1200,900);

  c0->Divide(3,3);
  c0->cd(1);
  htof_21->Draw("");  
  c0->cd(2);
  htof_31->Draw("");
  c0->cd(3);
  htof_32->Draw("");
  c0->cd(4);
  hThrough_12_1comp_1->ProjectionY()->Draw();
  c0->cd(5);
  hThrough_13_1comp_1->ProjectionY()->Draw();
  c0->cd(6);
  hThrough_23_1comp_1->ProjectionY()->Draw();
  c0->cd(7);
  hThrough_12_2comp_1->ProjectionY()->Draw();
  hThrough_12_2comp_1->ProjectionY()->Fit(gaus12,"Q");
  c0->cd(8);
  hThrough_13_2comp_1->ProjectionY()->Draw();
  hThrough_13_2comp_1->ProjectionY()->Fit(gaus13,"Q");
  c0->cd(9);
  hThrough_23_2comp_1->ProjectionY()->Draw();
  hThrough_23_2comp_1->ProjectionY()->Fit(gaus23, "Q");

  c1->Divide(2,3);
  c1->cd(1);
  hThrough_12_1->Draw("colz");
  c1->cd(2);
  hThrough_12_2->Draw("colz");
  c1->cd(3);
  hThrough_12_1comp_1->Draw("colz");
  c1->cd(4);
  hThrough_12_1comp_2->Draw("colz");
  c1->cd(5);
  hThrough_12_2comp_1->Draw("colz");
  c1->cd(6);
  hThrough_12_2comp_2->Draw("colz");

  c2->Divide(2,3);
  c2->cd(1);
  hThrough_13_1->Draw("colz");
  c2->cd(2);
  hThrough_13_2->Draw("colz");
  c2->cd(3);
  hThrough_13_1comp_1->Draw("colz");
  c2->cd(4);
  hThrough_13_1comp_2->Draw("colz");
  c2->cd(5);
  hThrough_13_2comp_1->Draw("colz");
  c2->cd(6);
  hThrough_13_2comp_2->Draw("colz");

  c3->Divide(2,3);
  c3->cd(1);
  hThrough_23_1->Draw("colz");
  c3->cd(2);
  hThrough_23_2->Draw("colz");
  c3->cd(3);
  hThrough_23_1comp_1->Draw("colz");
  c3->cd(4);
  hThrough_23_1comp_2->Draw("colz");
  c3->cd(5);
  hThrough_23_2comp_1->Draw("colz");
  c3->cd(6);
  hThrough_23_2comp_2->Draw("colz");

  
  double sigma12 = gaus12->GetParameter("Sigma");
  double sigma13 = gaus13->GetParameter("Sigma");
  double sigma23 = gaus23->GetParameter("Sigma");
  cout << "sigma_1^2 + sigma_2^2 =" << sigma12*sigma12 << endl;
  cout << "sigma_1^2 + sigma_3^2 =" << sigma13*sigma13 << endl;
  cout << "sigma_1^2 + sigma_3^2 =" << sigma23*sigma23 << endl;
  cout << endl;
  cout << "sigma_1 = " << sqrt((sigma12*sigma12 + sigma13*sigma13 - sigma23*sigma23)/2) << "[ns]" << endl;
  cout << "sigma_2 = " << sqrt((sigma12*sigma12 + sigma23*sigma23 - sigma13*sigma13)/2) << "[ns]" << endl;
  cout << "sigma_3 = " << sqrt((sigma13*sigma13 + sigma23*sigma23 - sigma12*sigma12)/2) << "[ns]" << endl;


}
