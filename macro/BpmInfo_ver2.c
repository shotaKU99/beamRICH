#include "TStyle.h"
#include "TCanvas.h"
#include "stdbool.h"



void BpmInfo_ver2( int run){

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

  // RF setting
  //int RFCh = 47;
  //int TagCh = 45;
  // double interval = 1.966 * 12;
  //int Nbunch = 420;


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




# if 0
  for( int i=0; i<ChMapX.size(); i++){
    cout << ChMapX[i] << endl;
  }
#endif



  
  



  
  TFile *fin = new TFile(Form("histo/run%04d.root",run));
  cout << Form("root/run00%2d",run) << " was opened" << endl;
  TTree *tree = (TTree*)fin->Get("tree");

  std::vector<std::vector<double>> *ltdc = 0;
  std::vector<std::vector<double>> *ttdc = 0;


    
  TBranch *Bltdc = 0;
  TBranch *Bttdc = 0;



  tree->SetBranchAddress("ltdc", &ltdc, &Bltdc );
  tree->SetBranchAddress("ttdc", &ttdc, &Bttdc );


  //definitions of histograms
  TH1F *hHitChX = new TH1F("hHitChX","X Hit Map; Ch No.",16,0,16);
  TH1F *hHitChY = new TH1F("hHitChY","Y Hit Map; Ch No.",16,0,16);
  TH2F *hCh_LvsR = new TH2F("Beam Profile from Downstream","; X Ch No.; Y Ch No. times (-1)",16,0,16,16,-16,0);//16,0,16,16,0,16);
  
  
  ////leading histograms
  TH1F *hltdcX[ChMapX.size()];
  for( int i=0; i<ChMapX.size(); i++){
    hltdcX[i] = new TH1F(Form("hltdcX%d",i),Form("ltdc X%d; ltdc [ns]",i),100,290,320);
  }
  
  TH1F *hltdcY[ChMapY.size()];
  for( int i=0; i<ChMapY.size(); i++){
    hltdcY[i] = new TH1F(Form("hltdcY%d",i),Form("ltdc Y%d; ltdc [ns]",i),100,290,320);
  }

  ////TOT histograms
  TH1F *htotX[ChMapX.size()];
  for( int i=0; i<ChMapX.size(); i++){
    htotX[i] = new TH1F(Form("htotX%d",i),Form("TOT X%d; TOT [ns]",i),100,-1,20);
  }

  TH1F *htotY[ChMapY.size()];
  for( int i=0; i<ChMapY.size(); i++){
    htotY[i] = new TH1F(Form("htotY%d",i),Form("TOT Y%d; TOT [ns]",i),100,-1,20);
  }


  ////ltdc of trigger
  TH1F *hltdc_trig[3]; 
  for( int i=0; i<3; i++){
    hltdc_trig[i]= new TH1F(Form("hltdc_trig[%d]",i), Form("ltdc_trig%d; ltdc [ns]",i+1),50,210,260);
  }

  //// hist of multiplicity
  TH1F *hMltptyX = new TH1F("hist of multiplicity X","hist of multiplicity X",5,-0.5,4.5);
  TH1F *hMltptyY = new TH1F("hist of multiplicity Y","hist of multiplicity Y",5,-0.5,4.5);

  //definitions of variables
  double AmpX;
  double AmpY;
  int MaxChAmpX, MaxAmpChY;
  double MaxAmpX, MaxAmpY;
  double leadL, leadR, leadTrig;
  double trailL, trailR, trailTrig;
  bool eventflagL, eventflagR;
  bool evflagL, evflagR;
  int ev=0, evL=0, evR=0;
  int sevL = 0, sevR = 0;
  int trev=0;

  
  const int n = tree->GetEntries();
  //const int n = 10000;

  
  for(int i=0; i<n; i++){
    if(i%10000 == 0) cout << "event num = " << i << endl;
    tree->GetEntry(i);
    std::vector<int> HitChX;
    std::vector<int> HitChY;

    if( ltdc->at(trig1ch).size()<1 || ttdc->at(trig1ch).size()<1 ) continue;
    if( ltdc->at(trig2ch).size()<1 || ttdc->at(trig2ch).size()<1 ) continue;
    if( ltdc->at(trig3ch).size()<1 || ttdc->at(trig3ch).size()<1 ) continue;
    //if( ltdc->at(trig4ch).size()<1 || ttdc->at(trig4ch).size()<1 ) continue;

    //if( ltdc->at(TagCh).size()<1 || ttdc->at(TagCh).size()<1 ) continue;
    //if( ltdc->at(TagCh).at(0)<80 || ltdc->at(TagCh).at(0)>140 ) continue;


    //    cout << i << endl;
    MaxAmpX = 0;
    MaxAmpY = 0;
    MaxChAmpX = 100;
    MaxAmpChY = 100;
    eventflagL = false;
    eventflagR = false;
    evflagL=false;
    evflagR=false;
    trev++;
      
    for(int j=0; j<ChMapX.size(); j++){
      if( ltdc->at(ChMapX[j]).size() < 1 || ttdc->at(ChMapX[j]).size() < 1 ) continue;

      leadL = 100;
      trailL = 100;
      leadL = ltdc->at(ChMapX[j]).at(0);
      trailL = ttdc->at(ChMapX[j]).at(0);
      AmpX = leadL - trailL;
      hltdcX[j]->Fill(leadL);
      if( AmpX<0 ) continue;
      //RF solution
      //      for(int k=0; k<Nbunch; k++){
      //        if( fabs(RF-leadL-fpeakL - interval*k)<interval/2 ){
      htotX[j]->Fill(AmpX);
      HitChX.push_back(j);
      if(j>3&&j<9){
	if(leadL>tdcrangelowL&&leadL<tdcrangehighL) eventflagL=true;
      }
      if(j==6&&leadL>tdcrangelowL&&leadL<tdcrangehighL) evflagL=true;
	  //  break;
	  //        }
	  //      }
      if( AmpX > MaxAmpX){
	MaxAmpX = AmpX;
	MaxChAmpX = j;
      }
    }

    for(int j=0; j<ChMapY.size(); j++){
      if( ltdc->at(ChMapY[j]).size() < 1 || ttdc->at(ChMapY[j]).size() < 1 ) continue;
      
      leadR = 100;
      trailR = 100;
      leadR = ltdc->at(ChMapY[j]).at(0);
      trailR = ttdc->at(ChMapY[j]).at(0);
      AmpY = leadR-trailR;
      hltdcY[j]->Fill(leadR);

      //RF solution
      //      for(int k=0; k<Nbunch; k++){
      //        if( fabs(RF-leadR-fpeakR - interval*k)<interval/2 ){
      htotY[j]->Fill(AmpY);
      if(j>3&&j<9){
	if(leadR>tdcrangelowR&&leadR<tdcrangehighR) eventflagR=true;
      }
      HitChY.push_back(j);
      if(j==6&&leadR>tdcrangelowR&&leadR<tdcrangehighR) evflagR=true;
      //  break;
      //        }
      //      }
      
      if( AmpY > MaxAmpY){
	MaxAmpY = AmpY;
	MaxAmpChY = j;
      }
    }
    

    for(int j=0; j<HitChX.size(); j++){
      hHitChX->Fill(HitChX[j]);
    }
    for(int j=0; j<HitChY.size(); j++){
      hHitChY->Fill(HitChY[j]);
    }




    ////fill hltdc_trig
    if( ltdc->at(trig1ch).size()  > 0 && ttdc->at(trig1ch).size() > 0 ){
      leadTrig = 100;
      trailTrig = 100;
      leadTrig = ltdc->at(trig1ch).at(0);
      trailTrig = ttdc->at(trig1ch).at(0);
      hltdc_trig[0]->Fill(leadTrig);
    };
    
    if( ltdc->at(trig2ch).size() > 0 && ttdc->at(trig2ch).size() > 0 ){
      leadTrig = 100;
      trailTrig = 100;
      leadTrig = ltdc->at(trig2ch).at(0);
      trailTrig = ttdc->at(trig2ch).at(0);
      hltdc_trig[1]->Fill(leadTrig);
    };

    if( ltdc->at(trig3ch).size() > 0 && ttdc->at(trig3ch).size() > 0 ){
    leadTrig = 100;
    trailTrig = 100;
    leadTrig = ltdc->at(trig3ch).at(0);
    trailTrig = ttdc->at(trig3ch).at(0);
    hltdc_trig[2]->Fill(leadTrig);
    };

    //// fill hist if multiplicity
    hMltptyX->Fill(HitChX.size());
    hMltptyY->Fill(HitChY.size());



    if(eventflagL) evL++;
    if(eventflagR) evR++;
    if(eventflagL && eventflagR) ev++;
    if(evflagL) sevL++;
    if(evflagR) sevR++;
    if(HitChX.size()==1 && HitChY.size()==1 ) hCh_LvsR->Fill(MaxChAmpX, -MaxAmpChY);
  }

  // ltdc left canvas
  TCanvas *c0 = new TCanvas("c0","ltdc of BPM X",3,2,1200,900);
  c0->Divide(4,3);  
  for( int i=0; i<ChMapX.size(); i++){
    c0->cd(i+1);
    hltdcX[i]->Draw("");
    hltdcX[i]->GetYaxis()->SetLabelSize(0.05);
    hltdcX[i]->GetXaxis()->SetLabelSize(0.07);
    hltdcX[i]->GetXaxis()->SetTitleSize(0.06);
    hltdcX[i]->GetXaxis()->SetTitleOffset(0.9);
  }
  // ltdc right canvas
  TCanvas *c1 = new TCanvas("c1","ltdc of BPM Y",33,32,1200,900);
  c1->Divide(4,3);  
  for( int i=0; i<ChMapY.size(); i++){
    c1->cd(i+1);
    hltdcY[i]->Draw("");
    hltdcY[i]->GetYaxis()->SetLabelSize(0.05);
    hltdcY[i]->GetXaxis()->SetLabelSize(0.05);
    hltdcY[i]->GetXaxis()->SetTitleSize(0.06);
    hltdcY[i]->GetXaxis()->SetTitleOffset(0.9);
  }

  // tot left canvas
  TCanvas *c2 = new TCanvas("c2","tot of BPM X",63,62,1200,900);
  c2->Divide(4,3);  
  for( int i=0; i<ChMapX.size(); i++){
    c2->cd(i+1);
    htotX[i]->Draw("");
    htotX[i]->GetYaxis()->SetLabelSize(0.05);
    htotX[i]->GetXaxis()->SetLabelSize(0.07);
    htotX[i]->GetXaxis()->SetTitleSize(0.06);
    htotX[i]->GetXaxis()->SetTitleOffset(0.9);
  }

  // tot right canvas
  TCanvas *c3 = new TCanvas("c3","tot of BPM Y",93,92,1200,900);
  c3->Divide(4,3);  
  for( int i=0; i<ChMapY.size(); i++){
    c3->cd(i+1);
    htotY[i]->Draw("");
    htotY[i]->GetYaxis()->SetLabelSize(0.05);
    htotY[i]->GetXaxis()->SetLabelSize(0.07);
    htotY[i]->GetXaxis()->SetTitleSize(0.06);
    htotY[i]->GetXaxis()->SetTitleOffset(0.9);
  }
  TCanvas *c4 = new TCanvas("c4","c4",3,2,1200,600);
  c4->Divide(4,2);

  c4->cd(1);
  hHitChX->Draw("");
  hHitChX->GetYaxis()->SetLabelSize(0.05);
  hHitChX->GetXaxis()->SetLabelSize(0.07); 
  hHitChX->GetXaxis()->SetTitleSize(0.06);
  hHitChX->GetXaxis()->SetTitleOffset(0.9);  

  c4->cd(2);
  hHitChY->Draw("");
  hHitChY->GetYaxis()->SetLabelSize(0.05);
  hHitChY->GetXaxis()->SetLabelSize(0.07); 
  hHitChY->GetXaxis()->SetTitleSize(0.06);
  hHitChY->GetXaxis()->SetTitleOffset(0.9);    

  c4->cd(3);
  hCh_LvsR->Draw("colz");
  hCh_LvsR->GetYaxis()->SetLabelSize(0.05);
  hCh_LvsR->GetXaxis()->SetLabelSize(0.07); 
  hCh_LvsR->GetXaxis()->SetTitleSize(0.06);
  hCh_LvsR->GetXaxis()->SetTitleOffset(0.9);
  hCh_LvsR->SetStats(0);
  hCh_LvsR->SetTitle("from Downstream");

  c4->cd(4);
  hMltptyY->Draw("box");
  hMltptyX->Draw("same");
  hMltptyX->SetStats(0);
  hMltptyY->SetLineColor(2);
  hMltptyX->GetXaxis()->SetNdivisions(5);
  TLegend *lg = new TLegend(0.7,0.7,0.9,0.9);
  lg->AddEntry(hMltptyX, "BPM X", "l");
  lg->AddEntry(hMltptyY, "BPM Y", "l");
  lg->Draw();
  hMltptyX->SetTitle("Hist of Multiplicity");
  

  c4->cd(5);
  hltdc_trig[0]->Draw("box");
  hltdc_trig[0]->GetYaxis()->SetLabelSize(0.05);
  hltdc_trig[0]->GetXaxis()->SetLabelSize(0.07);
  hltdc_trig[0]->GetXaxis()->SetTitleSize(0.06);
  hltdc_trig[0]->GetXaxis()->SetTitleOffset(0.9);
  hltdc_trig[0]->GetXaxis()->SetNdivisions(5);

  c4->cd(6);
  hltdc_trig[1]->Draw("box");
  hltdc_trig[1]->GetYaxis()->SetLabelSize(0.05);
  hltdc_trig[1]->GetXaxis()->SetLabelSize(0.07);
  hltdc_trig[1]->GetXaxis()->SetTitleSize(0.06);
  hltdc_trig[1]->GetXaxis()->SetTitleOffset(0.9);
  hltdc_trig[1]->GetXaxis()->SetNdivisions(5);

  c4->cd(7);
  hltdc_trig[2]->Draw("box");
  hltdc_trig[2]->GetYaxis()->SetLabelSize(0.05);
  hltdc_trig[2]->GetXaxis()->SetLabelSize(0.07);
  hltdc_trig[2]->GetXaxis()->SetTitleSize(0.06);
  hltdc_trig[2]->GetXaxis()->SetTitleOffset(0.9);
  hltdc_trig[2]->GetYaxis()->SetNdivisions(10);
  hltdc_trig[2]->GetXaxis()->SetNdivisions(5);

  cout << "trigger event = " << trev << endl;
  cout << "event num = " << ev << endl;
  cout << "left event num = " << evL << endl;
  cout << "right event num = " << evR << endl;
  cout << "single strip efficiency L = " << (double) sevL/trev*100 << endl;
  cout << "efficiencyL = " << (double) evL/trev*100 << endl;
  cout << "single strip efficiency R = " << (double) sevR/trev*100 << endl;
  cout << "efficiencyR = " << (double) evR/trev*100 << endl;




}//macro end
