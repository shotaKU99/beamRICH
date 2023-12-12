#include <math.h>

#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraphErrors.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TLegend.h"
#include "TMath.h"
#include "TStyle.h"
#include "TTree.h"
#include "TMath.h"
#include "TChain.h"

void SetTH1Style(TH1* h) {
    double fontsize = 0.05;
    double tickfontsize = 0.04;
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetTitleSize(fontsize);
    h->GetXaxis()->SetLabelSize(tickfontsize);

    h->GetYaxis()->SetTitleOffset(1.);
    h->GetYaxis()->SetTitleSize(fontsize);
    h->GetYaxis()->SetLabelSize(tickfontsize);

    h->GetXaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetMaxDigits(3);
}
void SetTH2Style(TH2* h) {
    double fontsize = 0.05;
    double tickfontsize = 0.04;
    h->GetXaxis()->SetTitleOffset(0.9);
    h->GetXaxis()->SetTitleSize(fontsize);
    h->GetXaxis()->SetLabelSize(tickfontsize);

    h->GetYaxis()->SetTitleOffset(1.);
    h->GetYaxis()->SetTitleSize(fontsize);
    h->GetYaxis()->SetLabelSize(tickfontsize);

    h->GetXaxis()->SetMaxDigits(3);
    h->GetYaxis()->SetMaxDigits(3);
}
void SetCanvasGrid(TCanvas* cv, int npads){
    if(npads==1){
        cv->SetGrid(1);
        cv->SetTicks(1);
    }else{
        for(int i=1;i<=npads;i++){
            cv->GetPad(i)->SetGrid(1);
            cv->GetPad(i)->SetTicks(1);
        }
    }
}

void ltdc_ttdc4(int run_no) {
    double fontsize = 0.05;
    gStyle->SetTitleFontSize(fontsize);

    const Int_t NRGBs = 5;
    const Int_t NCont = 255;
    Double_t stops[NRGBs] = {0.00, 0.34, 0.61, 0.84, 1.00};
    Double_t Red[NRGBs] = {0.00, 0.00, 0.87, 1.00, 0.51};
    Double_t Green[NRGBs] = {0.00, 0.81, 1.00, 0.20, 0.00};
    Double_t Blue[NRGBs] = {0.51, 1.00, 0.12, 0.00, 0.00};
    TColor::CreateGradientColorTable(NRGBs, stops, Red, Green, Blue, NCont);
    gStyle->SetNumberContours(NCont);
    //  gStyle->SetOptStat(0);
    // gStyle->SetPadGridX(true);
    // gStyle->SetPadGridY(true);
    gStyle->SetPadLeftMargin(0.15);
    gStyle->SetPadRightMargin(0.15);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetTextFont(132);
    gStyle->SetTitleFont(132, "XYZ");
    gStyle->SetLabelFont(132, "XYZ");
    gStyle->SetStatH(4);
    gStyle->SetStatW(0.20);
    gStyle->SetStatFontSize(0.08);
    gStyle->SetLabelSize(0.05, "XY");
    gStyle->SetTitleSize(0.05, "XY");
    gStyle->SetTitleOffset(0.7, "X");
    gStyle->SetTitleOffset(1.3, "X");
    gStyle->SetTitleSize(0.05);

    TFile* fin = new TFile(Form("/mnt/c/Users/shota/ELPHroot/data/run%04d.root", run_no));
    TTree* tree = (TTree*)fin->Get("tree");

    //TChain* tree = new TChain("tree", "");
    //tree->Add("/mnt/c/Users/shota/ELPHroot/data/run0184.root");
    //tree->Add("/mnt/c/Users/shota/ELPHroot/data/run0185.root");

    std::vector<std::vector<double> >* tdcl = 0;  // NIM-EASIROC
    std::vector<std::vector<double> >* tdct = 0;  // NIM-EASIROC
    std::vector<std::vector<double> >* ltdc = 0;  // HR-TDC
    std::vector<std::vector<double> >* ttdc = 0;  // HR-TDC

    tree->SetBranchAddress("tdcl", &tdcl);
    tree->SetBranchAddress("tdct", &tdct);
    tree->SetBranchAddress("ltdc", &ltdc);
    tree->SetBranchAddress("ttdc", &ttdc);
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("tdcl", 1);
    tree->SetBranchStatus("tdct", 1);
    tree->SetBranchStatus("ltdc", 1);
    tree->SetBranchStatus("ttdc", 1);

    TH1F* hsizediff = new TH1F("hsizediff", "Size Diffrence;size diff;events", 13, -6.5, 6.5);
    SetTH1Style(hsizediff);

    TH1F* htdcl = new TH1F("htdcl", "MPPC ltdc;MPPC ltdc [ns];counts", 400, 0, 400);
    SetTH1Style(htdcl);

    TH1F* hltdc = new TH1F("hltdc", "Trigger2 ltdc;Trigger2 ltdc [ns];counts", 400, 0, 400);
    SetTH1Style(hltdc);

    TH1F* htot_mppc = new TH1F("htot_mppc", "MPPC TOT;TOT [ns];counts", 100, 0, 100);
    SetTH1Style(htot_mppc);

    TH1F* htot_trig2 = new TH1F("htot_trig2", "Trigger2 TOT;TOT [ns];counts", 100, -20, 80);
    SetTH1Style(htot_trig2);

    TH1F* htof = new TH1F("htof", "TOF;TOF [ns];counts", 1000, -500, 500);
    SetTH1Style(htof);

    TH2F* htof_tot_trig2 =
        new TH2F("htof_tot_trig2", "Time walk;Trig2 TOT [ns];TOF [ns]", 50, 0, 50, 100, -80, 20);
    htof_tot_trig2->SetStats(0);
    SetTH2Style(htof_tot_trig2);

    TH2F* htof_tot_mppc =
        new TH2F("htof_tot_mppc", "Time walk;MPPC TOT [ns];TOF [ns]", 80, 0, 80, 100, -80, 20);
    SetTH2Style(htof_tot_mppc);
    htof_tot_mppc->SetStats(0);

    //TF1 *fcalibMPPC = new TF1("f1","pol2",0.000001,80);
    //fcalibMPPC->SetParameters(-13.9014, -0.720844, 0.00832389);

    int size_tdcl, size_tdct, sizediff;
    int select_MPPC = 0; //MPPC
    int MPPC_channel = 8;
    int trig1ch = 33;
    int trig2ch = 37;
    int trig3ch = 20;
    double tdcCut_MPPClow = 200.0;
    double tdcCut_MPPCHigh = 400.0;
    double calib_hrtdc = 1.04002;

    int select_ch; //MPPC
    double tof_2_mppc;
    double tdcl_ch, tdct_ch, tot_mppc;  // MPPC
    double ltdc_ch, ttdc_ch, tot_trig;  // Trigger

    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);

        if (ltdc->at(trig2ch).size() == 1 && ttdc->at(trig2ch).size() == 1) {
            ltdc_ch = ltdc->at(trig2ch).at(0)*calib_hrtdc;
            ttdc_ch = ttdc->at(trig2ch).at(0)*calib_hrtdc;
            tot_trig = ltdc_ch - ttdc_ch;
            hltdc->Fill(ltdc_ch);
            htot_trig2->Fill(tot_trig);

            for(int k=0;k<MPPC_channel;k++){
                select_ch = select_MPPC*MPPC_channel+k;
                size_tdcl = tdcl->at(select_ch).size();
                size_tdct = tdct->at(select_ch).size();
                sizediff = size_tdcl - size_tdct;
                hsizediff->Fill(sizediff);
                if (sizediff == 0 || sizediff==1) {
                    for (int j = 0; j < size_tdct; j++) {
                        tdcl_ch = tdcl->at(select_ch).at(j + sizediff);
                        tdct_ch = tdct->at(select_ch).at(j);
                        tot_mppc = tdcl_ch - tdct_ch;
                        htdcl->Fill(tdcl->at(select_ch).at(j + sizediff));
                        htot_mppc->Fill(tot_mppc);

                        tof_2_mppc = ltdc_ch - tdcl->at(select_ch).at(j + sizediff);
                        htof->Fill(tof_2_mppc);
                        if (tdcCut_MPPClow < tdcl_ch && tdcl_ch < tdcCut_MPPCHigh) {
                            htof_tot_mppc->Fill(tot_mppc, tof_2_mppc);
                            htof_tot_trig2->Fill(tot_trig, tof_2_mppc);
                            //htof_tot_mppc_calib1->Fill(tot_mppc, tof_2_mppc-fcalibMPPC->Eval(tot_mppc));
                        }
                    }
                } else if(sizediff==-1){
                    for (int j = 0; j < size_tdcl; j++) {
                        tdcl_ch = tdcl->at(select_ch).at(j);
                        tdct_ch = tdct->at(select_ch).at(j);
                        tot_mppc = tdcl_ch - tdct_ch;
                        htdcl->Fill(tdcl->at(select_ch).at(j));
                        htot_mppc->Fill(tot_mppc);

                        tof_2_mppc = ltdc_ch - tdcl->at(select_ch).at(j);
                        htof->Fill(tof_2_mppc);
                        if (tdcCut_MPPClow < tdcl_ch && tdcl_ch < tdcCut_MPPCHigh) {
                            htof_tot_mppc->Fill(tot_mppc, tof_2_mppc);
                            htof_tot_trig2->Fill(tot_trig, tof_2_mppc);
                            //htof_tot_mppc_calib1->Fill(tot_mppc, tof_2_mppc-fcalibMPPC->Eval(tot_mppc));
                        }
                    }
                }
            }
        }
    }
    double peak_height1pe, peak_mean1pe, peak_mean2pe;
    TF1* dgaus_tot = new TF1("dgaus_tot","gaus(0)+gaus(3)", 20, 50);
    peak_height1pe = htot_mppc->GetMaximum();
    peak_mean1pe = htot_mppc->GetBinCenter(htot_mppc->GetMaximumBin());
    peak_mean2pe = 1.4*peak_mean1pe;
    dgaus_tot->SetParameters(peak_height1pe, peak_mean1pe, 3, peak_height1pe/8.0, peak_mean2pe, 4);
    dgaus_tot->SetParLimits(1, peak_mean1pe-2.0, peak_mean1pe+2.0);
    dgaus_tot->SetParLimits(2, 0, 5);
    dgaus_tot->SetParLimits(4, peak_mean2pe-3.0, peak_mean2pe+3.0);
    dgaus_tot->SetParLimits(5, 0, 5);
    dgaus_tot->SetNpx(500);
    htot_mppc->Fit(dgaus_tot, "","", peak_mean1pe-5.0, peak_mean2pe+5.0);
 

    //-------------------------------------//
    //-- slewing correction for trigger2 --//
    //-------------------------------------//

    double tot_min_trig = 0.0;
    double tot_max_trig = 24.0;
    double tot_width_trig = 3.0;
    int num_proj_trig = (tot_max_trig-tot_min_trig)/tot_width_trig;

    std::vector<TH1D*> htof_proj_trig(num_proj_trig);
    std::vector<TF1*> gaus_tof_trig(num_proj_trig);
    std::vector<double> tot_center_trig(num_proj_trig);
    std::vector<double> tot_sigma_trig(num_proj_trig, tot_width_trig/std::sqrt(12));
    std::vector<double> mean_gausfit_trig(num_proj_trig);
    std::vector<double> sigma_gausfit_trig(num_proj_trig);

    double tot_low_trig, tot_high_trig;
    for(int i=0; i<num_proj_trig; i++){
        tot_low_trig = tot_min_trig + i*tot_width_trig;
        tot_high_trig = tot_low_trig + tot_width_trig;
        htof_proj_trig.at(i) = htof_tot_trig2->ProjectionY(Form("htof_proj_trig%d",i), tot_low_trig, tot_high_trig);
        
        //gaus_tof_trig.at(i) = new TF1(Form("gaus_trig%d", i), "gaus(0)+pol0(3)", -80, 20);
        //gaus_tof_trig.at(i)->SetParameters(1000,-25,3,10);
        gaus_tof_trig.at(i) = new TF1(Form("gaus_trig%d", i), "gaus", -80, 20);
        gaus_tof_trig.at(i)->SetParameters(100,-28,2);
        gaus_tof_trig.at(i)->SetNpx(500);
        htof_proj_trig.at(i)->Fit(gaus_tof_trig.at(i),"Q", "", -40,-20);
        
        tot_center_trig.at(i) = tot_low_trig + tot_width_trig/2.0;
        mean_gausfit_trig.at(i) = gaus_tof_trig.at(i)->GetParameter(1);
        sigma_gausfit_trig.at(i) = gaus_tof_trig.at(i)->GetParameter(2);
    }

    TGraphErrors* gTOF_trig = new TGraphErrors(num_proj_trig, &tot_center_trig[0], &mean_gausfit_trig[0], &tot_sigma_trig[0], &sigma_gausfit_trig[0]);
    gTOF_trig->SetLineColor(kRed);
    gTOF_trig->SetMarkerColor(kRed);
    gTOF_trig->SetLineWidth(2);

    TF1 *f3 = new TF1("f3","[0] + [1]/sqrt(x-[2])",0.000001,80);
    f3->SetParameters(20.0, -680, -300);
    gTOF_trig->Fit("f3");


    //------------------------------------------------------------------------------//
    //-- make histograms TOF vs MPPC TOT after slewing correction of trigger2 TOT --//
    //------------------------------------------------------------------------------//

    TH2F* htof_tot_trig2_calib =
        new TH2F("htof_tot_trig2_calib", "Time walk;Trigger2 TOT [ns];TOF [ns]", 50, 0, 50, 125, -80, 20);
    SetTH2Style(htof_tot_trig2_calib);
    htof_tot_trig2_calib->SetStats(0);

    TH2F* htof_tot_mppc_calib1 =
        new TH2F("htof_tot_mppc_calib1", "Time walk;MPPC TOT [ns];TOF [ns]", 80, 0, 80, 100, -40, 40);
    SetTH2Style(htof_tot_mppc_calib1);
    htof_tot_mppc_calib1->SetStats(0);


    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);

        if (ltdc->at(trig2ch).size() == 1 && ttdc->at(trig2ch).size() == 1) {
            ltdc_ch = ltdc->at(trig2ch).at(0)*calib_hrtdc;
            ttdc_ch = ttdc->at(trig2ch).at(0)*calib_hrtdc;
            tot_trig = ltdc_ch - ttdc_ch;
            //hltdc->Fill(ltdc_ch);
            //htot_trig2->Fill(tot_trig);

            for(int k=0;k<MPPC_channel;k++){
                select_ch = select_MPPC*MPPC_channel+k;
                size_tdcl = tdcl->at(select_ch).size();
                size_tdct = tdct->at(select_ch).size();
                sizediff = size_tdcl - size_tdct;
                //hsizediff->Fill(sizediff);
                if (sizediff == 0 || sizediff==1) {
                    for (int j = 0; j < size_tdct; j++) {
                        tdcl_ch = tdcl->at(select_ch).at(j + sizediff);
                        tdct_ch = tdct->at(select_ch).at(j);
                        tot_mppc = tdcl_ch - tdct_ch;
                        //htdcl->Fill(tdcl->at(select_ch).at(j + sizediff));
                        //htot_mppc->Fill(tot_mppc);

                        tof_2_mppc = ltdc_ch - tdcl->at(select_ch).at(j + sizediff);
                        //htof->Fill(tof_2_mppc);
                        if (tdcCut_MPPClow < tdcl_ch && tdcl_ch < tdcCut_MPPCHigh) {
                            htof_tot_trig2_calib->Fill(tot_trig, tof_2_mppc - f3->Eval(tot_trig));
                            htof_tot_mppc_calib1->Fill(tot_mppc, tof_2_mppc - f3->Eval(tot_trig));
                            //htof_tot_mppc_calib1->Fill(tot_mppc, tof_2_mppc-fcalibMPPC->Eval(tot_mppc));
                        }
                    }
                } else if(sizediff==-1){
                    for (int j = 0; j < size_tdcl; j++) {
                        tdcl_ch = tdcl->at(select_ch).at(j);
                        tdct_ch = tdct->at(select_ch).at(j);
                        tot_mppc = tdcl_ch - tdct_ch;
                        //htdcl->Fill(tdcl->at(select_ch).at(j));
                        //htot_mppc->Fill(tot_mppc);

                        tof_2_mppc = ltdc_ch - tdcl->at(select_ch).at(j);
                        //htof->Fill(tof_2_mppc);
                        if (tdcCut_MPPClow < tdcl_ch && tdcl_ch < tdcCut_MPPCHigh) {
                            htof_tot_trig2_calib->Fill(tot_trig, tof_2_mppc-f3->Eval(tot_trig));
                            htof_tot_mppc_calib1->Fill(tot_mppc, tof_2_mppc-f3->Eval(tot_trig));
                            //htof_tot_mppc_calib1->Fill(tot_mppc, tof_2_mppc-fcalibMPPC->Eval(tot_mppc));
                        }
                    }
                }
            }
        }
    }


    //---------------------------------------//
    //--- TOF projection fitting for MPPC ---//
    //---------------------------------------//

    double tot_min = 10.0;
    double tot_max = 55.0;
    double tot_width = 5.0;
    int num_proj = (tot_max-tot_min)/tot_width;
    double tot_cut = 35.0;
    int num_doublegaus = (tot_max-tot_cut)/tot_width;
    int num_singlegaus = num_proj-num_doublegaus;

    std::vector<TH1D*> htof_proj(num_proj);
    std::vector<TF1*> gaus_tof(num_proj);
    std::vector<double> tot_center1(num_proj);
    std::vector<double> tot_center2;
    tot_center2.reserve(num_doublegaus+5);
    //std::vector<double> tot_sigma(num_proj, tot_width/std::sqrt(12));
    std::vector<double> tot_sigma1(num_proj, tot_width/std::sqrt(12));
    std::vector<double> tot_sigma2(num_doublegaus, tot_width/std::sqrt(12));
    std::vector<double> mean_gausfit1(num_proj);
    std::vector<double> mean_gausfit2;
    mean_gausfit2.reserve(num_doublegaus+5);
    //std::vector<double> mean_gausfit1;
    //mean_gausfit1.reserve(num_proj);
    std::vector<double> sigma_gausfit1(num_proj);
    std::vector<double> sigmaErr_gausfit1(num_proj);
    std::vector<double> sigma_gausfit2;
    std::vector<double> sigmaErr_gausfit2;
    sigma_gausfit2.reserve(num_doublegaus+5);
    sigmaErr_gausfit2.reserve(num_doublegaus+5);
    //std::vector<double> sigma_gausfit1;
    //sigma_gausfit1.reserve(num_proj);

    std::vector<double> fitrange_low{1, 0, -2, -4, -4, -6, -8, -8, -7, -7};
    std::vector<double> fitrange_high{15, 11, 11, 8, 8, 5, 5, 5, 5, 5};

    double tot_low, tot_high;
    for(int i=0; i<num_proj; i++){
        tot_low = tot_min + i*tot_width;
        tot_high = tot_low + tot_width;
        htof_proj.at(i) = htof_tot_mppc_calib1->ProjectionY(Form("htof_proj%d",i), tot_low, tot_high);
        if(i<num_singlegaus){
            gaus_tof.at(i) = new TF1(Form("gaus%d", i), "gaus", -80, 20);
            gaus_tof.at(i)->SetParameters(htof_proj.at(i)->GetMaximum(), 4, 2);
            gaus_tof.at(i)->SetParLimits(1, -2,10);
        }else{
            gaus_tof.at(i) = new TF1(Form("gaus%d", i), "gaus(0)+gaus(3)+pol0(6)", -80, 20);
            gaus_tof.at(i)->SetParameters(200, -3.5, 1., 20, 1, 1., 50);
            gaus_tof.at(i)->SetParLimits(1, -6.0,-2);
            gaus_tof.at(i)->SetParLimits(2, 0.4, 3);
            gaus_tof.at(i)->SetParLimits(4, -2, 2);
            gaus_tof.at(i)->SetParLimits(5, 0.4, 3);
            gaus_tof.at(i)->SetParLimits(6, 0, 1000);
        }
        gaus_tof.at(i)->SetNpx(500);
        htof_proj.at(i)->Fit(gaus_tof.at(i), "", "", fitrange_low.at(i), fitrange_high.at(i));
        
        if(i<num_singlegaus){
            tot_center1.at(i) = tot_low + tot_width/2.0;
            mean_gausfit1.at(i) = gaus_tof.at(i)->GetParameter(1);
            sigma_gausfit1.at(i) = gaus_tof.at(i)->GetParameter(2);
            sigmaErr_gausfit1.at(i) = gaus_tof.at(i)->GetParError(2);
        }else{
            //tot_center1.at(i) = tot_low + tot_width/2.0;
            //mean_gausfit1.at(i) = gaus_tof.at(i)->GetParameter(4);
            //sigma_gausfit1.at(i) = gaus_tof.at(i)->GetParameter(5);
            //sigmaErr_gausfit1.at(i) = gaus_tof.at(i)->GetParError(5);
            //
            //tot_center2.push_back(tot_low + tot_width/2.0);
            //mean_gausfit2.push_back(gaus_tof.at(i)->GetParameter(1));
            //sigma_gausfit2.push_back(gaus_tof.at(i)->GetParameter(2));
            //sigmaErr_gausfit2.push_back(gaus_tof.at(i)->GetParError(2));
            
            tot_center1.at(i) = tot_low + tot_width/2.0;
            mean_gausfit1.at(i) = gaus_tof.at(i)->GetParameter(1);
            sigma_gausfit1.at(i) = gaus_tof.at(i)->GetParameter(2);
            sigmaErr_gausfit1.at(i) = gaus_tof.at(i)->GetParError(2);
            
            tot_center2.push_back(tot_low + tot_width/2.0);
            mean_gausfit2.push_back(gaus_tof.at(i)->GetParameter(4));
            sigma_gausfit2.push_back(gaus_tof.at(i)->GetParameter(5));
            sigmaErr_gausfit2.push_back(gaus_tof.at(i)->GetParError(5));
        }
    }


    //---------------------------------//
    //-- slewing correction for MPPC --//
    //---------------------------------//

    TGraphErrors* gTOF1 = new TGraphErrors(num_proj, &tot_center1[0],  &mean_gausfit1[0], &tot_sigma1[0], &sigma_gausfit1[0]);
    gTOF1->SetLineColor(kRed);
    gTOF1->SetMarkerColor(kRed);
    gTOF1->SetLineWidth(2);
    
    TGraphErrors* gTOF2 = new TGraphErrors(num_doublegaus, &tot_center2[0],  &mean_gausfit2[0], &tot_sigma2[0], &sigma_gausfit2[0]);
    gTOF2->SetLineColor(kViolet);
    gTOF2->SetMarkerColor(kViolet);
    gTOF2->SetLineWidth(2);

    //definitions of fitting functions
    TF1 *f1pe = new TF1("f1pe","pol2",0.000001,80);
    //f1pe->SetParameters(-10.0, -0.8, 0.001);
    f1pe->SetLineColor(kViolet);
    //TF1 *f2 = new TF1("f2","[0] + [1]/sqrt([2]*(x-[3]))",0.000001,80);
    //f2->SetParameters(-7, -10, 1, -100);
    TF1 *f2pe = new TF1("f2pe","pol2",0.000001,80);
    TF1 *f2pelow = new TF1("f2pelow","pol2",0.000001,80);
    //f2pe->SetParameters(-10, -0.8, 0.001);

    //htof_tot_mppc->Fit("f1");
    gTOF1->Fit(f2pe, "","",tot_cut, tot_max);
    gTOF1->Fit(f2pelow, "+","",tot_min, tot_cut);
    //htof_tot_trig2->Fit("f2");
    gTOF2->Fit(f1pe);

    TGraphErrors* gsigma_gaus1 = new TGraphErrors(num_proj, &tot_center1[0], &sigma_gausfit1[0], &tot_sigma1[0], &sigmaErr_gausfit1[0]);
    gsigma_gaus1->SetLineWidth(2);
    gsigma_gaus1->SetLineColor(kRed);
    TF1* pol0_forf1sigma = new TF1("pol0_forf1sigma", "pol0", 0, 80);
    pol0_forf1sigma->SetParameter(0, 2);
    gsigma_gaus1->Fit(pol0_forf1sigma, "Q","",35, 60);

    TGraphErrors* gsigma_gaus2 = new TGraphErrors(num_doublegaus, &tot_center2[0], &sigma_gausfit2[0], &tot_sigma2[0], &sigmaErr_gausfit2[0]);
    gsigma_gaus2->SetLineWidth(2);
    gsigma_gaus2->SetLineColor(kRed);
    TF1* pol0_forf2sigma = new TF1("pol0_forf2sigma", "pol0", 0, 80);
    pol0_forf2sigma->SetParameter(0, 2);
    gsigma_gaus2->Fit(pol0_forf2sigma, "Q","",35, 60);



    //----------------------------------------------------//
    //-- make histograms of slewing correction MPPC TOT --//
    //----------------------------------------------------//

    double tot_1pe_mean = dgaus_tot->GetParameter(1);
    double tot_1pe_sigma = dgaus_tot->GetParameter(2);
    double tot_2pe_mean = dgaus_tot->GetParameter(4);
    double tot_2pe_sigma = dgaus_tot->GetParameter(5);
    //double tot_cutRange_low = tot_1pe_mean + tot_1pe_sigma*(tot_2pe_mean-tot_1pe_mean)/(tot_1pe_sigma+tot_2pe_sigma);
    htot_mppc->GetXaxis()->SetRangeUser(tot_1pe_mean, tot_2pe_mean);
    double tot_cutRange_low = htot_mppc->GetBinCenter(htot_mppc->GetMinimumBin())-htot_mppc->GetBinWidth(1)/2.0;
    htot_mppc->GetXaxis()->SetRange(0,0);
    std::cout << "tot_cutRange_low = " << tot_cutRange_low << std::endl;
    double tot_cutRange_high = tot_max;
    //double sigma_1pe = pol0_forf1sigma->GetParameter(0);
    //double sigma_2pe = pol0_forf2sigma->GetParameter(0);
    double sigma_1pe = pol0_forf1sigma->GetParameter(0);
    double sigma_2pe = pol0_forf2sigma->GetParameter(0);
    double sigma_add = sigma_1pe + sigma_2pe;

    TH2F* htof_tot_mppc_calib2 = // after all slewing correction
        new TH2F("htof_tot_mppc_calib2", "Time walk;MPPC TOT [ns];TOF [ns]", 80, 0, 80, 100, -40, 40);
    SetTH2Style(htof_tot_mppc_calib2);
    htof_tot_mppc_calib2->SetStats(0);


    double tof_diff_1pe, tof_diff_2pe;
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);

        if (ltdc->at(trig2ch).size() == 1 && ttdc->at(trig2ch).size() == 1) {
            ltdc_ch = ltdc->at(trig2ch).at(0)*calib_hrtdc;
            ttdc_ch = ttdc->at(trig2ch).at(0)*calib_hrtdc;
            tot_trig = ltdc_ch - ttdc_ch;
            //hltdc->Fill(ltdc_ch);
            //htot_trig2->Fill(tot_trig);

            for(int k=0;k<MPPC_channel;k++){
                select_ch = select_MPPC*MPPC_channel+k;
                size_tdcl = tdcl->at(select_ch).size();
                size_tdct = tdct->at(select_ch).size();
                sizediff = size_tdcl - size_tdct;
                //hsizediff->Fill(sizediff);
                if (sizediff == 0 || sizediff==1) {
                    for (int j = 0; j < size_tdct; j++) {
                        tdcl_ch = tdcl->at(select_ch).at(j + sizediff);
                        tdct_ch = tdct->at(select_ch).at(j);
                        tot_mppc = tdcl_ch - tdct_ch;
                        //htdcl->Fill(tdcl->at(select_ch).at(j + sizediff));
                        //htot_mppc->Fill(tot_mppc);

                        tof_2_mppc = ltdc_ch - tdcl->at(select_ch).at(j + sizediff);
                        //htof->Fill(tof_2_mppc);
                        if (tdcCut_MPPClow < tdcl_ch && tdcl_ch < tdcCut_MPPCHigh) {
                            // f1pe > f2pe (35 <= MPPC TOT <= 55ns)
                            //tof_diff_1pe = tof_2_mppc-f3->Eval(tot_trig)-f1pe->Eval(tot_mppc);
                            //tof_diff_2pe = tof_2_mppc-f3->Eval(tot_trig)-f2pe->Eval(tot_mppc);
                            tof_diff_1pe = tof_2_mppc-f3->Eval(tot_trig)-f1pe->Eval(tot_mppc);
                            tof_diff_2pe = tof_2_mppc-f3->Eval(tot_trig)-f2pe->Eval(tot_mppc);
                            if(tot_mppc<tot_cutRange_low){
                                //htof_tot_mppc_calib2->Fill(tot_mppc, tof_diff_1pe);
                                htof_tot_mppc_calib2->Fill(tot_mppc, tof_2_mppc-f3->Eval(tot_trig)-f2pelow->Eval(tot_mppc));
                            }else if(tot_cutRange_low<=tot_mppc && tot_mppc <= tot_cutRange_high){
                                if(tof_diff_1pe>0){
                                    htof_tot_mppc_calib2->Fill(tot_mppc, tof_diff_1pe);
                                }else if(tof_diff_2pe<0){
                                    htof_tot_mppc_calib2->Fill(tot_mppc, tof_diff_2pe);
                                }else{
                                    //if(tof_diff_1pe/sigma_1pe < tof_diff_2pe/sigma_2pe){
                                    if(std::abs(tof_diff_1pe)*sigma_2pe/sigma_add < std::abs(tof_diff_2pe)*sigma_1pe/sigma_add){
                                        htof_tot_mppc_calib2->Fill(tot_mppc, tof_diff_1pe);
                                    }else{
                                        htof_tot_mppc_calib2->Fill(tot_mppc, tof_diff_2pe);
                                    }
                                }
                            }else{
                                htof_tot_mppc_calib2->Fill(tot_mppc, tof_2_mppc-f3->Eval(tot_trig));
                            }
                        }
                    }
                } else if(sizediff==-1){
                    for (int j = 0; j < size_tdcl; j++) {
                        tdcl_ch = tdcl->at(select_ch).at(j);
                        tdct_ch = tdct->at(select_ch).at(j);
                        tot_mppc = tdcl_ch - tdct_ch;
                        //htdcl->Fill(tdcl->at(select_ch).at(j));
                        //htot_mppc->Fill(tot_mppc);

                        tof_2_mppc = ltdc_ch - tdcl->at(select_ch).at(j);
                        //htof->Fill(tof_2_mppc);
                        if (tdcCut_MPPClow < tdcl_ch && tdcl_ch < tdcCut_MPPCHigh) {
                            // f1pe > f2pe (35 <= MPPC TOT <= 55ns)
                            //tof_diff_1pe = tof_2_mppc-f3->Eval(tot_trig)-f1pe->Eval(tot_mppc);
                            //tof_diff_2pe = tof_2_mppc-f3->Eval(tot_trig)-f2pe->Eval(tot_mppc);
                            tof_diff_1pe = tof_2_mppc-f3->Eval(tot_trig)-f1pe->Eval(tot_mppc);
                            tof_diff_2pe = tof_2_mppc-f3->Eval(tot_trig)-f2pe->Eval(tot_mppc);
                            if(tot_mppc<tot_cutRange_low){
                                //htof_tot_mppc_calib2->Fill(tot_mppc, tof_diff_1pe);
                                htof_tot_mppc_calib2->Fill(tot_mppc, tof_diff_2pe);
                            }else if(tot_cutRange_low<=tot_mppc && tot_mppc <= tot_cutRange_high){
                                if(tof_diff_1pe>0){
                                    htof_tot_mppc_calib2->Fill(tot_mppc, tof_diff_1pe);
                                }else if(tof_diff_2pe<0){
                                    htof_tot_mppc_calib2->Fill(tot_mppc, tof_diff_2pe);
                                }else{
                                    //if(tof_diff_1pe/sigma_1pe < tof_diff_2pe/sigma_2pe){
                                    if(std::abs(tof_diff_1pe)*sigma_2pe/sigma_add < std::abs(tof_diff_2pe)*sigma_1pe/sigma_add){
                                        htof_tot_mppc_calib2->Fill(tot_mppc, tof_diff_1pe);
                                    }else{
                                        htof_tot_mppc_calib2->Fill(tot_mppc, tof_diff_2pe);
                                    }
                                }
                            }else{
                                htof_tot_mppc_calib2->Fill(tot_mppc, tof_2_mppc-f3->Eval(tot_trig));
                            }
                        }
                    }
                }
            }
        }
    }


    //----------------------------------//
    //-- Evaluate TOF time resolution --//
    //----------------------------------//

    TH1D* htof_proj_all = htof_tot_mppc_calib1->ProjectionY();
    TF1* gaus_proj_all = new TF1("gaus_proj_all","gaus(0)+pol0(3)", -80, 20);
    gaus_proj_all->SetParameters(5000, 0, 1, 200);
    gaus_proj_all->SetNpx(500);
    htof_proj_all->SetStats(0);
    htof_proj_all->Fit(gaus_proj_all);
    
    TH1D* htof_proj_all_calib2 = htof_tot_mppc_calib2->ProjectionY();
    TF1* gaus_proj_all_calib2 = new TF1("gaus_proj_all_calib2","gaus(0)+pol0(3)", -80, 20);
    gaus_proj_all_calib2->SetNpx(500);
    gaus_proj_all_calib2->SetParameters(5000, 0, 1, 200);
    htof_proj_all_calib2->SetStats(0);
    htof_proj_all_calib2->Fit(gaus_proj_all_calib2);
    //htof_proj_all_calib->Fit(gaus_proj_all_calib, "","", -5,5);




    //--------------------//
    //-- histogram draw --//
    //--------------------//

    TCanvas* cv0 = new TCanvas("cv0", "TOF vs TOT MPPC after trig2 slewing correction", 800, 600);
    SetCanvasGrid(cv0, 1);
    cv0->SetLogz(1);
    htof_tot_mppc_calib1->Draw("colz 0");
    //htof_proj_all->Draw("");
    gTOF1->Draw("same P");    
    gTOF2->Draw("same P");    

    TCanvas* cv0_1 = new TCanvas("cv0_1", "TOF vs TOT trigger2", 800, 600);
    SetCanvasGrid(cv0_1, 1);
    cv0_1->SetLogz(1);
    htof_tot_trig2->Draw("colz 0");
    gTOF_trig->Draw("same P");

    TCanvas* cv0_1_1 = new TCanvas("cv0_1_1", "TOF vs TOT trigger2 after trig2 correction", 800, 600);
    SetCanvasGrid(cv0_1_1, 1);
    cv0_1_1->SetLogz(1);
    htof_tot_trig2_calib->Draw("colz 0");

    TCanvas* cv0_2 = new TCanvas("cv0_2", "TOF vs MPPC TOT after all slewing correction", 800, 600);
    SetCanvasGrid(cv0_2, 1);
    cv0_2->SetLogz(1);
    htof_tot_mppc_calib2->Draw("colz 0");
    //htof_proj_all_calib2->Draw();
    
    //cv0->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231211/tof_tot_mppc_calib_run%d_%d.pdf",run_no, select_MPPC));
    //cv0_1->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231211/tof_tot_trig2_run%d_%d.pdf",run_no, select_MPPC));
    //cv0_1_1->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231211/tof_tot_trig2_calib_run%d_%d.pdf",run_no, select_MPPC));
    //cv0_2->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231211/tof_tot_mppc_calib2_run%d_%d.pdf",run_no, select_MPPC));

    TCanvas* cv = new TCanvas("cv", "Slewing Correction", 750, 1000);
    cv->Divide(2,3);
    SetCanvasGrid(cv, 6);
    //cv->SetLeftMargin(0.14);
    //cv->SetRightMargin(0.1);
    //cv->SetTopMargin(0.08);
    //cv->SetBottomMargin(0.16);
    //cv->SetLogy(0);

    cv->cd(1);
    cv->GetPad(1)->SetLogz(1);
    htof_tot_trig2->Draw("colz 0");
    gTOF_trig->Draw("same P");

    cv->cd(2);
    cv->GetPad(2)->SetLogz(1);
    htof_tot_mppc->Draw("colz 0");
    //htof_proj_all->Draw();
    
    cv->cd(3);
    cv->GetPad(3)->SetLogz(1);
    htof_tot_trig2_calib->Draw("colz 0");
    //htof_tot_mppc_calib1->ProjectionY()->Draw();
    
    cv->cd(4);
    cv->GetPad(4)->SetLogz(1);
    htof_tot_mppc_calib1->Draw("colz 0");
    gTOF1->Draw("same P");
    gTOF2->Draw("same P");
    f2pe->Draw("same");
    
    cv->cd(5);
    cv->GetPad(5)->SetLogz(1);
    htof_tot_mppc_calib2->Draw("colz 0");
    
    cv->cd(6);
    //htof_proj_all_calib->Draw("");
    htof_proj_all_calib2->Draw();

    TCanvas* cv2 = new TCanvas("cv2", "Basic information", 800, 800);
    cv2->Divide(3, 3);
    SetCanvasGrid(cv2, 9);
    // cv2->SetLeftMargin(0.14);
    // cv2->SetRightMargin(0.1);
    // cv2->SetTopMargin(0.08);
    // cv2->SetBottomMargin(0.16);
    
    TH1I* htemp1 = new TH1I("htemp1", "Sigma of 1p.e. TOF timing;MPPC TOT [ns];sigma [ns]", 80, 0, 80);
    htemp1->SetStats(0);
    htemp1->GetYaxis()->SetRangeUser(0, 5);
    TH1I* htemp2 = new TH1I("htemp2", "Sigma of 2p.e. TOF timing;MPPC TOT [ns];sigma [ns]", 80, 0, 80);
    htemp2->SetStats(0);
    htemp2->GetYaxis()->SetRangeUser(0, 5);

    cv2->cd(1)->SetLogy(1);
    hsizediff->Draw();
    cv2->cd(2);
    htof->Draw();
    cv2->cd(3);
    htdcl->Draw();
    cv2->cd(4);
    htot_mppc->Draw();
    cv2->cd(5);
    hltdc->Draw();
    cv2->cd(6);
    htot_trig2->Draw();
    cv2->cd(7);
    htemp1->Draw();
    gsigma_gaus1->Draw("same P");
    cv2->cd(8);
    htemp2->Draw();
    gsigma_gaus2->Draw("same P");

    TCanvas* cv0_3 = new TCanvas("cv0_3", "Sigma 1p.e.", 800, 600);
    SetCanvasGrid(cv0_3, 1);
    //htemp1->Draw();
    //gsigma_gaus1->Draw("same P");
    htof_proj_all_calib2->Draw();
    //cv0_3->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231211/sigma1pe_run%d_%d.pdf",run_no, select_MPPC));
    cv0_3->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231211/tof_proj_all_calib2_run%d_%d.pdf",run_no, select_MPPC));
    
    TCanvas* cv0_4 = new TCanvas("cv0_4", "Sigma 2p.e.", 800, 600);
    SetCanvasGrid(cv0_4, 1);
    htemp2->Draw();
    gsigma_gaus2->Draw("same P");
    //cv0_4->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231211/sigma2pe_run%d_%d.pdf",run_no, select_MPPC));
    




    TCanvas* cv3 = new TCanvas("cv3", "TOF fitting for MPPC TOT", 1200, 800);
    cv3->Divide(int(std::ceil(num_proj/2.0)),2);
    SetCanvasGrid(cv3, num_proj);

    for(int i=0;i<num_proj;i++){
        cv3->cd(i+1);
        htof_proj.at(i)->Draw();
    }
    //cv3->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231211/tof_mppc_run%d_%d.pdf",run_no, select_MPPC));
   

    TCanvas* cv4 = new TCanvas("cv4", "TOF fitting for Trigger TOT", 1200, 800);
    cv4->Divide(num_proj_trig/2,2);
    SetCanvasGrid(cv4, num_proj_trig);
   
    for(int i=0;i<num_proj_trig;i++){
        cv4->cd(i+1);
        htof_proj_trig.at(i)->Draw();
    }
    //cv4->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231211/tof_trig_run%d_%d.pdf",run_no, select_MPPC));

    

}