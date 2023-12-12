/*
MPPC 8Ch 毎に histgram に詰めて Time walk correction
補正関数は2本、TOT 1p.e.と真のTOT 2p.e.を１つの補正関数で補正したバージョン
*/


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
#include "TObjArray.h"
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

void slewingCorrection5(int run_no, int event_cal=-1) {
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

    //-------------//
    //-- Setting --//
    //-------------//

    int Num_MPPC = 48; // number of mppc
    int MPPC_channel = 8;
    int trig1ch = 33;
    int trig2ch = 37;
    int trig3ch = 20;
    double calib_hrtdc = 1.04002;

    int select_trigch = trig2ch;

    int nloop_event;
    if(event_cal<0){
        nloop_event = tree->GetEntries();
    }else if(event_cal < tree->GetEntries()+1){
        nloop_event = event_cal;
    }else{
        nloop_event = tree->GetEntries();
    }

    
    TH1F* htot_trig2 = new TH1F("htot_trig2", "Trigger2 TOT;TOT [ns];counts", 100, -20, 80);
    SetTH1Style(htot_trig2);

    std::vector<TH1F*> htot_mppc(Num_MPPC);
    std::vector<TH1F*> htof(Num_MPPC);
    std::vector<TH2F*> htof_tot_trig2(Num_MPPC);
    std::vector<TH2F*> htof_tot_mppc(Num_MPPC);
    std::vector<TH2F*> htof_tot_trig2_calib(Num_MPPC);
    std::vector<TH2F*> htof_tot_mppc_calib1(Num_MPPC);
    std::vector<TH2F*> htof_tot_mppc_calib2(Num_MPPC);

    for(int i=0; i<Num_MPPC; i++){
        htot_mppc.at(i) = new TH1F(Form("htot_mppc_%d", i), "MPPC TOT;TOT [ns];counts", 100, 0, 100);
        SetTH1Style(htot_mppc.at(i));
        htof.at(i) = new TH1F(Form("htof_%d", i), "TOF(trigger2-mppc);TOF [ns];counts", 100, 0, 100);
        SetTH1Style(htof.at(i));
        htof_tot_trig2.at(i) = new TH2F(Form("htof_tot_trig2_%d", i), "TOF vs TOT; Trigger2 TOT [ns];TOF [ns]", 50, 0, 50, 100, -80, 20);
        SetTH2Style(htof_tot_trig2.at(i));
        htof_tot_mppc.at(i) = new TH2F(Form("htof_tot_mppc_%d", i), "TOF vs TOT; MPPC TOT [ns];TOF [ns]", 80, 0, 80, 100, -80, 20);
        SetTH2Style(htof_tot_mppc.at(i));
        htof_tot_trig2_calib.at(i) = new TH2F(Form("htof_tot_trig2_calib_%d", i), "TOF vs TOT; Trigger2 TOT [ns];TOF [ns]", 50, 0, 50, 100, -50, 50);
        SetTH2Style(htof_tot_trig2_calib.at(i));
        htof_tot_mppc_calib1.at(i) = new TH2F(Form("htof_tot_mppc_calib1_%d", i), "TOF vs TOT; MPPC TOT [ns];TOF [ns]", 80, 0, 80, 100, -50, 50);
        SetTH2Style(htof_tot_mppc_calib1.at(i));
        htof_tot_mppc_calib2.at(i) = new TH2F(Form("htof_tot_mppc_calib2_%d", i), "TOF vs TOT; MPPC TOT [ns];TOF [ns]", 80, 0, 80, 100, -50, 50);
        SetTH2Style(htof_tot_mppc_calib2.at(i));
    }


    std::cout << "##########" << std::endl;
    std::cout << "-- Setting end --" << std::endl;
    std::cout << "##########" << std::endl;


    int size_tdcl, size_tdct, sizediff;
    //int select_MPPC = 10; //MPPC

    double tdcCut_MPPClow = 200.0;
    double tdcCut_MPPCHigh = 400.0;

    int select_ch; //MPPC
    double tof_2_mppc;
    double tdcl_ch, tdct_ch, tot_mppc;  // MPPC
    double ltdc_ch, ttdc_ch, tot_trig;  // Trigger

    std::vector<TF1*> dgaus_tot(Num_MPPC);

    for (int i = 0; i < nloop_event; i++) { // i event loop
        tree->GetEntry(i);

        if (ltdc->at(select_trigch).size() == 1 && ttdc->at(select_trigch).size() == 1) {
            ltdc_ch = ltdc->at(select_trigch).at(0)*calib_hrtdc;
            ttdc_ch = ttdc->at(select_trigch).at(0)*calib_hrtdc;
            tot_trig = ltdc_ch - ttdc_ch;
            htot_trig2->Fill(tot_trig);

            for(int m=0; m<Num_MPPC; m++){ // m number of MPPC (24) loop
                for(int k=0;k<MPPC_channel;k++){ // k MPPC channel (16) loop
                    select_ch = m*MPPC_channel+k;
                    size_tdcl = tdcl->at(select_ch).size();
                    size_tdct = tdct->at(select_ch).size();
                    sizediff = size_tdcl - size_tdct;
                    if (sizediff == 0 || sizediff==1) {
                        for (int j = 0; j < size_tdct; j++) { // j mppc hit loop of 1 channel
                            tdcl_ch = tdcl->at(select_ch).at(j + sizediff);
                            tdct_ch = tdct->at(select_ch).at(j);
                            tot_mppc = tdcl_ch - tdct_ch;
                            htot_mppc.at(m)->Fill(tot_mppc);

                            tof_2_mppc = ltdc_ch - tdcl->at(select_ch).at(j + sizediff);
                            htof.at(m)->Fill(tof_2_mppc);
                            if (tdcCut_MPPClow < tdcl_ch && tdcl_ch < tdcCut_MPPCHigh) {
                                htof_tot_mppc.at(m)->Fill(tot_mppc, tof_2_mppc);
                                htof_tot_trig2.at(m)->Fill(tot_trig, tof_2_mppc);
                            }
                        } // j loop end
                    } else if(sizediff==-1){
                        for (int j = 0; j < size_tdcl; j++) { // j mppc hit loop of 1 channel
                            tdcl_ch = tdcl->at(select_ch).at(j);
                            tdct_ch = tdct->at(select_ch).at(j);
                            tot_mppc = tdcl_ch - tdct_ch;
                            htot_mppc.at(m)->Fill(tot_mppc);

                            tof_2_mppc = ltdc_ch - tdcl->at(select_ch).at(j);
                            htof.at(m)->Fill(tof_2_mppc);
                            if (tdcCut_MPPClow < tdcl_ch && tdcl_ch < tdcCut_MPPCHigh) {
                                htof_tot_mppc.at(m)->Fill(tot_mppc, tof_2_mppc);
                                htof_tot_trig2.at(m)->Fill(tot_trig, tof_2_mppc);
                            }
                        } // j loop end
                    } // if sizediff end
                } // k loop end
            } // m loop end
        } // if trigger is one hit end
    } // i loop end

    double peak_height1pe, peak_mean1pe, peak_mean2pe;
    for(int m=0; m<Num_MPPC; m++){
        dgaus_tot.at(m) = new TF1(Form("dgaus_tot_%d",m),"gaus(0)+gaus(3)", 20, 50);
        peak_height1pe = htot_mppc.at(m)->GetMaximum();
        peak_mean1pe = htot_mppc.at(m)->GetBinCenter(htot_mppc.at(m)->GetMaximumBin());
        peak_mean2pe = 1.4*peak_mean1pe;
        dgaus_tot.at(m)->SetParameters(peak_height1pe, peak_mean1pe, 3, peak_height1pe/8.0, peak_mean2pe, 4);
        dgaus_tot.at(m)->SetParLimits(1, peak_mean1pe-2.0, peak_mean1pe+2.0);
        dgaus_tot.at(m)->SetParLimits(2, 0, 5);
        dgaus_tot.at(m)->SetParLimits(4, peak_mean2pe-3.0, peak_mean2pe+3.0);
        dgaus_tot.at(m)->SetParLimits(5, 0, 5);
        dgaus_tot.at(m)->SetNpx(500);
        htot_mppc.at(m)->Fit(dgaus_tot.at(m), "Q","", peak_mean1pe-5.0, peak_mean2pe+5.0);
    }


    std::cout << "##########" << std::endl;
    std::cout << "-- First Data Fill end --" << std::endl;
    std::cout << "##########" << std::endl;


    //-------------------------------------//
    //-- slewing correction for trigger2 --//
    //-------------------------------------//

    double tot_min_trig = 0.0;
    double tot_max_trig = 24.0;
    double tot_width_trig = 3.0;
    int num_proj_trig = (tot_max_trig-tot_min_trig)/tot_width_trig;

    std::vector<TH1D*> htof_proj_trig(num_proj_trig);
    std::vector<TF1*> gaus_tof_trig(num_proj_trig);
    std::vector<std::vector<double>> tot_center_trig(Num_MPPC , std::vector<double>(num_proj_trig));
    std::vector<std::vector<double>> tot_sigma_trig(Num_MPPC,std::vector<double>(num_proj_trig, tot_width_trig/std::sqrt(12)));
    std::vector<std::vector<double>> mean_gausfit_trig(Num_MPPC, std::vector<double>(num_proj_trig));
    std::vector<std::vector<double>> sigma_gausfit_trig(Num_MPPC, std::vector<double>(num_proj_trig));

    double tot_low_trig, tot_high_trig;
    std::vector<TGraphErrors*> gTOF_trig2(Num_MPPC);
    std::vector<TF1*> f3(Num_MPPC);
    for(int j=0; j<Num_MPPC; j++){
        for(int i=0; i<num_proj_trig; i++){
            tot_low_trig = tot_min_trig + i*tot_width_trig;
            tot_high_trig = tot_low_trig + tot_width_trig;
            htof_proj_trig.at(i) = htof_tot_trig2.at(j)->ProjectionY(Form("htof_proj_trig%d",i), tot_low_trig, tot_high_trig);

            //gaus_tof_trig.at(i) = new TF1(Form("gaus_trig%d", i), "gaus(0)+pol0(3)", -80, 20);
            gaus_tof_trig.at(i) = new TF1(Form("gaus_trig%d", i), "gaus", -80, 20);
            //gaus_tof_trig.at(i)->SetParameters(1000,-25,3,50);
            gaus_tof_trig.at(i)->SetParameters(1000,-25,3);
            gaus_tof_trig.at(i)->SetNpx(500);
            htof_proj_trig.at(i)->Fit(gaus_tof_trig.at(i),"Q", "", -40, -20);

            tot_center_trig.at(j).at(i) = tot_low_trig + tot_width_trig/2.0;
            mean_gausfit_trig.at(j).at(i) = gaus_tof_trig.at(i)->GetParameter(1);
            sigma_gausfit_trig.at(j).at(i) = gaus_tof_trig.at(i)->GetParameter(2);
        }
        
        gTOF_trig2.at(j) = new TGraphErrors(num_proj_trig, &tot_center_trig.at(j)[0], &mean_gausfit_trig.at(j)[0], &tot_sigma_trig.at(j)[0], &sigma_gausfit_trig.at(j)[0]);
        gTOF_trig2.at(j)->SetLineColor(kRed);
        gTOF_trig2.at(j)->SetMarkerColor(kRed);
        gTOF_trig2.at(j)->SetLineWidth(2);

        f3.at(j) = new TF1(Form("f3_%d",j),"[0] + [1]/sqrt(x-[2])",0.000001,80);
        f3.at(j)->SetParameters(20.0, -680, -300);
        gTOF_trig2.at(j)->Fit(f3.at(j), "Q");
    }


    //------------------------------------------------------------------------------//
    //-- make histograms TOF vs MPPC TOT after slewing correction of trigger2 TOT --//
    //------------------------------------------------------------------------------//


    for (int i = 0; i < nloop_event; i++) {
        tree->GetEntry(i);

        if (ltdc->at(select_trigch).size() == 1 && ttdc->at(select_trigch).size() == 1) {
            ltdc_ch = ltdc->at(select_trigch).at(0)*calib_hrtdc;
            ttdc_ch = ttdc->at(select_trigch).at(0)*calib_hrtdc;
            tot_trig = ltdc_ch - ttdc_ch;
            htot_trig2->Fill(tot_trig);

            for(int m=0; m<Num_MPPC; m++){ // m number of MPPC (24) loop
                for(int k=0;k<MPPC_channel;k++){
                    select_ch = m*MPPC_channel+k;
                    size_tdcl = tdcl->at(select_ch).size();
                    size_tdct = tdct->at(select_ch).size();
                    sizediff = size_tdcl - size_tdct;
                    if (sizediff == 0 || sizediff==1) {
                        for (int j = 0; j < size_tdct; j++) {
                            tdcl_ch = tdcl->at(select_ch).at(j + sizediff);
                            tdct_ch = tdct->at(select_ch).at(j);
                            tot_mppc = tdcl_ch - tdct_ch;

                            tof_2_mppc = ltdc_ch - tdcl->at(select_ch).at(j + sizediff);
                            if (tdcCut_MPPClow < tdcl_ch && tdcl_ch < tdcCut_MPPCHigh) {
                                htof_tot_trig2_calib.at(m)->Fill(tot_trig, tof_2_mppc - f3.at(m)->Eval(tot_trig));
                                htof_tot_mppc_calib1.at(m)->Fill(tot_mppc, tof_2_mppc - f3.at(m)->Eval(tot_trig));
                            }
                        }
                    } else if(sizediff==-1){
                        for (int j = 0; j < size_tdcl; j++) {
                            tdcl_ch = tdcl->at(select_ch).at(j);
                            tdct_ch = tdct->at(select_ch).at(j);
                            tot_mppc = tdcl_ch - tdct_ch;

                            tof_2_mppc = ltdc_ch - tdcl->at(select_ch).at(j);
                            if (tdcCut_MPPClow < tdcl_ch && tdcl_ch < tdcCut_MPPCHigh) {
                                htof_tot_trig2_calib.at(m)->Fill(tot_trig, tof_2_mppc-f3.at(m)->Eval(tot_trig));
                                htof_tot_mppc_calib1.at(m)->Fill(tot_mppc, tof_2_mppc-f3.at(m)->Eval(tot_trig));
                            }
                        } // j loop end
                    } // if sizediff end
                } // k loop end
            } // m loop end
        } // if trigger is one hit end
    } // i loop end

    std::cout << "##########" << std::endl;
    std::cout << "-- Second Data Fill (trigger slewing correction) end --" << std::endl;
    std::cout << "##########" << std::endl;

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
    std::vector<std::vector<double>> tot_center1(Num_MPPC, std::vector<double>(num_proj));
    std::vector<std::vector<double>> tot_center2(Num_MPPC, std::vector<double>(num_doublegaus));
    //tot_center2.reserve(num_doublegaus+5);
    //std::vector<std::vector<double>> tot_sigma(num_proj, tot_width/std::sqrt(12));
    std::vector<std::vector<double>> tot_sigma1(Num_MPPC, std::vector<double>(num_proj ,tot_width/std::sqrt(12)));
    std::vector<std::vector<double>> tot_sigma2(Num_MPPC, std::vector<double>(num_doublegaus, tot_width/std::sqrt(12)));
    std::vector<std::vector<double>> mean_gausfit1(Num_MPPC, std::vector<double>(num_proj));
    std::vector<std::vector<double>> mean_gausfit2(Num_MPPC, std::vector<double>(num_doublegaus));
    //mean_gausfit2.reserve(num_doublegaus+5);

    std::vector<std::vector<double>> sigma_gausfit1(Num_MPPC, std::vector<double>(num_proj));
    std::vector<std::vector<double>> sigmaErr_gausfit1(Num_MPPC, std::vector<double>(num_proj));
    std::vector<std::vector<double>> sigma_gausfit2(Num_MPPC, std::vector<double>(num_doublegaus));
    std::vector<std::vector<double>> sigmaErr_gausfit2(Num_MPPC, std::vector<double>(num_doublegaus));
    //sigma_gausfit2.reserve(num_doublegaus+5);
    //sigmaErr_gausfit2.reserve(num_doublegaus+5);


    std::vector<double> fitrange_low{1, 0, -2, -4, -4, -6, -8, -8, -7, -7};
    std::vector<double> fitrange_high{15, 11, 11, 8, 8, 5, 5, 5, 5, 5};

    double tot_low, tot_high;
    for(int j=0;j<Num_MPPC;j++){
        for(int i=0; i<num_proj; i++){
            tot_low = tot_min + i*tot_width;
            tot_high = tot_low + tot_width;
            htof_proj.at(i) = htof_tot_mppc_calib1.at(j)->ProjectionY(Form("htof_proj%d",i), tot_low, tot_high);
            if(i<num_singlegaus){
                gaus_tof.at(i) = new TF1(Form("gaus%d", i), "gaus", -80, 20);
                gaus_tof.at(i)->SetParameters(htof_proj.at(i)->GetMaximum(), 4, 2);
                gaus_tof.at(i)->SetParLimits(1, -2,10);
            }else{
                gaus_tof.at(i) = new TF1(Form("gaus%d", i), "gaus(0)+gaus(3)+pol0(6)", -80, 20);
                gaus_tof.at(i)->SetParameters(200, -3.5, 1., 20, 1, 1., 10);
                gaus_tof.at(i)->SetParLimits(1, -6.0,-2);
                gaus_tof.at(i)->SetParLimits(2, 0.4, 3);
                gaus_tof.at(i)->SetParLimits(4, -2, 2);
                gaus_tof.at(i)->SetParLimits(5, 0.4, 3);
                gaus_tof.at(i)->SetParLimits(6, 0, 1000);            
            }
            gaus_tof.at(i)->SetNpx(500);
            htof_proj.at(i)->Fit(gaus_tof.at(i), "Q", "", fitrange_low.at(i), fitrange_high.at(i));

            if(i<num_singlegaus){
                tot_center1.at(j).at(i) = tot_low + tot_width/2.0;
                mean_gausfit1.at(j).at(i) = gaus_tof.at(i)->GetParameter(1);
                sigma_gausfit1.at(j).at(i) = gaus_tof.at(i)->GetParameter(2);
                sigmaErr_gausfit1.at(j).at(i) = gaus_tof.at(i)->GetParError(2);
            }else{
                tot_center1.at(j).at(i) = tot_low + tot_width/2.0;
                mean_gausfit1.at(j).at(i) = gaus_tof.at(i)->GetParameter(1);
                sigma_gausfit1.at(j).at(i) = gaus_tof.at(i)->GetParameter(2);
                sigmaErr_gausfit1.at(j).at(i) = gaus_tof.at(i)->GetParError(2);

                tot_center2.at(j).at(i-num_singlegaus) = tot_low + tot_width/2.0;
                mean_gausfit2.at(j).at(i-num_singlegaus) = gaus_tof.at(i)->GetParameter(4);
                sigma_gausfit2.at(j).at(i-num_singlegaus) = gaus_tof.at(i)->GetParameter(5);
                sigmaErr_gausfit2.at(j).at(i-num_singlegaus) = gaus_tof.at(i)->GetParError(5);
            }
        }
    }

    //---------------------------------//
    //-- slewing correction for MPPC --//
    //---------------------------------//

    std::vector<TGraphErrors*> gTOF1(Num_MPPC);
    std::vector<TGraphErrors*> gTOF2(Num_MPPC);
    std::vector<TF1*> f1pe(Num_MPPC);
    std::vector<TF1*> f2pelow(Num_MPPC);
    std::vector<TF1*> f2pe(Num_MPPC);

    std::vector<TGraphErrors*> gsigma_gaus1(Num_MPPC); 
    std::vector<TGraphErrors*> gsigma_gaus2(Num_MPPC); 
    std::vector<TF1*> pol0_forf1sigma(Num_MPPC);
    std::vector<TF1*> pol0_forf2sigma(Num_MPPC);

    for(int i=0; i<Num_MPPC; i++){
        gTOF1.at(i) = new TGraphErrors(num_proj, &tot_center1.at(i)[0],  &mean_gausfit1.at(i)[0], &tot_sigma1.at(i)[0], &sigma_gausfit1.at(i)[0]);
        gTOF1.at(i)->SetLineColor(kRed);
        gTOF1.at(i)->SetMarkerColor(kRed);
        gTOF1.at(i)->SetLineWidth(2);

        gTOF2.at(i) = new TGraphErrors(num_doublegaus, &tot_center2.at(i)[0],  &mean_gausfit2.at(i)[0], &tot_sigma2.at(i)[0], &sigma_gausfit2.at(i)[0]);
        gTOF2.at(i)->SetLineColor(kViolet);
        gTOF2.at(i)->SetMarkerColor(kViolet);
        gTOF2.at(i)->SetLineWidth(2);
        
        f1pe.at(i) = new TF1(Form("f1pe_%d",i),"pol2",0.000001,80);
        f2pelow.at(i) = new TF1(Form("f1pelow_%d",i),"pol2",0.000001,80);
        //f1pe.at(i)->SetParameters(-10.0, -0.8, 0.001);
        f2pe.at(i) = new TF1(Form("f2pe_%d",i),"pol2",0.000001,80);
        //f2pe.at(i)->SetParameters(-10.0, -0.8, 0.001);
        f2pe.at(i)->SetLineColor(kViolet);

        //gTOF1.at(i)->Fit(f1pe.at(i), "Q");
        //gTOF2.at(i)->Fit(f2pe.at(i), "Q");
        gTOF1.at(i)->Fit(f2pe.at(i), "Q", "", tot_cut, tot_max);
        gTOF1.at(i)->Fit(f2pelow.at(i), "Q+", "", tot_min, tot_cut);
        gTOF2.at(i)->Fit(f1pe.at(i), "Q");
    
        gsigma_gaus2.at(i) = new TGraphErrors(num_proj, &tot_center1.at(i)[0], &sigma_gausfit1.at(i)[0], &tot_sigma1.at(i)[0], &sigmaErr_gausfit1.at(i)[0]);
        pol0_forf2sigma.at(i) = new TF1(Form("pol0_forf2sigma_%d",i), "pol0", 0, 80);
        pol0_forf2sigma.at(i)->SetParameter(0, 2);
        gsigma_gaus2.at(i)->Fit(pol0_forf2sigma.at(i), "Q","",35, 60);

        gsigma_gaus1.at(i) = new TGraphErrors(num_doublegaus, &tot_center2.at(i)[0], &sigma_gausfit2.at(i)[0], &tot_sigma2.at(i)[0], &sigmaErr_gausfit2.at(i)[0]);
        pol0_forf1sigma.at(i) = new TF1(Form("pol0_forf1sigma_%d",i), "pol0", 0, 80);
        pol0_forf1sigma.at(i)->SetParameter(0, 2);
        gsigma_gaus1.at(i)->Fit(pol0_forf1sigma.at(i), "Q","",35, 60);
    
    }


    //----------------------------------------------------//
    //-- make histograms of slewing correction MPPC TOT --//
    //----------------------------------------------------//

    double tot_1pe_mean, tot_2pe_mean;
    double tot_1pe_sigma, tot_2pe_sigma;
    double tot_cutRange_low;
    double tot_cutRange_high = tot_max;
    double sigma_1pe;
    double sigma_2pe;
    double sigma_add;
    double tof_value;

    double tof_cutrange_low = -5.0;// ns
    double tof_cutrange_high = 5.0;// ns


    //const char* outfilename = Form("/mnt/c/Users/shota/ELPHroot/data_corr/run%04d_correct_4.root",run_no);
    //const char* outfilename = "/mnt/c/Users/shota/ELPHroot/data_corr/run184_185_correct_4.root";
    std::string outfilename = "/mnt/c/Users/shota/ELPHroot/data_corr/run" + std::to_string(run_no) + "_correct_5.root";
    TFile* fout = TFile::Open(outfilename.c_str(),"RECREATE");
    TTree* tout = new TTree("tree", "Slewing correction");

    std::vector<std::vector<double>>* ltdc1 = 0;
    std::vector<std::vector<double>>* ttdc1 = 0;
    std::vector<std::vector<double>>* tdcl1 = 0;
    std::vector<std::vector<double>>* tdct1 = 0;
    std::vector<std::vector<double>>* width1 = 0;
    std::vector<std::vector<double>>* width_cut1 = 0;
    std::vector<std::vector<double>>* photonNum1 = 0;
    std::vector<std::vector<double>>* tof1 = 0;
    std::vector<double>* seg1 = 0;
    int numCh_hrtdc = 64;
    int numCh_EASI = 384;
    //int numCh_tof = 24; // = Num_MPPC
    int default_vectorsize = 30;

    tout->Branch("ltdc", &ltdc1);
    tout->Branch("ttdc", &ttdc1);
    tout->Branch("tdcl", &tdcl1);
    tout->Branch("tdct", &tdct1);
    tout->Branch("width", &width1); 
    tout->Branch("width_cut", &width_cut1); 
    tout->Branch("photonNum", &photonNum1); 
    tout->Branch("tof", &tof1);
    tout->Branch("seg", &seg1);

    ltdc1->reserve(numCh_hrtdc*default_vectorsize);
    ttdc1->reserve(numCh_hrtdc*default_vectorsize);
    tdcl1->reserve(numCh_EASI*default_vectorsize);
    tdct1->reserve(numCh_EASI*default_vectorsize);
    width1->reserve(numCh_EASI*default_vectorsize);
    width_cut1->reserve(numCh_EASI*default_vectorsize);
    photonNum1->reserve(numCh_EASI*default_vectorsize);
    tof1->reserve(numCh_EASI*default_vectorsize);
    seg1->reserve(numCh_EASI*default_vectorsize);
    
    ltdc1->resize(numCh_hrtdc);
    ttdc1->resize(numCh_hrtdc);
    tdcl1->resize(numCh_EASI);
    tdct1->resize(numCh_EASI);
    width1->resize(numCh_EASI);
    width_cut1->resize(numCh_EASI);
    photonNum1->resize(numCh_EASI);
    tof1->resize(numCh_EASI);
    
    //for(int n=0;n<numCh_hrtdc; n++){
    //    ltdc1->at(n).resize(default_vectorsize);
    //    ttdc1->at(n).resize(default_vectorsize);
    //}
    //for(int n=0;n<numCh_EASI; n++){
    //    tdcl1->at(n).resize(default_vectorsize);
    //    tdct1->at(n).resize(default_vectorsize);
    //    width1->at(n).resize(default_vectorsize);
    //    width_cut1->at(n).resize(default_vectorsize);
    //    tof1->at(n).resize(default_vectorsize);
    //}



    std::cout << "##########" << std::endl;
    std::cout << "-- Third Data Fill (MPPC slewing correction) start --" << std::endl;
    std::cout << "##########" << std::endl;

    double tof_diff_1pe, tof_diff_2pe;
    for (int i = 0; i < nloop_event; i++) {
        tree->GetEntry(i);
        if(i%5000==0){
            std::cout << "event num = " << i << std::endl;
        }

        //ltdc1->reserve(numCh_hrtdc*default_vectorsize);
        //ttdc1->reserve(numCh_hrtdc*default_vectorsize);
        //tdcl1->reserve(numCh_EASI*default_vectorsize);
        //tdct1->reserve(numCh_EASI*default_vectorsize);
        //width1->reserve(numCh_EASI*default_vectorsize);
        //tof1->reserve(numCh_EASI*default_vectorsize);
        //seg1->reserve(numCh_EASI*default_vectorsize);
        //
        //ltdc1->resize(numCh_hrtdc);
        //ttdc1->resize(numCh_hrtdc);
        //tdcl1->resize(numCh_EASI);
        //tdct1->resize(numCh_EASI);
        //width1->resize(numCh_EASI);
        //tof1->resize(numCh_EASI);
        
        for(int n=0;n<numCh_hrtdc; n++){
            //ltdc1->at(n).resize(default_vectorsize);
            //ttdc1->at(n).resize(default_vectorsize);
            ltdc1->at(n).clear();
            ttdc1->at(n).clear();
        }
        for(int n=0;n<numCh_EASI; n++){
            //tdcl1->at(n).resize(default_vectorsize);
            //tdct1->at(n).resize(default_vectorsize);
            //width1->at(n).resize(default_vectorsize);
            //tof1->at(n).resize(default_vectorsize);
            tdcl1->at(n).clear();
            tdct1->at(n).clear();
            width1->at(n).clear();
            width_cut1->at(n).clear();
            photonNum1->at(n).clear();
            tof1->at(n).clear();
        }
        seg1->clear();



        if (ltdc->at(select_trigch).size() == 1 && ttdc->at(select_trigch).size() == 1) {
            ltdc_ch = ltdc->at(select_trigch).at(0)*calib_hrtdc;
            ttdc_ch = ttdc->at(select_trigch).at(0)*calib_hrtdc;
            tot_trig = ltdc_ch - ttdc_ch;

            for(int n=0;n<numCh_hrtdc;n++){
                for(int vecsize=0;vecsize<ltdc->at(n).size();vecsize++){
                    ltdc1->at(n).push_back(ltdc->at(n).at(vecsize));
                }
                for(int vecsize=0;vecsize<ttdc->at(n).size();vecsize++){
                    ttdc1->at(n).push_back(ttdc->at(n).at(vecsize));
                }
            }

            for(int m=0; m<Num_MPPC; m++){ // m number of MPPC (24) loop
                tot_1pe_mean = dgaus_tot.at(m)->GetParameter(1);
                tot_1pe_sigma = dgaus_tot.at(m)->GetParameter(2);
                tot_2pe_mean = dgaus_tot.at(m)->GetParameter(4);
                tot_2pe_sigma = dgaus_tot.at(m)->GetParameter(5);
                //tot_1pe_mean = gaus_tot1pe.at(m)->GetParameter(1);
                //tot_1pe_sigma = gaus_tot1pe.at(m)->GetParameter(2);
                //tot_2pe_mean = gaus_tot2pe.at(m)->GetParameter(1);
                //tot_2pe_sigma = gaus_tot2pe.at(m)->GetParameter(2);
                htot_mppc.at(m)->GetXaxis()->SetRangeUser(tot_1pe_mean, tot_2pe_mean);
                tot_cutRange_low = htot_mppc.at(m)->GetBinCenter(htot_mppc.at(m)->GetMinimumBin())-htot_mppc.at(m)->GetBinWidth(1)/2.0;
                htot_mppc.at(m)->GetXaxis()->SetRange(0,0);
                //tot_cutRange_low = tot_1pe_mean+tot_1pe_sigma*(tot_2pe_mean-tot_1pe_mean)/(tot_1pe_sigma+tot_2pe_sigma);
                //std::cout << "tot_cut_range_low (MPPC " << m << ") = " << tot_cutRange_low << std::endl;
                for(int k=0;k<MPPC_channel;k++){ // k MPPC array ch (16) loop
                    select_ch = m*MPPC_channel+k;
                    size_tdcl = tdcl->at(select_ch).size();
                    size_tdct = tdct->at(select_ch).size();
                    sizediff = size_tdcl - size_tdct;
                    if (sizediff == 0 || sizediff==1) {
                        for (int j = 0; j < size_tdct; j++) { // j MPPC tdc hit loop
                            tdcl_ch = tdcl->at(select_ch).at(j + sizediff);
                            tdct_ch = tdct->at(select_ch).at(j);
                            tot_mppc = tdcl_ch - tdct_ch;


                            tof_2_mppc = ltdc_ch - tdcl->at(select_ch).at(j + sizediff);
                            if (tdcCut_MPPClow < tdcl_ch && tdcl_ch < tdcCut_MPPCHigh) {

                                tdcl1->at(select_ch).push_back(tdcl_ch);
                                tdct1->at(select_ch).push_back(tdct_ch);
                                width1->at(select_ch).push_back(tot_mppc);
                                width_cut1->at(select_ch).push_back(tot_cutRange_low);

                                tof_diff_1pe = tof_2_mppc-f3.at(m)->Eval(tot_trig)-f1pe.at(m)->Eval(tot_mppc);
                                tof_diff_2pe = tof_2_mppc-f3.at(m)->Eval(tot_trig)-f2pe.at(m)->Eval(tot_mppc);
                                sigma_1pe = pol0_forf1sigma.at(m)->GetParameter(0);
                                sigma_2pe = pol0_forf2sigma.at(m)->GetParameter(0);
                                sigma_add = sigma_1pe + sigma_2pe;

                                if(tot_mppc<tot_cutRange_low){
                                    tof_value = tof_2_mppc-f3.at(m)->Eval(tot_trig)-f2pelow.at(m)->Eval(tot_mppc);
                                    htof_tot_mppc_calib2.at(m)->Fill(tot_mppc, tof_value);
                                    tof1->at(select_ch).push_back(tof_value);
                                    photonNum1->at(select_ch).push_back(1);
                                    if(tof_cutrange_low<=tof_value && tof_value <= tof_cutrange_high){
                                        seg1->push_back(select_ch);
                                    }
                                }else if(tot_cutRange_low <= tot_mppc && tot_mppc <= tot_cutRange_high){
                                    if(tof_diff_1pe>0){
                                        tof_value = tof_diff_1pe;
                                        htof_tot_mppc_calib2.at(m)->Fill(tot_mppc, tof_value);
                                        tof1->at(select_ch).push_back(tof_value);
                                        if(std::abs(tof_diff_1pe-sigma_1pe)<0){
                                            photonNum1->at(select_ch).push_back(1);
                                        }else{
                                            photonNum1->at(select_ch).push_back(2);
                                        }
                                        if(tof_cutrange_low<=tof_value && tof_value <= tof_cutrange_high){
                                            seg1->push_back(select_ch);
                                        }
                                    }else if(tof_diff_2pe<0){
                                        tof_value = tof_diff_2pe;
                                        htof_tot_mppc_calib2.at(m)->Fill(tot_mppc,tof_value);
                                        tof1->at(select_ch).push_back(tof_value);
                                        photonNum1->at(select_ch).push_back(2);
                                        if(tof_cutrange_low<=tof_value && tof_value <= tof_cutrange_high){
                                            seg1->push_back(select_ch);
                                        }
                                    }else{
                                        if(std::abs(tof_diff_1pe)*sigma_2pe/sigma_add < std::abs(tof_diff_2pe)*sigma_1pe/sigma_add){
                                            tof_value = tof_diff_1pe;
                                            htof_tot_mppc_calib2.at(m)->Fill(tot_mppc, tof_value);
                                            tof1->at(select_ch).push_back(tof_value);
                                            photonNum1->at(select_ch).push_back(1);
                                            if(tof_cutrange_low<=tof_value && tof_value <= tof_cutrange_high){
                                                seg1->push_back(select_ch);
                                            }
                                        }else{
                                            tof_value = tof_diff_2pe;
                                            htof_tot_mppc_calib2.at(m)->Fill(tot_mppc,tof_value);
                                            tof1->at(select_ch).push_back(tof_value);
                                            photonNum1->at(select_ch).push_back(2);
                                            if(tof_cutrange_low<=tof_value && tof_value <= tof_cutrange_high){
                                                seg1->push_back(select_ch);
                                            }
                                        }
                                    }   
                                }else{
                                    tof_value = tof_2_mppc-f3.at(m)->Eval(tot_trig);
                                    htof_tot_mppc_calib2.at(m)->Fill(tot_mppc, tof_value);
                                    tof1->at(select_ch).push_back(tof_value);
                                    photonNum1->at(select_ch).push_back(2);
                                    if(tof_cutrange_low<=tof_value && tof_value <= tof_cutrange_high){
                                        seg1->push_back(select_ch);
                                    }
                                }
                            }
                        }
                    } else if(sizediff==-1){
                        for (int j = 0; j < size_tdcl; j++) {
                            tdcl_ch = tdcl->at(select_ch).at(j);
                            tdct_ch = tdct->at(select_ch).at(j);
                            tot_mppc = tdcl_ch - tdct_ch;

                            tof_2_mppc = ltdc_ch - tdcl->at(select_ch).at(j);
                            if (tdcCut_MPPClow < tdcl_ch && tdcl_ch < tdcCut_MPPCHigh) {

                                tdcl1->at(select_ch).push_back(tdcl_ch);
                                tdct1->at(select_ch).push_back(tdct_ch);
                                width1->at(select_ch).push_back(tot_mppc);
                                width_cut1->at(select_ch).push_back(tot_mppc);

                                tof_diff_1pe = tof_2_mppc-f3.at(m)->Eval(tot_trig)-f1pe.at(m)->Eval(tot_mppc);
                                tof_diff_2pe = tof_2_mppc-f3.at(m)->Eval(tot_trig)-f2pe.at(m)->Eval(tot_mppc);
                                sigma_1pe = pol0_forf1sigma.at(m)->GetParameter(0);
                                sigma_2pe = pol0_forf2sigma.at(m)->GetParameter(0);
                                sigma_add = sigma_1pe + sigma_2pe;

                                if(tot_mppc<tot_cutRange_low){
                                    tof_value = tof_diff_2pe;
                                    htof_tot_mppc_calib2.at(m)->Fill(tot_mppc, tof_value);
                                    tof1->at(select_ch).push_back(tof_value);
                                    photonNum1->at(select_ch).push_back(1);
                                    if(tof_cutrange_low<=tof_value && tof_value <= tof_cutrange_high){
                                        seg1->push_back(select_ch);
                                    }
                                }else if(tot_cutRange_low <= tot_mppc && tot_mppc <= tot_cutRange_high){
                                    if(tof_diff_1pe>0){
                                        tof_value = tof_diff_1pe;
                                        htof_tot_mppc_calib2.at(m)->Fill(tot_mppc, tof_value);
                                        tof1->at(select_ch).push_back(tof_value);
                                        if(std::abs(tof_diff_1pe-sigma_1pe)<0){
                                            photonNum1->at(select_ch).push_back(1);
                                        }else{
                                            photonNum1->at(select_ch).push_back(2);
                                        }                                        
                                        if(tof_cutrange_low<=tof_value && tof_value <= tof_cutrange_high){
                                            seg1->push_back(select_ch);
                                        }
                                    }else if(tof_diff_2pe<0){
                                        tof_value = tof_diff_2pe;
                                        htof_tot_mppc_calib2.at(m)->Fill(tot_mppc,tof_value);
                                        tof1->at(select_ch).push_back(tof_value);
                                        photonNum1->at(select_ch).push_back(2);
                                        if(tof_cutrange_low<=tof_value && tof_value <= tof_cutrange_high){
                                            seg1->push_back(select_ch);
                                        }
                                    }else{
                                        if(std::abs(tof_diff_1pe)*sigma_2pe/sigma_add < std::abs(tof_diff_2pe)*sigma_1pe/sigma_add){
                                            tof_value = tof_diff_1pe;
                                            htof_tot_mppc_calib2.at(m)->Fill(tot_mppc, tof_value);
                                            tof1->at(select_ch).push_back(tof_value);
                                            photonNum1->at(select_ch).push_back(1);
                                            if(tof_cutrange_low<=tof_value && tof_value <= tof_cutrange_high){
                                                seg1->push_back(select_ch);
                                            }
                                        }else{
                                            tof_value = tof_diff_2pe;
                                            htof_tot_mppc_calib2.at(m)->Fill(tot_mppc,tof_value);
                                            tof1->at(select_ch).push_back(tof_value);
                                            photonNum1->at(select_ch).push_back(2);
                                            if(tof_cutrange_low<=tof_value && tof_value <= tof_cutrange_high){
                                                seg1->push_back(select_ch);
                                            }
                                        }
                                    }  
                                }else{
                                    tof_value = tof_2_mppc-f3.at(m)->Eval(tot_trig);
                                    htof_tot_mppc_calib2.at(m)->Fill(tot_mppc, tof_value);
                                    tof1->at(select_ch).push_back(tof_value);
                                    photonNum1->at(select_ch).push_back(2);
                                    if(tof_cutrange_low<=tof_value && tof_value <= tof_cutrange_high){
                                        seg1->push_back(select_ch);
                                }
                                }
                            }
                        } // j loop end
                    } // if sizediff end
                } // k loop end
            } // m loop end
            tout->Fill();
        } // if trigger is one hit end
    } // i loop end

    std::cout << "##########" << std::endl;
    std::cout << "-- Thrid Data Fill (MPPC slewing correction) end --" << std::endl;
    std::cout << "##########" << std::endl;

    //----------------------------------//
    //-- Evaluate TOF time resolution --//
    //----------------------------------//
    
    std::vector<TH1D*> htof_proj_calib2(Num_MPPC);
    std::vector<TF1*> gaus_proj_calib2(Num_MPPC);
    for(int i=0;i<Num_MPPC;i++){
        htof_proj_calib2.at(i) = htof_tot_mppc_calib2.at(i)->ProjectionY(Form("htof_proj_calib2_%d",i));
        gaus_proj_calib2.at(i) = new TF1(Form("gaus_proj_calib2_%d",i),"gaus(0)+pol0(3)", -80, 20);
        gaus_proj_calib2.at(i)->SetNpx(500);
        gaus_proj_calib2.at(i)->SetParameters(5000, 0, 1, 200);
        htof_proj_calib2.at(i)->Fit(gaus_proj_calib2.at(i), "Q");
    }



    //----------------------//
    //-- histogram Output --//
    //----------------------//

    TObjArray* hist_tot_mppc = new TObjArray();
    TObjArray* hist_tof = new TObjArray();
    TObjArray* hist_tof_tot_trig2 = new TObjArray();
    TObjArray* hist_tof_tot_mppc = new TObjArray();
    TObjArray* hist_tof_tot_trig2_calib = new TObjArray();
    TObjArray* hist_tof_tot_mppc_calib1 = new TObjArray();
    TObjArray* hist_tof_tot_mppc_calib2 = new TObjArray();
    
    for(int i=0; i<Num_MPPC;i++){
        hist_tot_mppc->Add(htot_mppc.at(i));
        hist_tof->Add(htof.at(i));
        hist_tof_tot_trig2->Add(htof_tot_trig2.at(i));
        hist_tof_tot_mppc->Add(htof_tot_mppc.at(i));
        hist_tof_tot_trig2_calib->Add(htof_tot_trig2_calib.at(i));
        hist_tof_tot_mppc_calib1->Add(htof_tot_mppc_calib1.at(i));
        hist_tof_tot_mppc_calib2->Add(htof_tot_mppc_calib2.at(i));
    }
    hist_tot_mppc->Write();
    hist_tof->Write();
    hist_tof_tot_trig2->Write();
    hist_tof_tot_mppc->Write();
    hist_tof_tot_trig2_calib->Write();
    hist_tof_tot_mppc_calib1->Write();
    hist_tof_tot_mppc_calib2->Write();


    tout->Write();
    fout->Close();

    TCanvas* cv = new TCanvas("cv","cv");
    cv -> SetGrid(1);
    cv -> SetTicks(1);
    cv -> SetLogz(1);
    htof_tot_mppc_calib1.at(19)->Draw("colz 0");
    gTOF1.at(19)->Draw("same P");
    gTOF2.at(19)->Draw("same P");


    std::cout << "------- Analysis info --------" << std::endl;
    std::cout << "Run No. : " << run_no << std::endl;
    std::cout << "# of all events : " << tree->GetEntries() << std::endl;
    std::cout << "# of calcuration events : " << nloop_event << std::endl;
    std::cout << "Output root file name : " << outfilename << std::endl;
    std::cout << "------- Analysis info --------" << std::endl;


}