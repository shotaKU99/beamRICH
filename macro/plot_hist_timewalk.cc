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
void SetCanvasGrid(TCanvas* cv, int npads) {
    if (npads == 1) {
        cv->SetGrid(1);
        cv->SetTicks(1);
    } else {
        for (int i = 1; i <= npads; i++) {
            cv->GetPad(i)->SetGrid(1);
            cv->GetPad(i)->SetTicks(1);
        }
    }
}

void plot_hist_timewalk(int run_no, int drawhist) {
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

    TFile* fin = new TFile(Form("/mnt/c/Users/shota/ELPHroot/data_corr/run%03d_correct_5.root",
     run_no));
    //TFile* fin =
    //    new TFile(Form("/mnt/c/Users/shota/ELPHroot/data_corr/run184_185_correct_4.root", run_no));
    TTree* tree = (TTree*)fin->Get("tree");

    int num_hist = 48;
    std::vector<TH1F*> htot_mppc(num_hist);
    std::vector<TH2F*> htof_tot_trig(num_hist);
    std::vector<TH2F*> htof_tot_calib0(num_hist);
    std::vector<TH2F*> htof_tot_calib1(num_hist);
    std::vector<TH2F*> htof_tot_calib2(num_hist);

    // -1: htot_mppc
    // 0 : htof_tot_trig2_calib
    // 1 : htof_tot_mppc_calib1
    // 2 : htof_tot_mppc_calib2
    int draw_hist = drawhist;

    if (draw_hist == -1) {
        TCanvas* cv_1_0 =
            new TCanvas("cv_1_0", "MPPC TOT (0-11)", 1200, 800);
        cv_1_0->Divide(4, 3);
        SetCanvasGrid(cv_1_0, 12);
        for (int i = 0; i < 12; i++) {
            cv_1_0->cd(i + 1);
            cv_1_0->GetPad(i + 1)->SetLogz(1);
            htot_mppc.at(i) = (TH1F*)fin->Get(Form("htot_mppc_%d", i));
            htot_mppc.at(i)->Draw();
        }

        TCanvas* cv_1_1 =
            new TCanvas("cv_1_1", "MPPC TOT (12-23)", 1200, 800);
        cv_1_1->Divide(4, 3);
        SetCanvasGrid(cv_1_1, 12);
        for (int i = 0; i < 12; i++) {
            cv_1_1->cd(i + 1);
            cv_1_1->GetPad(i + 1)->SetLogz(1);
            htot_mppc.at(i + 12) = (TH1F*)fin->Get(Form("htot_mppc_%d", i + 12));
            htot_mppc.at(i + 12)->Draw();
        }

        TCanvas* cv_1_2 =
            new TCanvas("cv_1_2", "MPPC TOT (24-35)", 1200, 800);
        cv_1_2->Divide(4, 3);
        SetCanvasGrid(cv_1_2, 12);
        for (int i = 0; i < 12; i++) {
            cv_1_2->cd(i + 1);
            cv_1_2->GetPad(i + 1)->SetLogz(1);
            htot_mppc.at(i + 24) = (TH1F*)fin->Get(Form("htot_mppc_%d", i + 24));
            htot_mppc.at(i + 24)->Draw();
        }

        TCanvas* cv_1_3 =
            new TCanvas("cv_1_3", "MPPC TOT (36-47)", 1200, 800);
        cv_1_3->Divide(4, 3);
        SetCanvasGrid(cv_1_3, 12);
        for (int i = 0; i < 12; i++) {
            cv_1_3->cd(i + 1);
            cv_1_3->GetPad(i + 1)->SetLogz(1);
            htot_mppc.at(i + 36) = (TH1F*)fin->Get(Form("htot_mppc_%d", i + 36));
            htot_mppc.at(i + 36)->Draw();
        }
    }else if (draw_hist == 0) {
        TCanvas* cv0_0 =
            new TCanvas("cv0_0", "TOF vs Trigger TOT (0-11) after Trigger calibration", 1200, 800);
        cv0_0->Divide(4, 3);
        SetCanvasGrid(cv0_0, 12);
        for (int i = 0; i < 12; i++) {
            cv0_0->cd(i + 1);
            cv0_0->GetPad(i + 1)->SetLogz(1);
            htof_tot_calib0.at(i) = (TH2F*)fin->Get(Form("htof_tot_trig2_calib_%d", i));
            htof_tot_calib0.at(i)->Draw("colz 0");
        }

        TCanvas* cv0_1 =
            new TCanvas("cv0_1", "TOF vs Trigger TOT (12-23) after Trigger calibration", 1200, 800);
        cv0_1->Divide(4, 3);
        SetCanvasGrid(cv0_1, 12);
        for (int i = 0; i < 12; i++) {
            cv0_1->cd(i + 1);
            cv0_1->GetPad(i + 1)->SetLogz(1);
            htof_tot_calib0.at(i + 12) = (TH2F*)fin->Get(Form("htof_tot_trig2_calib_%d", i + 12));
            htof_tot_calib0.at(i + 12)->Draw("colz 0");
        }

        TCanvas* cv0_2 =
            new TCanvas("cv0_2", "TOF vs Trigger TOT (24-35) after Trigger calibration", 1200, 800);
        cv0_2->Divide(4, 3);
        SetCanvasGrid(cv0_2, 12);
        for (int i = 0; i < 12; i++) {
            cv0_2->cd(i + 1);
            cv0_2->GetPad(i + 1)->SetLogz(1);
            htof_tot_calib0.at(i + 24) = (TH2F*)fin->Get(Form("htof_tot_trig2_calib_%d", i + 24));
            htof_tot_calib0.at(i + 24)->Draw("colz 0");
        }

        TCanvas* cv0_3 =
            new TCanvas("cv0_3", "TOF vs Trigger TOT (36-47) after Trigger calibration", 1200, 800);
        cv0_3->Divide(4, 3);
        SetCanvasGrid(cv0_3, 12);
        for (int i = 0; i < 12; i++) {
            cv0_3->cd(i + 1);
            cv0_3->GetPad(i + 1)->SetLogz(1);
            htof_tot_calib0.at(i + 36) = (TH2F*)fin->Get(Form("htof_tot_trig2_calib_%d", i + 36));
            htof_tot_calib0.at(i + 36)->Draw("colz 0");
        }
    } else if (draw_hist == 1) {
        TCanvas* cv1_0 =
            new TCanvas("cv1_0", "TOF vs TOT (0-11) after Trigger calibration", 1200, 800);
        cv1_0->Divide(4, 3);
        SetCanvasGrid(cv1_0, 12);
        for (int i = 0; i < 12; i++) {
            cv1_0->cd(i + 1);
            cv1_0->GetPad(i + 1)->SetLogz(1);
            htof_tot_calib1.at(i) = (TH2F*)fin->Get(Form("htof_tot_mppc_calib1_%d", i));
            htof_tot_calib1.at(i)->Draw("colz 0");
        }
        TCanvas* cv1_1 =
            new TCanvas("cv1_1", "TOF vs TOT (12-23) after Trigger calibration", 1200, 800);
        cv1_1->Divide(4, 3);
        SetCanvasGrid(cv1_1, 12);
        for (int i = 0; i < 12; i++) {
            cv1_1->cd(i + 1);
            cv1_1->GetPad(i + 1)->SetLogz(1);
            htof_tot_calib1.at(i + 12) = (TH2F*)fin->Get(Form("htof_tot_mppc_calib1_%d", i + 12));
            htof_tot_calib1.at(i + 12)->Draw("colz 0");
        }
        TCanvas* cv1_2 =
            new TCanvas("cv1_2", "TOF vs TOT (24-35) after Trigger calibration", 1200, 800);
        cv1_2->Divide(4, 3);
        SetCanvasGrid(cv1_2, 12);
        for (int i = 0; i < 12; i++) {
            cv1_2->cd(i + 1);
            cv1_2->GetPad(i + 1)->SetLogz(1);
            htof_tot_calib1.at(i + 24) = (TH2F*)fin->Get(Form("htof_tot_mppc_calib1_%d", i + 24));
            htof_tot_calib1.at(i + 24)->Draw("colz 0");
        }
        TCanvas* cv1_3 =
            new TCanvas("cv1_3", "TOF vs TOT (36-47) after Trigger calibration", 1200, 800);
        cv1_3->Divide(4, 3);
        SetCanvasGrid(cv1_3, 12);
        for (int i = 0; i < 12; i++) {
            cv1_3->cd(i + 1);
            cv1_3->GetPad(i + 1)->SetLogz(1);
            htof_tot_calib1.at(i + 36) = (TH2F*)fin->Get(Form("htof_tot_mppc_calib1_%d", i + 36));
            htof_tot_calib1.at(i + 36)->Draw("colz 0");
        }
    } else if (draw_hist == 2) {
        // TCanvas* cv2 = new TCanvas("cv2", "TOF vs TOT (12-23) after Trigger calibration", 1200,
        // 800); cv2->Divide(4,3); SetCanvasGrid(cv2, 12); for(int i=0;i<12;i++){
        //     cv2->cd(i+1);
        //     cv2->GetPad(i+1)->SetLogz(1);
        //     htof_tot_calib1.at(i+12) = (TH2F*)fin->Get(Form("htof_tot_mppc_calib1_%d",i+12));
        //     htof_tot_calib1.at(i+12) -> Draw("colz 0");
        // }

        TCanvas* cv3 = new TCanvas("cv3", "TOF vs TOT (0-11)", 1200, 800);
        cv3->Divide(4, 3);
        SetCanvasGrid(cv3, 12);
        for (int i = 0; i < 12; i++) {
            cv3->cd(i + 1);
            cv3->GetPad(i + 1)->SetLogz(1);
            htof_tot_calib2.at(i) = (TH2F*)fin->Get(Form("htof_tot_mppc_calib2_%d", i));
            htof_tot_calib2.at(i)->Draw("colz 0");
        }

        TCanvas* cv4 = new TCanvas("cv4", "TOF vs TOT (12-23)", 1200, 800);
        cv4->Divide(4, 3);
        SetCanvasGrid(cv4, 12);
        for (int i = 0; i < 12; i++) {
            cv4->cd(i + 1);
            cv4->GetPad(i + 1)->SetLogz(1);
            htof_tot_calib2.at(i + 12) = (TH2F*)fin->Get(Form("htof_tot_mppc_calib2_%d", i + 12));
            htof_tot_calib2.at(i + 12)->Draw("colz 0");
        }

        TCanvas* cv5 = new TCanvas("cv5", "TOF vs TOT (24-35)", 1200, 800);
        cv5->Divide(4, 3);
        SetCanvasGrid(cv5, 12);
        for (int i = 0; i < 12; i++) {
            cv5->cd(i + 1);
            cv5->GetPad(i + 1)->SetLogz(1);
            htof_tot_calib2.at(i + 24) = (TH2F*)fin->Get(Form("htof_tot_mppc_calib2_%d", i + 24));
            htof_tot_calib2.at(i + 24)->Draw("colz 0");
        }

        TCanvas* cv6 = new TCanvas("cv6", "TOF vs TOT (36-47)", 1200, 800);
        cv6->Divide(4, 3);
        SetCanvasGrid(cv6, 12);
        for (int i = 0; i < 12; i++) {
            cv6->cd(i + 1);
            cv6->GetPad(i + 1)->SetLogz(1);
            htof_tot_calib2.at(i + 36) = (TH2F*)fin->Get(Form("htof_tot_mppc_calib2_%d", i + 36));
            htof_tot_calib2.at(i + 36)->Draw("colz 0");
        }
    }
    // TCanvas* cv0_3 = new TCanvas("cv0_3", "TOF vs TOT", 800, 600);
    // SetCanvasGrid(cv0_3, 1);
    // cv0_3->SetLogz(1);
    // htof_tot_calib2.at(15)->SetStats(0);
    // htof_tot_calib2.at(15)->Draw("colz 0");
    // cv0_3->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231211/tof_tot_mppc_calib2_run%d_%d.pdf",run_no,
    // 15));
}