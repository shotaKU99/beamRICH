#include <math.h>

#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TMath.h"
#include "TStyle.h"
#include "TTree.h"
#include "TLegend.h"
#include "TGraphErrors.h"

void Vov_vs_nump(){
    
    //----------------------------------------//
    //-- Multiplicity, dark current, Photon --//
    //----------------------------------------//

    // run no. {107, 109, 111, 112, 113, 114}
    double Vov[6] = {3.0, 3.5, 4.0, 4.5, 5.0, 5.5};
    double multiplicity[6] = {11.3444, 12.215, 13.2026, 13.8056, 14.3501, 14.8853};
    double multiplicity_err[6] = {0.0162725, 0.0158371, 0.0174299, 0.00177925, 0.0181182, 0.0170698};

    double darkC[6] = {2.02025, 2.23565, 2.47295, 2.65609, 2.84631, 3.02689};
    double darkC_err[6] = {0.00666217, 0.0064475, 0.00710975, 0.00739246, 0.00771733, 0.00723895};

    double nump[6] = {9.32413, 9.97938, 10.7297, 11.1495, 11.5038, 11.8584};
    double nump_err[6] = {0.0175835, 0.0170992, 0.0188242, 0.0192671, 0.0196933, 0.0185413};

    TGraphErrors* gmulti = new TGraphErrors(6, Vov, multiplicity, 0, multiplicity_err);
    gmulti -> SetMarkerSize(1.5);
    gmulti -> SetMarkerColor(kRed);
    gmulti -> SetLineColor(kRed);
    gmulti -> SetMarkerStyle(8);

    TGraphErrors* gdark = new TGraphErrors(6, Vov, darkC, 0, darkC_err);
    gdark -> SetMarkerSize(1.5);
    gdark -> SetMarkerColor(kBlue);
    gdark -> SetLineColor(kBlue);
    gdark -> SetMarkerStyle(8);

    TGraphErrors* gphoton = new TGraphErrors(6, Vov, nump, 0, nump_err);
    gphoton -> SetMarkerSize(1.5);
    gphoton -> SetMarkerColor(kViolet);
    gphoton -> SetLineColor(kViolet);
    gphoton -> SetMarkerStyle(8);

    TH1I* htemp = new TH1I("htemp", ";V_{ov} [V];count", 3, 2.8, 5.8);
    htemp->SetStats(0);
    htemp->GetYaxis()->SetRangeUser(0, 17);
    double fontsize = 0.05;
    double tickfontsize = 0.04;
    htemp->GetXaxis()->SetTitleOffset(0.9);
    htemp->GetXaxis()->SetTitleSize(fontsize);
    htemp->GetXaxis()->SetLabelSize(tickfontsize);
    htemp->GetYaxis()->SetTitleOffset(1.);
    htemp->GetYaxis()->SetTitleSize(fontsize);
    htemp->GetYaxis()->SetLabelSize(tickfontsize);


    TCanvas* cv = new TCanvas("cv", "cv", 800, 600);
    cv->SetGrid(1);
    cv->SetTicks(1);
    cv->SetLeftMargin(0.14);
    cv->SetRightMargin(0.1);
    cv->SetTopMargin(0.08);
    cv->SetBottomMargin(0.16);


    htemp->Draw();
    gmulti->Draw("same P");
    gdark->Draw("same P");
    gphoton->Draw("same P");

    TLegend *legend = new TLegend(0.54, 0.4, 0.85, 0.6) ; //（）の中は位置の指定（左下の x , y 、右上の x , y ）
    legend->AddEntry( gmulti, "photon+darkcurrent" , "p") ; // AddEntry( pointer , "interpretation" , "option" )
    //legend->AddEntry( htheta2bg, "Dark Current" , "f") ; // option は　"f"=box, "l"="L"=line, "p"=marker
    legend->AddEntry(gdark, "darkcurrent", "p");
    legend->AddEntry(gphoton, "photon", "p");
    legend->SetFillColor(0);
    legend->Draw();

    //cv->Print("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231129/Vov_depend.pdf");


    //-------------------------------------//
    //--- Photon / (Photon+darkcurrent) ---//
    //-------------------------------------//

    double photonOvermulti[6];
    for(int i=0; i<6; i++){
        photonOvermulti[i] = 100.0*nump[i]/multiplicity[i];
    }

    TGraph* gphotonOvermulti = new TGraph(6, Vov, photonOvermulti);
    gphotonOvermulti -> SetMarkerSize(1.5);
    gphotonOvermulti -> SetMarkerColor(kRed);
    gphotonOvermulti -> SetLineColor(kRed);
    gphotonOvermulti -> SetMarkerStyle(8);

    TCanvas* cv2 = new TCanvas("cv2", "cv2", 800, 600);
    cv2->SetGrid(1);
    cv2->SetTicks(1);
    cv2->SetLeftMargin(0.14);
    cv2->SetRightMargin(0.1);
    cv2->SetTopMargin(0.08);
    cv2->SetBottomMargin(0.16);
    cv2->cd();

    TH1I* htemp2 = new TH1I("htemp2", ";V_{ov} [V];Number of photon / Multiplicity [%]", 3, 2.8, 5.8);
    htemp2->SetStats(0);
    htemp2->GetYaxis()->SetRangeUser(70, 100);
    htemp2->GetXaxis()->SetTitleOffset(0.9);
    htemp2->GetXaxis()->SetTitleSize(fontsize);
    htemp2->GetXaxis()->SetLabelSize(tickfontsize);
    htemp2->GetYaxis()->SetTitleOffset(1.);
    htemp2->GetYaxis()->SetTitleSize(fontsize);
    htemp2->GetYaxis()->SetLabelSize(tickfontsize);

    htemp2->Draw();
    gphotonOvermulti->Draw("same P");

    //cv2->Print("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231129/Vov_photonRatio.pdf");

    //------------------------//
    //-- Quantum Efficiency --//
    //------------------------//

    double queff[6];
    for(int i=0; i<6;i++){
        queff[i] = nump[i]/nump[0];
    }

    TGraph* gqueff = new TGraph(6, Vov, queff);
    gqueff -> SetMarkerSize(1.5);
    gqueff -> SetMarkerColor(kRed);
    gqueff -> SetLineColor(kRed);
    gqueff -> SetMarkerStyle(8);

    TCanvas* cv3 = new TCanvas("cv3", "cv3", 800, 600);
    cv3->SetGrid(1);
    cv3->SetTicks(1);
    cv3->SetLeftMargin(0.14);
    cv3->SetRightMargin(0.1);
    cv3->SetTopMargin(0.08);
    cv3->SetBottomMargin(0.16);
    cv3->cd();

    TH1I* htemp3 = new TH1I("htemp3", ";V_{ov} [V];# of photon / # of photon (V_{ov}=3 V)", 3, 2.8, 5.8);
    htemp3->SetStats(0);
    htemp3->GetYaxis()->SetRangeUser(0.8, 1.4);
    htemp3->GetXaxis()->SetTitleOffset(0.9);
    htemp3->GetXaxis()->SetTitleSize(fontsize);
    htemp3->GetXaxis()->SetLabelSize(tickfontsize);
    htemp3->GetYaxis()->SetTitleOffset(1.);
    htemp3->GetYaxis()->SetTitleSize(fontsize);
    htemp3->GetYaxis()->SetLabelSize(tickfontsize);

    htemp3->Draw();
    gqueff->Draw("same P");

    //cv3->Print("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231129/Vov_queff.pdf");
 



}