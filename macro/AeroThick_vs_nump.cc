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

void AeroThick_vs_nump(){
    
    //----------------------------------------//
    //-- Multiplicity, dark current, Photon --//
    //----------------------------------------//

    // run no. {175, 176, 177, 178, 180, 181, 183}
    double AeroThick[7] = {10.0, 8.0, 6.0, 5.0, 3.0, 2.0, 0.0};
    double multiplicity[7] = {14.2439, 12.7538, 11.0608, 10.4933, 7.99183, 6.56744, 2.58145};
    double multiplicity_err[7] = {0.0387916, 0.0357671, 0.0346669, 0.0213655, 0.0291573, 0.0255236, 0.0160014};

    double darkC[7] = {2.69721, 2.623, 2.61932, 2.6906, 2.64072, 2.75585, 2.58272};
    double darkC_err[7] = {0.0164439, 0.0164018, 0.0162028, 0.0104122, 0.0163016, 0.0159184, 0.0157694};

    double nump[7] = {11.44457, 10.1308, 8.44146, 7.80272, 5.35111, 3.81159, -0.0012751};
    double nump_err[7] = {0.042133, 0.0393485, 0.0382665, 0.0237676, 0.0334049, 0.0300807, 0.0224659};

    TGraphErrors* gmulti = new TGraphErrors(7, AeroThick, multiplicity, 0, multiplicity_err);
    gmulti -> SetMarkerSize(1.5);
    gmulti -> SetMarkerColor(kRed);
    gmulti -> SetLineColor(kRed);
    gmulti -> SetMarkerStyle(8);

    TGraphErrors* gdark = new TGraphErrors(7, AeroThick, darkC, 0, darkC_err);
    gdark -> SetMarkerSize(1.5);
    gdark -> SetMarkerColor(kBlue);
    gdark -> SetLineColor(kBlue);
    gdark -> SetMarkerStyle(8);

    TGraphErrors* gphoton = new TGraphErrors(7, AeroThick, nump, 0, nump_err);
    gphoton -> SetMarkerSize(1.5);
    gphoton -> SetMarkerColor(kViolet);
    gphoton -> SetLineColor(kViolet);
    gphoton -> SetMarkerStyle(8);

    TH1I* htemp = new TH1I("htemp", ";Aerogel Thickness [cm];count", 11, 0, 11);
    htemp->SetStats(0);
    htemp->GetYaxis()->SetRangeUser(-0.5, 17);
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

    TLegend *legend = new TLegend(0.15, 0.7, 0.45, 0.9) ; //（）の中は位置の指定（左下の x , y 、右上の x , y ）
    legend->AddEntry( gmulti, "photon+darkcurrent" , "p") ; // AddEntry( pointer , "interpretation" , "option" )
    //legend->AddEntry( htheta2bg, "Dark Current" , "f") ; // option は　"f"=box, "l"="L"=line, "p"=marker
    legend->AddEntry(gdark, "darkcurrent", "p");
    legend->AddEntry(gphoton, "photon", "p");
    legend->SetFillColor(0);
    legend->Draw();

    cv->Print("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231129/AeroThick_depend.pdf");


    //-------------------------------------//
    //--- Photon / (Photon+darkcurrent) ---//
    //-------------------------------------//

    double photonOvermulti[6];
    for(int i=0; i<6; i++){
        photonOvermulti[i] = 100.0*nump[i]/multiplicity[i];
    }

    TGraph* gphotonOvermulti = new TGraph(6, AeroThick, photonOvermulti);
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

    TH1I* htemp2 = new TH1I("htemp2", ";Aerogel Thickness [cm];Number of photon / Multiplicity [%]", 11, 0, 11);
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

    cv2->Print("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231129/AeroThick_photonRatio.pdf");


}