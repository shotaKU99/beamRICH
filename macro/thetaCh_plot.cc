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

void thetaCh_plot(std::string inputfilename) {
    const int ch = 16 * 2 * 12;
    const int trig1 = 33;
    const int trig2 = 37;
    const int trig3 = 20;

    TFile* fin = new TFile(inputfilename.c_str());
    //TFile* fin = new TFile(Form("/mnt/c/Users/shota/ELPHroot/data/run%04d.root", run_no));
    TTree* tree = (TTree*)fin->Get("tree");

    // Get Run Number from root file name
    int run_no = -1;
    std::string find_word = "run";
    std::string::size_type word_pos = inputfilename.rfind(find_word);
    //std::cout << "Input file name = " << inputfilename << std::endl;
    //std::cout << "find word = " << find_word << std::endl;
    //std::cout << "find word position = " << word_pos << std::endl;
    //std::cout << "Word of finc word position = " << inputfilename.at(word_pos) << std::endl;
    //std::cout << "size_type npos = " << std::string::npos << std::endl;
    if(word_pos!=std::string::npos){
        run_no = std::stoi(inputfilename.substr(word_pos+find_word.length(), 3));
        //std::cout << "Run number = " << run_no << std::endl;
    }

    std::vector<double>* thetaCh = 0;
    //std::vector<std::vector<double> >* ltdc = 0;

    //TBranch* btdcl = 0;
    //TBranch* bltdc = 0;

    tree->SetBranchAddress("thetaCh", &thetaCh);
    //tree->SetBranchAddress("ltdc", &ltdc, &bltdc);
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("thetaCh", 1);
    //tree->SetBranchStatus("ltdc", 1);

    TH1F* h1 = new TH1F("h1", Form("Run No. %d;Cherenkov Angle [mrad];count",run_no), 160, 160, 240);
    //h1->SetLineColor(kBlack);
    //h1->SetFillColorAlpha(kRed, 0.3);
    
    //TH1F* h2 = new TH1F("h2", Form("Run No. %d;leading edge [ns];count",run_no), 150, 250, 400);
    //h2->SetLineColor(kBlack);
    //h2->SetFillColorAlpha(kRed, 0.3);
    //TH1F* h3 = new TH1F("h3", Form("Run No. %d;leading edge [ns];count",run_no), 150, 250, 400);
    //h3->SetLineColor(kBlack);
    //h3->SetFillColorAlpha(kBlue, 0.3);

    double fontsize = 0.05;
    double tickfontsize = 0.04;
    h1->GetXaxis()->SetTitleOffset(0.9);
    h1->GetXaxis()->SetTitleSize(fontsize);
    h1->GetXaxis()->SetLabelSize(tickfontsize);
    
    h1->GetYaxis()->SetTitleOffset(1.);
    h1->GetYaxis()->SetTitleSize(fontsize);
    h1->GetYaxis()->SetLabelSize(tickfontsize);

    gStyle->SetTitleFontSize(fontsize);

    h1->GetXaxis()->SetMaxDigits(3);
    h1->GetYaxis()->SetMaxDigits(3);


    // TString f1 = "../image/multiplicity.pdf";

    int counter = 0;
    int dccounter = 0;
    bool hit = false;
    bool dchit = false;

    gStyle->SetOptStat(0);

    for (int i = 0; i < tree->GetEntries(); ++i) {
        tree->GetEntry(i);
        //counter = 0;
        //dccounter = 0;
        if (!thetaCh->empty()) {
            for (int k = 0; k < thetaCh->size(); ++k) {
                if (0 < thetaCh->at(k) && thetaCh->at(k) < 2000) {
                    h1->Fill(1000.0*thetaCh->at(k));
                    //h2->Fill(tdcl->at(j).at(k));
                }
                //if (!dchit && 330 < tdcl->at(j).at(k) && tdcl->at(j).at(k) < 370) {
                //    h3->Fill(tdcl->at(j).at(k));
                //}
            }  // k hit loop of mppc j ch
        }
    }  // i event loop

    TF1* gaus_theta = new TF1("gaus_theta","gaus", 160, 240);
    //TF1* gaus_theta = new TF1("gaus_theta","gaus(0)+pol0(3)", 160, 240);
    //gaus_theta->SetParameters(14000, 201, 4.7, 2000);
    //gaus_theta->SetParameters(14000, 201, 4.7, 0, 0.01);
    gaus_theta->SetNpx(500);
    h1->Fit(gaus_theta, "","", 195, 212);

    //TF1* poli1 = new TF1("poli1", "pol0", 160, 240);
    //poli1->SetParameter(0,gaus_theta->GetParameter(3));//, gaus_theta->GetParameter(4));
    //poli1->SetLineColor(kViolet);

    TCanvas* cv = new TCanvas("cv", "cv", 800, 600);
    cv->SetGrid(1);
    cv->SetTicks(1);
    cv->SetLeftMargin(0.14);
    cv->SetRightMargin(0.1);
    cv->SetTopMargin(0.08);
    cv->SetBottomMargin(0.16);

    h1->Draw("");
    //poli1->Draw("same");
    //h2->Draw("same");
    //h3->Draw("same");

    //cv->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231204/thetaCh-1mm_run%d.pdf", run_no));
    //cv->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231129/thetaCh_run%d_mirror+1mm.pdf", run_no));
}