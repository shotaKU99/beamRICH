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

void tofRange_vs_nump3(std::string inputfilename){

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

    TFile* fin = new TFile(inputfilename.c_str());
    //TFile* fin = new TFile(Form("/mnt/c/Users/shota/ELPHroot/data/run%04d.root", run_no));
    TTree* tree = (TTree*)fin->Get("tree");

    int run_no = -1;
    std::string find_word = "run";
    std::string::size_type word_pos = inputfilename.rfind(find_word);
    if(word_pos!=std::string::npos){
        run_no = std::stoi(inputfilename.substr(word_pos+find_word.length(), 3));
    }

    std::vector<std::vector<double>>* tof = 0;
    std::vector<std::vector<double>>* photonNum = 0;
    //std::vector<std::vector<double>>* width = 0;
    //std::vector<std::vector<double>>* width_cut = 0;


    tree->SetBranchAddress("tof", &tof);
    tree->SetBranchAddress("photonNum", &photonNum);
    //tree->SetBranchAddress("width", &width);
    //tree->SetBranchAddress("width_cut", &width_cut);
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("tof", 1);
    tree->SetBranchStatus("photonNum", 1);
    //tree->SetBranchStatus("width", 1);
    //tree->SetBranchStatus("width_cut", 1);

    TH1F* htof = new TH1F("htof", "TOF;TOF [ns];counts", 200,-25,25);
    htof->GetXaxis()->SetMaxDigits(3);
    htof->GetYaxis()->SetMaxDigits(3);
    htof->SetStats(0);

    int ncut = 10;
    double tof_cutWidth = 1.0; //ns
    double tof_DCrange_low = 20;
    double tof_DCrange_high;// = tof_DCrange_low + ncut*tof_cutWidth;

    std::vector<TH1F*> hphotonDC(ncut);
    std::vector<TH1F*> hDC(ncut);
    std::vector<TLegend*> legend(ncut);

    std::vector<TF1*> fphotonDC(ncut);
    std::vector<TF1*> fDC(ncut);

    double min_bin = -0.5;
    double max_bin = 34.5;
    int nbin = (int)(max_bin-min_bin);
    for(int j=0;j<ncut;j++){
        hphotonDC.at(j) = new TH1F(Form("hphotonDC_%d",j), ";Multiplicity;events", nbin, min_bin, max_bin);
        hphotonDC.at(j)->GetXaxis()->SetMaxDigits(3);
        hphotonDC.at(j)->GetYaxis()->SetMaxDigits(3);
        hphotonDC.at(j)->SetLineColor(kBlack);
        hphotonDC.at(j)->SetFillColorAlpha(kRed, 0.3);
        hphotonDC.at(j)->SetStats(0);
        hDC.at(j) = new TH1F(Form("hDC_%d",j), ";Multiplicity;events", nbin, min_bin, max_bin);
        hDC.at(j)->GetXaxis()->SetMaxDigits(3);
        hDC.at(j)->GetYaxis()->SetMaxDigits(3);
        hDC.at(j)->SetLineColor(kBlack);
        hDC.at(j)->SetFillColorAlpha(kBlue, 0.3);
        hDC.at(j)->SetStats(0);
        fphotonDC.at(j) = new TF1(Form("fphotonDC_%d",j),"[0]*TMath::Poisson(x,[1])", min_bin, max_bin);
        fphotonDC.at(j)->SetLineColor(kBlue);
        fphotonDC.at(j)->SetLineWidth(2);
        fphotonDC.at(j)->SetNpx(1000);
        fDC.at(j) = new TF1(Form("fDC_%d",j),"[0]*TMath::Poisson(x,[1])", min_bin, max_bin);
        fDC.at(j)->SetLineColor(kRed);
        fDC.at(j)->SetLineWidth(2);
        fDC.at(j)->SetNpx(1000);
    }

    std::vector<int> nphotonDC(ncut);
    std::vector<int> nDC(ncut);

    double tof_value;//, tot_value, tot_cut_low;
    for(int i=0; i<tree->GetEntries(); i++){
    //for(int i=0; i<1; i++){
        tree->GetEntry(i);
        nphotonDC = std::vector<int>(ncut, 0);
        nDC = std::vector<int>(ncut, 0);

        //std::cout << "tof vector size = " << tof->size() << std::endl;
        for(int j=0; j<tof->size();j++){ // MPPC ch
            if(!tof->at(j).empty()){
                for(int k=0; k<tof->at(j).size(); k++){

                    tof_value = tof->at(j).at(k);

                    //std::cout<< "tof = " << tof_value << std::endl;
                    for(int l=0; l<ncut; l++){
                        if(std::abs(tof_value)<(l+1)*tof_cutWidth){
                            for(int m=l;m<ncut;m++){
                                nphotonDC.at(m) += (int)photonNum->at(j).at(k);
                            }
                            break;
                        }
                        
                        tof_DCrange_high = tof_DCrange_low+2*(l+1)*tof_cutWidth;
                        if(tof_DCrange_low<tof_value && tof_value<tof_DCrange_high){
                            for(int m=l;m<ncut;m++){
                                nDC.at(m) += (int)photonNum->at(j).at(k);
                            }  
                            break;                     
                        }
                    } // l tof cut width loop
                } // k hit size loop of 1 MPPC ch
            }
        } // j MPPC ch (384) loop

        //std::cout << "fill data " << i << " event" << std::endl;  
        for(int nc=0; nc<ncut; nc++){
            hphotonDC.at(nc)->Fill(nphotonDC.at(nc));
            hDC.at(nc)->Fill(nDC.at(nc));
        }

    } // i event loop

    //-----------------//
    //-- Poisson fit --//
    //-----------------//

    
    std::vector<double> tof_cutrange(ncut);
    std::vector<double> fitresult_mean_phdc(ncut);
    std::vector<double> fitresult_mean_dc(ncut);
    std::vector<double> fitresult_mean_ph(ncut);
    std::vector<double> fitresult_mean_phdc_err(ncut);
    std::vector<double> fitresult_mean_dc_err(ncut);
    std::vector<double> fitresult_mean_ph_err(ncut);

    for(int i=0; i<ncut; i++){
        fphotonDC.at(i)->SetParameters(hphotonDC.at(i)->GetMaximum(), hphotonDC.at(i)->GetBinCenter(hphotonDC.at(i)->GetMaximumBin()));
        fDC.at(i)->SetParameters(hDC.at(i)->GetMaximum(), hDC.at(i)->GetBinCenter(hDC.at(i)->GetMaximumBin())+0.1);

        hphotonDC.at(i)->Fit(fphotonDC.at(i),"Q", "", 2.0, max_bin);
        hDC.at(i)->Fit(fDC.at(i),"Q");

        std::cout << "Chisquare/NDF" <<  i  << " = " << fphotonDC.at(i)->GetChisquare()/fphotonDC.at(i)->GetNDF() << std::endl;
        tof_cutrange.at(i) = 2*(i+1)*tof_cutWidth;

        fitresult_mean_phdc.at(i) = fphotonDC.at(i)->GetParameter(1);
        fitresult_mean_phdc_err.at(i) = fphotonDC.at(i)->GetParError(1);

        fitresult_mean_dc.at(i) = fDC.at(i)->GetParameter(1);
        fitresult_mean_dc_err.at(i) = fDC.at(i)->GetParError(1);

        fitresult_mean_ph.at(i) = fitresult_mean_phdc.at(i) - fitresult_mean_dc.at(i);
        fitresult_mean_ph_err.at(i) = std::sqrt(std::pow(fitresult_mean_phdc_err.at(i),2.0)+std::pow(fitresult_mean_dc_err.at(i),2.0));
    }

    TGraphErrors* gphdc = new TGraphErrors(ncut, &tof_cutrange[0], &fitresult_mean_phdc[0], 0, &fitresult_mean_phdc_err[0]);
    gphdc -> SetMarkerSize(1.5);
    gphdc -> SetMarkerColor(kRed);
    gphdc -> SetLineColor(kRed);
    gphdc -> SetMarkerStyle(8);    
    
    TGraphErrors* gdc = new TGraphErrors(ncut, &tof_cutrange[0], &fitresult_mean_dc[0], 0, &fitresult_mean_dc_err[0]);
    gdc -> SetMarkerSize(1.5);
    gdc -> SetMarkerColor(kBlue);
    gdc -> SetLineColor(kBlue);
    gdc -> SetMarkerStyle(8);    

    TGraphErrors* gph = new TGraphErrors(ncut, &tof_cutrange[0], &fitresult_mean_ph[0], 0, &fitresult_mean_ph_err[0]);
    gph -> SetMarkerSize(1.5);
    gph -> SetMarkerColor(kViolet);
    gph -> SetLineColor(kViolet);
    gph -> SetMarkerStyle(8);    

    TCanvas* cv = new TCanvas("cv", "Multiplicity Photon + DC", 1200, 800);
    cv->Divide(int(std::ceil(ncut/2.0)),2);
    SetCanvasGrid(cv, ncut);

    for(int i=0;i<ncut;i++){
        cv->cd(i+1);
        hphotonDC.at(i)->Draw();
    }

    TCanvas* cv1 = new TCanvas("cv1", "Multiplicity DC only", 1200, 800);
    cv1->Divide(int(std::ceil(ncut/2.0)),2);
    SetCanvasGrid(cv1, ncut);

    for(int i=0;i<ncut;i++){
        cv1->cd(i+1);
        hDC.at(i)->Draw();
    }

    TCanvas* cv2 = new TCanvas("cv2", "Multiplicity", 800,600);
    SetCanvasGrid(cv2, 1);
    cv2->SetLeftMargin(0.14);
    cv2->SetRightMargin(0.14);
    cv2->SetTopMargin(0.08);
    cv2->SetBottomMargin(0.16);
    for(int i=0; i<ncut; i++){
        hDC.at(i)->GetYaxis()->SetRangeUser(0, std::max(fDC.at(i)->GetMaximum()+200, hDC.at(i)->GetMaximum())+200);
        hDC.at(i)->Draw();
        hphotonDC.at(i)->Draw("same");
        legend.at(i) = new TLegend(0.54, 0.78, 0.85, 0.89);
        legend.at(i)->AddEntry( hphotonDC.at(i), "photon+darkcurrent" , "f") ; // AddEntry( pointer , "interpretation" , "option" )
        //legend.at(i)->AddEntry( htheta2bg, "Dark Current" , "f") ; // option は　"f"=box, "l"="L"=line, "p"=marker
        legend.at(i)->AddEntry(hDC.at(i), "darkcurrent", "f");
        legend.at(i)->SetFillColor(0);
        legend.at(i)->Draw();        
        
        //cv2->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231211/Multiplicity_2_run%d_%d.pdf", run_no, i));
    }

    TCanvas* cv3 = new TCanvas("cv3", "Multiplicity", 800,600);
    SetCanvasGrid(cv3, 1);
    cv3->SetLeftMargin(0.14);
    cv3->SetRightMargin(0.14);
    cv3->SetTopMargin(0.08);
    cv3->SetBottomMargin(0.16);    

    TH1I* htemp = new TH1I("htemp",";TOF cut range [ns]; multiplicity", 22, 0, 22);
    htemp->SetStats(0);
    htemp->GetYaxis()->SetRangeUser(0, 18);
    htemp->Draw();
    gphdc->Draw("same P");
    gph->Draw("same P");
    gdc->Draw("same P");

    TLegend *legend1 = new TLegend(0.54, 0.4, 0.85, 0.6) ; //（）の中は位置の指定（左下の x , y 、右上の x , y ）
    legend1->AddEntry( gphdc, "photon+darkcurrent" , "p") ; // AddEntry( pointer , "interpretation" , "option" )
    //legend->AddEntry( htheta2bg, "Dark Current" , "f") ; // option は "f"=box, "l"="L"=line, "p"=marker
    legend1->AddEntry(gdc, "darkcurrent", "p");
    legend1->AddEntry(gph, "photon", "p");
    legend1->SetFillColor(0);
    legend1->Draw();

    //cv3->Print(Form("/mnt/c/Users/shota/OneDrive/Desktop/nh/e50/beam-RICH/image/20231211/TGraph_photon_DC_run%d_3.pdf", run_no));


}