#include <time.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "Vector3D.hh"
// #include "ExpConfig.hh"
#include "Ch2Pos.hh"
#include "ReadYAML.hh"
#include "Rtypes.h"
#include "TApplication.h"
#include "TCanvas.h"
#include "TError.h"
#include "TFile.h"
#include "TH2.h"
#include "TH2Poly.h"
#include "TRandom2.h"
#include "TRootCanvas.h"
#include "TTree.h"

int main(int argc, char** argv) {
    /*std::string configfile = "./info.yml";
    std::string settingname = "geant4_sim2";
    if(argc==1){

    }else if(argc==2){
        settingname = argv[1];
    }else{
        std::cout << "=========================================" << std::endl;
        std::cout << "Too Many arguments !!" << std::endl;
        std::cout << "Use default cofigfile name & setting name : " << configfile << ", " <<
    settingname << std::endl; std::cout << "=========================================" << std::endl;
    }

    ExpConfig expconfig = ExpConfig(configfile);

    std::vector<std::string> keylist = expconfig.GetKeyList();

    auto construction_info = expconfig.GetConstruction(settingname);
    const char* inputfile_name = construction_info.file.input_file.c_str();
    */

    std::string inputfilename = "sample.root";
    std::string setting_name = "geant4_sim1";
    int calculation_event = -1;
    if (argc == 1) {
        std::cout << "Using default inputfile name :" << inputfilename << std::endl;
    } else if (argc == 2) {
        inputfilename = argv[1];
    } else if (argc == 3) {
        inputfilename = argv[1];
        setting_name = argv[2];
    } else if (argc == 4) {
        inputfilename = argv[1];
        setting_name = argv[2];
        calculation_event = atoi(argv[3]);
    } else if (argc > 4) {
        std::cout << "=========================================" << std::endl;
        std::cout << "Too Many arguments !!" << std::endl;
        std::cout << "Use default inputfile name : " << inputfilename << std::endl;
        std::cout << "=========================================" << std::endl;
    }

    ReadYAML readyaml = ReadYAML("mppc.yaml", setting_name);
    readyaml.GetConfigration();

    TFile* inputFile = TFile::Open(inputfilename.c_str());
    TTree* inputTree = (TTree*)inputFile->Get("tree");

    // std::vector<Double_t>* xVD1b = 0;
    // std::vector<Double_t>* yVD1b = 0;
    // std::vector<Double_t>* zVD1b = 0;
    // std::vector<Double_t>* xVD4b = 0;
    // std::vector<Double_t>* yVD4b = 0;
    // std::vector<Double_t>* zVD4b = 0;
    std::vector<Double_t>* wlenVD3 = 0;
    std::vector<Double_t>* xVD3 = 0;
    std::vector<Double_t>* yVD3 = 0;
    std::vector<Double_t>* zVD3 = 0;
    std::vector<Double_t>* qeffVD3 = 0;
    std::vector<Double_t>* hitChVD3 = 0;
    std::vector<Double_t>* parentIdVD3 = 0;
    std::vector<Double_t>* nhitchAllVD3 = 0;
    // std::vector<Double_t>* hitChXEdepVD4 = 0;
    // std::vector<Double_t>* hitChYEdepVD4 = 0;

    // inputTree -> SetBranchAddress("xVD1b", &xVD1b);
    // inputTree -> SetBranchAddress("yVD1b", &yVD1b);
    // inputTree -> SetBranchAddress("zVD1b", &zVD1b);
    // inputTree -> SetBranchAddress("xVD4b", &xVD4b);
    // inputTree -> SetBranchAddress("yVD4b", &yVD4b);
    // inputTree -> SetBranchAddress("zVD4b", &zVD4b);
    inputTree->SetBranchAddress("wlenVD3", &wlenVD3);
    inputTree->SetBranchAddress("xVD3", &xVD3);
    inputTree->SetBranchAddress("yVD3", &yVD3);
    inputTree->SetBranchAddress("zVD3", &zVD3);
    inputTree->SetBranchAddress("qeffVD3", &qeffVD3);
    inputTree->SetBranchAddress("hitChVD3", &hitChVD3);
    inputTree->SetBranchAddress("parentIdVD3", &parentIdVD3);
    inputTree->SetBranchAddress("nhitchAllVD3", &nhitchAllVD3);
    // inputTree -> SetBranchAddress("hitChXEdepVD4", &hitChXEdepVD4);
    // inputTree -> SetBranchAddress("hitChYEdepVD4", &hitChYEdepVD4);

    // std::cout << "test 1" << std::endl;
    Ch2Pos channel2position = Ch2Pos(readyaml);
    int nhit, nhitper1Ch, nhit_total;
    int num_channel;
    int nloop_event = 0;
    int hitChX, hitChY;
    int nsize0 = 0;
    Vector3D detectionPoint;

    clock_t start_time = clock();
    nhit_total = 0;
    if (calculation_event == -1) {
        nloop_event = inputTree->GetEntries();
    } else if (calculation_event < inputTree->GetEntries()) {
        nloop_event = calculation_event;
    } else {
        nloop_event = inputTree->GetEntries();
    }

    //---------------------
    // Define histogram (TH2Poly)
    //---------------------
    // TH2F* hhitmap = new TH2F("hhitmap", "hitmap;x [cm];y [cm]", 2000, -7, 13, 300, -15, 15);
    TH2Poly* hhitmap2 = new TH2Poly("hhitmap2", "Hit Map; x [cm];y [cm]", -4, 6, -10, 10);
    double size_mppc1Ch = 0.3;  // cm
    int NumCh_1mppc = 16;
    int NumMPPC = readyaml.GetNumMPPC();
    int NumChannel = NumMPPC * NumCh_1mppc;
    int parentIdOffset = readyaml.GetParentIdOffset();
    int hitChOffset = readyaml.GetHitChOffset();
    std::cout << "Number of channel " << NumChannel << std::endl;
    double xpos_center, ypos_center;
    for (int i = 0; i < NumChannel; i++) {
        xpos_center = channel2position.GetChPosX(i / NumChannel + parentIdOffset,
                                                 i % NumChannel + hitChOffset);
        ypos_center = channel2position.GetChPosY(i / NumChannel + parentIdOffset,
                                                 i % NumChannel + hitChOffset);
        hhitmap2->AddBin(xpos_center - 0.5 * size_mppc1Ch, ypos_center - 0.5 * size_mppc1Ch,
                         xpos_center + 0.5 * size_mppc1Ch, ypos_center + 0.5 * size_mppc1Ch);
    }

    //------------------
    // Fill data to histogram
    //------------------
    // std::cout << "test2 " << std::endl;
    // for (Int_t j = 0; j < nloop_event; j++) {
    //    // for(Int_t j=0; j<1; j++){
    //    inputTree->GetEntry(j);
    //    nhit = qeffVD3->size();
    //    for (Int_t i = 0; i < nhit; i++) {
    //        if (qeffVD3->at(i) == 1) {
    //            if (wlenVD3->at(i) < 2100) {
    //                detectionPoint =
    //                    channel2position.GetChPos3D(parentIdVD3->at(i), hitChVD3->at(i));
    //                nhit_total++;
    //                //hhitmap->Fill(detectionPoint.x, detectionPoint.y);
    //                hhitmap2->Fill(detectionPoint.x, detectionPoint.y);
    //            }
    //        }
    //    }
    //}
    for (Int_t j = 0; j < nloop_event; j++) {
        inputTree->GetEntry(j);
        num_channel = nhitchAllVD3->size();
        for (Int_t i = 0; i < num_channel; i++) {
            nhitper1Ch = nhitchAllVD3->at(i);
            for (int k = 0; k < nhitper1Ch; k++) {
                detectionPoint = channel2position.GetChPos3D(i / NumChannel + parentIdOffset,
                                                             i % NumChannel + hitChOffset);
                nhit_total++;
                // hhitmap->Fill(detectionPoint.x, detectionPoint.y);
                hhitmap2->Fill(detectionPoint.x, detectionPoint.y);
            }
        }
    }

    std::cout << "Underflow = " << hhitmap2->GetBinContent(0) << std::endl;
    std::cout << "Overflow = " << hhitmap2->GetBinContent(385) << std::endl;

    clock_t end_time = clock();

    std::cout << "------------ Analysis Info ------------" << std::endl;

    std::cout << "Input Root File name : " << inputfilename << std::endl;
    // std::cout << "Config File name : " << configfile << std::endl;
    // std::cout << "Setting name : " << settingname << std::endl;
    std::cout << "Total Event Size : " << inputTree->GetEntries() << std::endl;
    std::cout << "Total Hit Number : " << nhit_total << std::endl;
    std::cout << "Calculation time : " << (double)(end_time - start_time) / CLOCKS_PER_SEC << " sec"
              << std::endl;

    std::cout << "------------ Analysis Info ------------" << std::endl;
    std::cout << "Finish !!" << std::endl;

    TApplication app("app", &argc, argv);
    TCanvas* cv = new TCanvas("cv", "cv", 500, 1000);
    cv->SetGrid(1);
    cv->SetTicks(1);
    TRootCanvas* rc = (TRootCanvas*)cv->GetCanvasImp();
    rc->Connect("CloseWindow()", "TApplication", gApplication, "Terminate()");

    // hhitmap -> Draw("colz");
    hhitmap2->Draw("colz 0");
    cv->SetRealAspectRatio(2);
    cv->Draw();

    app.Run();

    return 0;
}
