#include <iostream>
#include "Vector3D.hh"
#include "ExpConfig.hh"

#include <vector>
#include <algorithm>

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "Rtypes.h"
#include "TF2.h"
#include "TCanvas.h"

double distance_func(double* xx, double* par){
    double x = xx[0];
    double y = xx[1];
    double z = 45.0;
    
    // Emission Point
    double xE = par[0];
    double yE = par[1];
    double zE = par[2];

    // Perpendicular of back of Radiator
    double nx = par[3];
    double ny = par[4];
    double nz = par[5];

    // Center & Radius of Spherical Mirror
    double x0Mirror = par[6];
    double y0Mirror = par[7];
    double z0Mirror = par[8];
    double Radius = par[9];

    // Refractive Index
    double n_Aero = par[10];
    double n_Air = par[11];

    // Detection Point
    double xD = par[12];
    double yD = par[13];
    double zD = par[14];


    Vector3D rP(x, y, z);
    Vector3D rE(xE, yE, zE);
    Vector3D n_perp(nx, ny, nz);
    Vector3D r0Mirror(x0Mirror, y0Mirror, z0Mirror);
    Vector3D rD(xD, yD, zD);

    // Refraction
    Vector3D d_refraction = Vector3D::refraction(rP-rE, n_perp, n_Aero, n_Air);

    // Cross-Point of refraction light and Spherical Mirror
    Vector3D rPM0 = rP - r0Mirror;
    double t_c = -1.0*Vector3D::dot(rPM0, d_refraction) 
                    + std::sqrt(
                            Radius*Radius-rPM0.norm2() 
                            + std::pow(Vector3D::dot(rPM0, d_refraction), 2.0)
                        );
    Vector3D rM = rP + t_c*d_refraction;

    // Reflect on Spherical Mirror
    Vector3D d_reflect = Vector3D::reflect(d_refraction, rM-r0Mirror);

    // Distance between detection point and reflect light
    Vector3D cross_tmp = Vector3D::cross(rD-rM, d_reflect);

    return cross_tmp.norm()/d_reflect.norm();
}

int main(int argc, char** argv){

    std::string configfile = "./info.yml"; 
    std::string settingname = "geant4_sim1"; 
    if(argc==1){

    }else if(argc==2){
        settingname = argv[1];
    }else{
        std::cout << "=========================================" << std::endl;
        std::cout << "Too Many arguments !!" << std::endl;
        std::cout << "Use default cofigfile name & setting name : " << configfile << ", " << settingname << std::endl;
        std::cout << "=========================================" << std::endl;
    }

    ExpConfig expconfig = ExpConfig(configfile);

    std::vector<std::string> keylist = expconfig.GetKeyList();
    
    auto construction_info = expconfig.GetConstruction(settingname);

    const char* inputfile_name = construction_info.file.input_file.c_str();
    TFile* inputFile = TFile::Open(inputfile_name);
    TTree* inputTree = (TTree*)inputFile -> Get("tree");

    std::vector<Double_t>* xVD3 = 0;
    std::vector<Double_t>* yVD3 = 0;
    std::vector<Double_t>* zVD3 = 0;

    inputTree -> SetBranchAddress("xVD3", &xVD3);
    inputTree -> SetBranchAddress("yVD3", &yVD3);
    inputTree -> SetBranchAddress("zVD3", &zVD3);
    
    double xE = construction_info.emit.x;
    double yE = construction_info.emit.y;
    double zE = construction_info.emit.z;
    double nx = construction_info.n_perp.x;
    double ny = construction_info.n_perp.y;
    double nz = construction_info.n_perp.z;
    double x0Mirror = construction_info.mirror.x0;
    double y0Mirror = construction_info.mirror.y0;
    double z0Mirror = construction_info.mirror.z0;
    double Radius = construction_info.mirror.Radius;
    double n_Aero = construction_info.rindex.Aero;
    double n_Air = construction_info.rindex.Air;
    //double xD = -29.7523;
    //double yD = 1.18987;
    //double zD = 0.12986;
    double xD = 0.0;
    double yD = 0.0;
    double zD = 0.0;

    double xmin = -5.0;
    double xmax = 5.0;
    double ymin = -5.0;
    double ymax = 5.0;

    TF2* func = new TF2("func", distance_func, xmin, xmax, ymin, ymax, 15, 2);
    //func->SetNpx(100);
    //func->SetNpy(100);

    func -> SetParameter(0, xE);
    func -> SetParameter(1, yE);
    func -> SetParameter(2, zE);
    func -> SetParameter(3, nx);
    func -> SetParameter(4, ny);
    func -> SetParameter(5, nz);
    func -> SetParameter(6, x0Mirror);
    func -> SetParameter(7, y0Mirror);
    func -> SetParameter(8, z0Mirror);
    func -> SetParameter(9, Radius);
    func -> SetParameter(10, n_Aero);
    func -> SetParameter(11, n_Air);
    func -> SetParameter(12, xD);
    func -> SetParameter(13, yD);
    func -> SetParameter(14, zD);

    const char* outFile_name = construction_info.file.output_file.c_str();
    TFile* outFile = new TFile(outFile_name, "recreate");
    TTree* outTree = new TTree("tree", "tree");

    int nhit = 0;
    std::vector<Double_t>* radius = 0;
    std::vector<Double_t>* xP = 0;
    std::vector<Double_t>* yP = 0;
    std::vector<Double_t>* minval = 0;
    std::vector<Double_t>* thetaCh = 0;
    outTree -> Branch("radius", &radius);
    outTree -> Branch("xP", &xP);
    outTree -> Branch("yP", &yP);
    outTree -> Branch("minval", &minval);
    outTree -> Branch("thetaCh", &thetaCh);


    double x_minimum = 0.0;
    double y_minimum = 0.0;
    int printVal = 0;
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-2);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);


    for(Int_t j=0; j<inputTree->GetEntries(); j++){
    //for(Int_t j=0; j<1; j++){
        inputTree -> GetEntry(j);
        nhit = xVD3 -> size();
        radius -> assign(nhit, -2222.0);
        xP -> assign(nhit, -2222.0);
        yP -> assign(nhit, -2222.0);
        minval -> assign(nhit, -2222.0);
        thetaCh -> assign(nhit, -2222.0);
        if(j%100==0) std::cout << "Vector Size (of " << j << " Event) = " << nhit << std::endl;


        for(Int_t i=0; i<nhit; i++){
        //for(Int_t i=0; i<10; i++){

            xD = xVD3->at(i);
            yD = yVD3->at(i);
            zD = zVD3->at(i);

            func -> SetParameter(12, xD);
            func -> SetParameter(13, yD);
            func -> SetParameter(14, zD);           

            minval->at(i) = func->GetMinimumXY(x_minimum, y_minimum);

            xP->at(i) = x_minimum;
            yP->at(i) = y_minimum;
            //minval->at(i) = func->Eval(x_minimum, y_minimum);
            
            radius->at(i) = std::sqrt(xD*xD + yD*yD);

            if(printVal==1){
                std::cout << "xE = " << func->GetParameter(0) << " +/- " << func->GetParError(0) << std::endl;
                std::cout << "yE = " << func->GetParameter(1) << " +/- " << func->GetParError(1) << std::endl;
                std::cout << "zE = " << func->GetParameter(2) << " +/- " << func->GetParError(2) << std::endl;
                std::cout << "xD = " << func->GetParameter(12) << std::endl;
                std::cout << "yD = " << func->GetParameter(13) << std::endl;
                std::cout << "zD = " << func->GetParameter(14) << std::endl;
                std::cout << "========================" << std::endl;
            }
            
            //std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << min->MinValue()  << std::endl;
            Vector3D emit(x_minimum-xE, y_minimum-yE, 45.0 - zE);
            thetaCh->at(i) = std::acos(Vector3D::dot(emit, Vector3D(0.0,0.0,1.0))/emit.norm());

        }
        outTree -> Fill();
    }

    outTree -> Write();
    outFile -> Close();

    std::cout << "---------------- Info -------------------" << std::endl;
    std::cout << "Input Root File name : " << inputfile_name << std::endl;
    std::cout << "Config File name : " << configfile << std::endl;
    std::cout << "Setting name : " << settingname << std::endl;
    std::cout << "Output Root File name :" << outFile_name << std::endl;
    std::cout << "---------------- Info -------------------" << std::endl;
    std::cout << "Finish !!" << std::endl;

    return 0;    
}