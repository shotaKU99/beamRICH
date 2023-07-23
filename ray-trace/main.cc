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


double distance_func(const double* xx){
    const double x = xx[0];
    const double y = xx[1];
    const double z = xx[2];
    
    // Emission Point
    const double xE = xx[3];
    const double yE = xx[4];
    const double zE = xx[5];

    // Perpendicular of back of Radiator
    const double nx = xx[6];
    const double ny = xx[7];
    const double nz = xx[8];

    // Center & Radius of Spherical Mirror
    const double x0Mirror = xx[9];
    const double y0Mirror = xx[10];
    const double z0Mirror = xx[11];
    const double Radius = xx[12];

    // Refractive Index
    const double n_Aero = xx[13];
    const double n_Air = xx[14];

    // Detection Point
    const double xD = xx[15];
    const double yD = xx[16];
    const double zD = xx[17];


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


    const char* minName="Minuit2";
    //const char* minName="GSLMultiMin";
    const char* algoName="";
    //const char* algoName="BFGS2";
    int randomSeed = -1;
    /*
    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    min->SetMaxFunctionCalls(10000000);
    min->SetMaxIterations(10000000);
    min->SetTolerance(1e-5);
    min->SetPrintLevel(0);

    ROOT::Math::Functor func(&distance_func, 18);
    
    double step[2] = {0.001,0.001};
    double variable[2] = {-4., 4.};

    min->SetFunction(func);
    min->SetVariable(0, "x", variable[0], step[0]);
    min->SetVariableLimits(0, -5.0, 5.0);
    

    min->SetVariable(1, "y", variable[1], step[1]);
    min->SetVariableLimits(1, -5.0, 5.0);
    
    double zP = 45.0;
    min->SetFixedVariable(2, "z", zP);

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
    */

/*    
    double xE = 0.0;
    double yE = 0.0;
    double zE = 40.0;
    double nx = 0.0;
    double ny = 0.0;
    double nz = 1.0;
    double x0Mirror = 0.0;
    double y0Mirror = 0.0;
    double z0Mirror = -150.0;
    double Radius = 300.0;
*/

/*    min->SetFixedVariable(3, "xE", xE);
    min->SetFixedVariable(4, "yE", yE);
    min->SetFixedVariable(5, "zE", zE);
    min->SetFixedVariable(6, "nx", nx);
    min->SetFixedVariable(7, "ny", ny);
    min->SetFixedVariable(8, "nz", nz);
    min->SetFixedVariable(9, "x0Mirror", x0Mirror);
    min->SetFixedVariable(10,"y0Mirror", y0Mirror);
    min->SetFixedVariable(11,"z0Mirror", z0Mirror);
    min->SetFixedVariable(12,"Radius", Radius);
    min->SetFixedVariable(13,"nAero", n_Aero);
    min->SetFixedVariable(14,"nAir", n_Air);
*/
    const char* outFile_name = construction_info.file.output_file.c_str();
    TFile* outFile = new TFile(outFile_name, "recreate");
    TTree* outTree = new TTree("tree", "tree");

    int nhit = 0;
    std::vector<Double_t>* radius = 0;
    std::vector<Int_t>* status = 0;
    std::vector<Double_t>* xP = 0;
    std::vector<Double_t>* yP = 0;
    std::vector<Double_t>* minval = 0;
    std::vector<Double_t>* edm = 0;
    std::vector<Double_t>* thetaCh = 0;
    outTree -> Branch("radius", &radius);
    outTree -> Branch("status", &status);
    outTree -> Branch("xP", &xP);
    outTree -> Branch("yP", &yP);
    outTree -> Branch("minval", &minval);
    outTree -> Branch("edm", &edm);
    outTree -> Branch("thetaCh", &thetaCh);

    double step[2] = {0.001,0.001};
    double variable[2] = {4., 4.};

    ROOT::Math::Minimizer* min = ROOT::Math::Factory::CreateMinimizer(minName, algoName);

    min->SetMaxFunctionCalls(10000000);
    min->SetMaxIterations(10000000);
    min->SetTolerance(1e-5);
    min->SetPrintLevel(0);
    //min->SetValidError(true);

    ROOT::Math::Functor func(&distance_func, 18);
    min->SetFunction(func);
    min->SetVariable(0, "x", variable[0], step[0]);
    //min->SetVariableInitialRange(0, -5.0, 5.0);
    min->SetVariableLimits(0, -5.0, 5.0);

    min->SetVariable(1, "y", variable[1], step[1]);
    //min->SetVariableInitialRange(1, -5.0, 5.0);
    min->SetVariableLimits(1, -5.0, 5.0);

    double zP = 45.0;
    min->SetFixedVariable(2, "z", zP);

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

    min->SetFixedVariable(3, "xE", xE);
    min->SetFixedVariable(4, "yE", yE);
    min->SetFixedVariable(5, "zE", zE);
    min->SetFixedVariable(6, "nx", nx);
    min->SetFixedVariable(7, "ny", ny);
    min->SetFixedVariable(8, "nz", nz);
    min->SetFixedVariable(9, "x0Mirror", x0Mirror);
    min->SetFixedVariable(10,"y0Mirror", y0Mirror);
    min->SetFixedVariable(11,"z0Mirror", z0Mirror);
    min->SetFixedVariable(12,"Radius", Radius);
    min->SetFixedVariable(13,"nAero", n_Aero);
    min->SetFixedVariable(14,"nAir", n_Air);

    double xD = 0.0;
    double yD = 0.0;
    double zD = 0.0;
    min->SetFixedVariable(15,"xD", xD);
    min->SetFixedVariable(16,"yD", yD);
    min->SetFixedVariable(17,"zD", zD);

    for(Int_t j=0; j<inputTree->GetEntries(); j++){
    //for(Int_t j=0; j<1; j++){

        inputTree -> GetEntry(j);
        nhit = xVD3 -> size();
        radius -> assign(nhit, -2222.0);
        status -> assign(nhit, -10.0);
        xP -> assign(nhit, -2222.0);
        yP -> assign(nhit, -2222.0);
        minval -> assign(nhit, -2222.0);
        edm -> assign(nhit, -2222.0);
        thetaCh -> assign(nhit, -2222.0);
        if(j%100==0) std::cout << "Vector Size (of " << j << " Event) = " << nhit << std::endl;


        for(Int_t i=0; i<nhit; i++){
        //for(Int_t i=0; i<2; i++){

            xD = xVD3->at(i);
            yD = yVD3->at(i);
            zD = zVD3->at(i);

            //min->SetVariableValue(0, 5.0);
            //min->SetVariableValue(1, 5.0);

            min->SetVariableValue(15, xD);
            min->SetVariableValue(16, yD);
            min->SetVariableValue(17, zD);

            min->Minimize();

            //min->SetVariable(0, "x", min->X()[0]+0.01, step[0]);
            //min->SetVariableLimits(0, -5.0, 5.0);
    
            //min->SetVariable(1, "y", min->X()[1]+0.01, step[1]);
            //min->SetVariableLimits(1, -5.0, 5.0);
            //min->SetLimitedVariable(0, "x", min->X()[0]+0.01, step[0], -5.0, 5.0);
            //min->SetLimitedVariable(1, "y", min->X()[1]+0.01, step[1], -5.0, 5.0);

            //min->Minimize();
            //min->Minimize();
            //min->Minimize();
            
            radius->at(i) = std::sqrt(xD*xD + yD*yD);

            status->at(i) = min->Status();    

            const double* xs = min->X();
            xP->at(i) = xs[0];
            yP->at(i) = xs[1];
            minval->at(i) = min->MinValue();
            edm->at(i) = min->Edm();
            //std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << min->MinValue()  << std::endl;
            Vector3D emit(xs[0]-xE, xs[1]-yE, zP - zE);
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