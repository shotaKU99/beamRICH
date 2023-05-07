#include <iostream>
#include "Vector3D.hh"
#include "ExpConfig.hh"
#include "ReadYAML.hh"
#include "Ch2Pos.hh"

#include <vector>
#include <algorithm>
#include <time.h>

#include "Math/Minimizer.h"
//#include "Math/Factory.h"
//#include "Math/Functor.h"
#include "TRandom2.h"
#include "TError.h"
#include "TFile.h"
#include "TTree.h"
#include "Rtypes.h"
#include "TF2.h"
//#include "TCanvas.h"

/*
#include "Math/Vector3D.h"
#include "Math/Vector3Dfwd.h"

ROOT::Math::XYZVector Reflect(ROOT::Math::XYZVector incident, ROOT::Math::XYZVector nperp){
    ROOT::Math::XYZVector inci_norm = incident.Unit();
    ROOT::Math::XYZVector nperp_norm = nperp.Unit();

    return (inci_norm - 2.0*inci_norm.Dot(nperp_norm)*nperp_norm);
}

ROOT::Math::XYZVector Refract(ROOT::Math::XYZVector incident, ROOT::Math::XYZVector nperp, double n1, double n2){
    ROOT::Math::XYZVector inci_norm = incident.Unit();
    ROOT::Math::XYZVector nperp_norm = nperp.Unit();
    double costheta = inci_norm.Dot(nperp_norm);
    double k = 1.0 - std::pow(n1/n2, 2.0)*(1.0-std::pow(costheta, 2.0));
    if(k<0.0){
        return ROOT::Math::XYZVector(0.0,0.0,1.0);
    }else{
        return n1*inci_norm/n2 + (std::sqrt(k)-n1*costheta/n2)*nperp_norm;
    }
}
*/


double distance_func(double* xx, double* par){
    double x = xx[0];
    double y = xx[1];
    double z = 15.0;
    
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

    /*
    ROOT::Math::XYZVector rP(x, y, z);
    ROOT::Math::XYZVector rE(xE, yE, zE);
    ROOT::Math::XYZVector n_perp(nx, ny, nz);
    ROOT::Math::XYZVector r0Mirror(x0Mirror, y0Mirror, z0Mirror);
    ROOT::Math::XYZVector rD(xD, yD, zD);

    // Refraction
    ROOT::Math::XYZVector d_refraction = Refract(rP-rE, n_perp, n_Aero, n_Air);

    // Cross-Point of refraction light and Spherical Mirror
    ROOT::Math::XYZVector rPM0 = rP - r0Mirror;
    double t_c = -1.0*rPM0.Dot(d_refraction)
                    + std::sqrt(
                            Radius*Radius-rPM0.Mag2() 
                            + std::pow(rPM0.Dot(d_refraction), 2.0)
                        );
    ROOT::Math::XYZVector rM = rP + t_c*d_refraction;

    // Reflect on Spherical Mirror
    ROOT::Math::XYZVector d_reflect = Reflect(d_refraction, rM-r0Mirror);

    // Distance between detection point and reflect light
    ROOT::Math::XYZVector rDM = rD - rM;
    ROOT::Math::XYZVector cross_tmp = rDM.Cross(d_reflect);

    return cross_tmp.R()/d_reflect.R();
    */
    
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

    /*std::string configfile = "./info.yml"; 
    std::string settingname = "geant4_sim2"; 
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
    */
    std::string configfile = "./info.yml"; 
    std::string settingname = "geant4_sim1";
    ExpConfig expconfig = ExpConfig(configfile);
    //std::vector<std::string> keylist = expconfig.GetKeyList();    
    auto expinfo = expconfig.GetConstruction(settingname);


    std::string inputfilename = "sample.root";
    if(argc==1){
        std::cout << "Using default inputfile name :" << inputfilename << std::endl;
    }
    else if(argc==2){
        inputfilename = argv[1];
    }
    else if(argc>2){
        std::cout << "=========================================" << std::endl;
        std::cout << "Too Many arguments !!" << std::endl;
        std::cout << "Use default inputfile name : " << inputfilename << std::endl;
        std::cout << "=========================================" << std::endl;
    }
    
    ReadYAML readyaml = ReadYAML(inputfilename + ".yaml");
    auto configuration = readyaml.GetConfigration();

    TFile* inputFile = TFile::Open(inputfilename.c_str());
    TTree* inputTree = (TTree*)inputFile -> Get("tree");

    std::vector<Double_t>* xVD1 = 0;
    std::vector<Double_t>* yVD1 = 0;
    std::vector<Double_t>* zVD1 = 0;
    std::vector<Double_t>* xVD2 = 0;
    std::vector<Double_t>* yVD2 = 0;
    std::vector<Double_t>* zVD2 = 0;
    std::vector<Double_t>* xVD3 = 0;
    std::vector<Double_t>* yVD3 = 0;
    std::vector<Double_t>* zVD3 = 0;
    std::vector<Double_t>* qeffVD3 = 0;
    std::vector<Double_t>* hitChVD3 = 0;
    std::vector<Double_t>* parentIdVD3 = 0;

    inputTree -> SetBranchAddress("xVD1", &xVD1);
    inputTree -> SetBranchAddress("yVD1", &yVD1);
    inputTree -> SetBranchAddress("zVD1", &zVD1);
    inputTree -> SetBranchAddress("xVD2", &xVD2);
    inputTree -> SetBranchAddress("yVD2", &yVD2);
    inputTree -> SetBranchAddress("zVD2", &zVD2);
    inputTree -> SetBranchAddress("xVD3", &xVD3);
    inputTree -> SetBranchAddress("yVD3", &yVD3);
    inputTree -> SetBranchAddress("zVD3", &zVD3);
    inputTree -> SetBranchAddress("qeffVD3", &qeffVD3);
    inputTree -> SetBranchAddress("hitChVD3", &hitChVD3);
    inputTree -> SetBranchAddress("parentIdVD3", &parentIdVD3);

    /*    
    double xE = 0.0; //configuration.emit.x;
    double yE = 0.0; //configuration.emit.y;
    double zE = configuration.emit.z;
    double nx = configuration.n_perp.x;
    double ny = configuration.n_perp.y;
    double nz = configuration.n_perp.z;
    double x0Mirror = configuration.mirror.x0;
    double y0Mirror = configuration.mirror.y0;
    double z0Mirror = configuration.mirror.z0;
    double Radius = configuration.mirror.Radius;
    double n_Aero = configuration.rindex.Aero;
    double n_Air = configuration.rindex.Air;
    //double xD = -29.7523;
    //double yD = 1.18987;
    //double zD = 0.12986;
    double xD = 0.0;
    double yD = 0.0;
    double zD = 0.0;
    */
    
    //std::string ROOTFilename = configuration.setting.output_file;
    //int nThread = configuration.setting.numberOfThread;

    int ParticleFlag = configuration.beam.particle;

    double MomSpreadType = configuration.beam.MomSpreadType;
    double BeamMomentum = configuration.beam.momentum.center;
    double BeamP_w = configuration.beam.momentum.width;
    double BeamTheta = configuration.beam.xprime.center;
    double BeamTheta_w = configuration.beam.xprime.width;
    double BeamPhi = configuration.beam.yprime.center;
    double BeamPhi_w = configuration.beam.yprime.width;
    
    double PosSpreadType = configuration.beam.PosSpreadType;
    double beamoffx = configuration.beam.x.center;
    double rasterx = configuration.beam.x.width;
    double beamoffy = configuration.beam.y.center;
    double rastery = configuration.beam.y.width;
    double beamoffz = configuration.beam.z.center;
    double rasterz = configuration.beam.z.width;

    double detec_sizeX = configuration.detector.size.x;
    double detec_sizeY = configuration.detector.size.y;
    double detec_sizeZ = configuration.detector.size.z;

    double detec_centerX = configuration.detector.center.x;
    double detec_centerY = configuration.detector.center.y;
    double detec_centerZ = configuration.detector.center.z;

    double detec_rotangleX = configuration.detector.rotangle.x;
    double detec_rotangleY = configuration.detector.rotangle.y;
    double detec_rotangleZ = configuration.detector.rotangle.z;

    double mirror_rotangleX = configuration.mirror.rotangle.x;
    double mirror_rotangleY = configuration.mirror.rotangle.y;
    double mirror_rotangleZ = configuration.mirror.rotangle.z;
    double mirror_radius = configuration.mirror.radius;

    double radiator_centerX = configuration.radiator.center.x;
    double radiator_centerY = configuration.radiator.center.y;
    double radiator_centerZ = configuration.radiator.center.z;
    double radiator_sizeX = configuration.radiator.size.x;
    double radiator_sizeY = configuration.radiator.size.y;
    double radiator_sizeZ = configuration.radiator.size.z;

    int NumMPPC = readyaml.GetNumMPPC();
    std::vector<double> pos_X_mppc = readyaml.GetMppcPosX();
    std::vector<double> pos_Y_mppc = readyaml.GetMppcPosY();
    std::vector<double> pos_Z_mppc = readyaml.GetMppcPosZ();
    std::vector<double> rot_Z_mppc = readyaml.GetMppcZRot();

    int NumMppc = pos_X_mppc.size();
    std::cout << "Number of MPPC = " << NumMppc << std::endl;
    
    std::cout << "test 1" << std::endl;
    Ch2Pos channel2position = Ch2Pos(readyaml);
    std::cout << "Vector3D (100,0) = " << channel2position.GetChPos3D(100,0).x << std::endl;
    std::cout << "test 2" << std::endl;
    
    double xmin = -5.0;
    double xmax = 5.0;
    double ymin = -5.0;
    double ymax = 5.0;

    TF2* func = new TF2("func", distance_func, xmin, xmax, ymin, ymax, 15, 2);
    //func->SetNpx(100);
    //func->SetNpy(100);

    double xE = 0.0; //configuration.emit.x;
    double yE = 0.0; //configuration.emit.y;
    double zE = radiator_centerZ; //expinfo.emit.z;
    double nx = expinfo.n_perp.x;   
    double ny = expinfo.n_perp.y;
    double nz = expinfo.n_perp.z;
    double x0Mirror = mirror_radius*std::sin(mirror_rotangleY*TMath::Pi()/180.0);//expinfo.mirror.x0;
    double y0Mirror = 0.0;//expinfo.mirror.y0;
    double z0Mirror = mirror_radius*(1.0-std::cos(mirror_rotangleY*TMath::Pi()/180.0));//expinfo.mirror.z0;
    double Radius = mirror_radius;//expinfo.mirror.Radius;
    double n_Aero = expinfo.rindex.Aero;
    double n_Air = expinfo.rindex.Air;
    //double xD = -29.7523;
    //double yD = 1.18987;
    //double zD = 0.12986;
    double xD = 0.0;
    double yD = 0.0;
    double zD = 0.0;

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

    const char* outFile_name = expinfo.file.output_file.c_str();
    TFile* outFile = new TFile(outFile_name, "recreate");
    TTree* outTree = new TTree("tree", "tree");

    int nhit = 0;
    std::vector<Double_t>* qeff = 0;
    std::vector<Double_t>* radius = 0;
    std::vector<Double_t>* distanceO2D = 0;
    std::vector<Double_t>* xP = 0;
    std::vector<Double_t>* yP = 0;
    std::vector<Double_t>* minval = 0;
    std::vector<Double_t>* thetaCh = 0;
    outTree -> Branch("qeff", &qeff);
    outTree -> Branch("radius", &radius);
    outTree -> Branch("distanceO2D", &distanceO2D);
    outTree -> Branch("xP", &xP);
    outTree -> Branch("yP", &yP);
    outTree -> Branch("minval", &minval);
    outTree -> Branch("thetaCh", &thetaCh);


    double x_minimum = 0.0;
    double y_minimum = 0.0;
    int printVal = 0;
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-2);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);

    int nhit_total = 0;
    int nhit_overThreshold = 0;
    double minval_threshold = 0.001;
    double minimum_value = 0.0;

    clock_t start_time = clock();

    // Variables for Tracker 1,2
    double x1, y1, z1;
    double x2, y2, z2;
    double s1;

    Vector3D detec2center;

    for(Int_t j=0; j<inputTree->GetEntries(); j++){
    //for(Int_t j=0; j<1; j++){
        inputTree -> GetEntry(j);
        nhit = xVD3 -> size();
        nhit_total += nhit;

        qeff -> assign(nhit, -2222.0);
        radius -> assign(nhit, -2222.0);
        distanceO2D -> assign(nhit, -2222.0);
        xP -> assign(nhit, -2222.0);
        yP -> assign(nhit, -2222.0);
        minval -> assign(nhit, -2222.0);
        thetaCh -> assign(nhit, -2222.0);
        if(j%100==0) std::cout << "Vector Size (of " << j << " Event) = " << nhit << std::endl;

        // Tracking
        if(xVD1->size()>0 && xVD2->size()>0){
            x1 = xVD1->at(0);
            y1 = yVD1->at(0);
            z1 = zVD1->at(0);
            x2 = xVD2->at(0);
            y2 = yVD2->at(0);
            z2 = zVD2->at(0);
    
            s1 = (zE-z1)/(z2-z1);
    
            xE = x1 + s1*(x2-x1);
            yE = y1 + s1*(y2-y1);
            
            func -> SetParameter(0, xE);
            func -> SetParameter(1, yE);
    
            for(Int_t i=0; i<nhit; i++){
            //for(Int_t i=0; i<10; i++){
    
                qeff->at(i) = qeffVD3->at(i);
                xD = xVD3->at(i);
                yD = yVD3->at(i);
                zD = zVD3->at(i);
    
                func -> SetParameter(12, xD);
                func -> SetParameter(13, yD);
                func -> SetParameter(14, zD);           
                
                minimum_value = func->GetMinimumXY(x_minimum, y_minimum);
                minval->at(i) = minimum_value;
                if(minimum_value>minval_threshold){
                    nhit_overThreshold++;
                }
    
    
                xP->at(i) = x_minimum;
                yP->at(i) = y_minimum;
                //minval->at(i) = func->Eval(x_minimum, y_minimum);
                
                radius->at(i) = std::sqrt(std::pow(xD-detec_centerX,2.0) + std::pow(yD-detec_centerY,2.0));
                
                detec2center = Vector3D(xD,yD,zD) - channel2position.GetChPos3D(parentIdVD3->at(i), hitChVD3->at(i));
                distanceO2D->at(i) = detec2center.norm();

                if(printVal==1){
                    std::cout << "xE = " << func->GetParameter(0) << " +/- " << func->GetParError(0) << std::endl;
                    std::cout << "yE = " << func->GetParameter(1) << " +/- " << func->GetParError(1) << std::endl;
                    std::cout << "zE = " << func->GetParameter(2) << " +/- " << func->GetParError(2) << std::endl;
                    std::cout << "xD = " << func->GetParameter(12) << std::endl;
                    std::cout << "yD = " << func->GetParameter(13) << std::endl;
                    std::cout << "zD = " << func->GetParameter(14) << std::endl;
                    std::cout << "xD = " << channel2position.GetChPosX(parentIdVD3->at(i), hitChVD3->at(i)) << std::endl;
                    std::cout << "yD = " << channel2position.GetChPosY(parentIdVD3->at(i), hitChVD3->at(i)) << std::endl;
                    std::cout << "zD = " << channel2position.GetChPosZ(parentIdVD3->at(i), hitChVD3->at(i)) << std::endl;
                    std::cout << "========================" << std::endl;
                }
                
                //std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << min->MinValue()  << std::endl;
                Vector3D emit(x_minimum-xE, y_minimum-yE, zE+radiator_sizeZ/2.0 - zE);
                Vector3D nbeam(x2-x1, y2-y1, z2-z1);
                thetaCh->at(i) = std::acos(Vector3D::dot(emit, nbeam.normalized())/emit.norm());
    
            }
            outTree -> Fill();
        }
    }
    clock_t end_time = clock();

    outTree -> Write();
    outFile -> Close();

    std::cout << "------------ Analysis Info ------------" << std::endl;

    std::cout << "Input Root File name : " << inputfilename << std::endl;
    //std::cout << "Config File name : " << configfile << std::endl;
    //std::cout << "Setting name : " << settingname << std::endl;
    std::cout << "Output Root File name :" << outFile_name << std::endl;

    std::cout << "Total Event Size : " << inputTree->GetEntries() << std::endl;
    std::cout << "Total Hit Number : " << nhit_total << std::endl;
    std::cout << "Hit Number over threshold (= " << minval_threshold << ") : " << nhit_overThreshold << std::endl;
    std::cout << "Ratio of over threshold to total : " << 100.0*nhit_overThreshold/nhit_total << " %" << std::endl;
    std::cout << "Calculation time : " << (double)(end_time-start_time)/CLOCKS_PER_SEC << " sec" << std::endl;

    std::cout << "------------ Analysis Info ------------" << std::endl;
    std::cout << "Finish !!" << std::endl;

    return 0;    
}