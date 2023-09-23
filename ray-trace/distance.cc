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
#include "TF1.h"
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
    double z = par[0];
    
    // Emission Point
    double xE = par[1];
    double yE = par[2];
    double zE = par[3];

    // Perpendicular of back of Radiator
    double nx = par[4];
    double ny = par[5];
    double nz = par[6];

    // Center & Radius of Spherical Mirror
    double x0Mirror = par[7];
    double y0Mirror = par[8];
    double z0Mirror = par[9];
    double Radius = par[10];

    // Refractive Index
    double n_Aero = par[11];
    double n_Air = par[12];

    // Detection Point
    double xD = par[13];
    double yD = par[14];
    double zD = par[15];

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

    
    //std::cout << "rP_norm = " << rP.normalized() << std::endl;
    //std::cout << "rE = " << rE << std::endl;
    //std::cout << "rE_norm = " << rE.normalized() << std::endl;
    //std::cout << "rE = " << rE << std::endl;
    //std::cout << "(rP-rE)_norm = " << (rP-rE).normalized() << std::endl;
    //double costheta =  std::abs(Vector3D::dot((rP-rE).normalized(), Vector3D(0,0,1.0)));
    //std::cout <<  "costheta = " << costheta << std::endl;
    //double k = 1.0 - std::pow(1.021 / 1.000273, 2.0) * (1.0 - std::pow(costheta, 2.0));
    //std::cout << "k = " << k << std::endl;
    
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

double refractive_index(double* wlen, double* par) {
    double wlength = wlen[0];  // nm
    double a0 = par[0];
    double wlength0 = par[1];  // nm

    return std::sqrt(1 + a0 * std::pow(wlength, 2.0) /
                             (std::pow(wlength, 2.0) - std::pow(wlength0, 2.0)));
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

    std::vector<Double_t>* xVD1b = 0;
    std::vector<Double_t>* yVD1b = 0;
    std::vector<Double_t>* zVD1b = 0;
    std::vector<Double_t>* xVD4b = 0;
    std::vector<Double_t>* yVD4b = 0;
    std::vector<Double_t>* zVD4b = 0;
    std::vector<Double_t>* wlenVD3 = 0;
    std::vector<Double_t>* xVD3 = 0;
    std::vector<Double_t>* yVD3 = 0;
    std::vector<Double_t>* zVD3 = 0;
    std::vector<Double_t>* qeffVD3 = 0;
    std::vector<Double_t>* hitChVD3 = 0;
    std::vector<Double_t>* parentIdVD3 = 0;
    std::vector<Double_t>* hitChXEdepVD4 = 0;
    std::vector<Double_t>* hitChYEdepVD4 = 0;

    inputTree -> SetBranchAddress("xVD1b", &xVD1b);
    inputTree -> SetBranchAddress("yVD1b", &yVD1b);
    inputTree -> SetBranchAddress("zVD1b", &zVD1b);
    inputTree -> SetBranchAddress("xVD4b", &xVD4b);
    inputTree -> SetBranchAddress("yVD4b", &yVD4b);
    inputTree -> SetBranchAddress("zVD4b", &zVD4b);
    inputTree -> SetBranchAddress("wlenVD3", &wlenVD3);
    inputTree -> SetBranchAddress("xVD3", &xVD3);
    inputTree -> SetBranchAddress("yVD3", &yVD3);
    inputTree -> SetBranchAddress("zVD3", &zVD3);
    inputTree -> SetBranchAddress("qeffVD3", &qeffVD3);
    inputTree -> SetBranchAddress("hitChVD3", &hitChVD3);
    inputTree -> SetBranchAddress("parentIdVD3", &parentIdVD3);
    inputTree -> SetBranchAddress("hitChXEdepVD4", &hitChXEdepVD4);
    inputTree -> SetBranchAddress("hitChYEdepVD4", &hitChYEdepVD4);

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

    double upstream_trig_centerX = configuration.uptrig.center.x;
    double upstream_trig_centerY = configuration.uptrig.center.y;
    double upstream_trig_centerZ = configuration.uptrig.center.z;

    double downstream_trig_centerX = configuration.downtrig.center.x;
    double downstream_trig_centerY = configuration.downtrig.center.y;
    double downstream_trig_centerZ = configuration.downtrig.center.z;


    int NumMPPC = readyaml.GetNumMPPC();
    std::vector<double> pos_X_mppc = readyaml.GetMppcPosX();
    std::vector<double> pos_Y_mppc = readyaml.GetMppcPosY();
    std::vector<double> pos_Z_mppc = readyaml.GetMppcPosZ();
    std::vector<double> rot_Z_mppc = readyaml.GetMppcZRot();

    int NumMppc = pos_X_mppc.size();
    //std::cout << "Number of MPPC = " << NumMppc << std::endl;
    
    //std::cout << "test 1" << std::endl;
    Ch2Pos channel2position = Ch2Pos(readyaml);
    //std::cout << "Vector3D (100,0) = " << channel2position.GetChPos3D(100,0).x << std::endl;
    //std::cout << "test 2" << std::endl;
    
    double xmin = -5.0;
    double xmax = 5.0;
    double ymin = -5.0;
    double ymax = 5.0;

    TF2* func = new TF2("func", distance_func, xmin, xmax, ymin, ymax, 16, 2);
    //func->SetNpx(100);
    //func->SetNpy(100);

    double zP = radiator_centerZ + radiator_sizeZ/2.0;
    double xE = 0.0; //configuration.emit.x;
    double yE = 0.0; //configuration.emit.y;
    double zE = radiator_centerZ + radiator_sizeZ/2.0 - 1.5; //expinfo.emit.z;
    //std::cout << radiator_centerZ + radiator_sizeZ/2.0 << std::endl;
    double nx = expinfo.n_perp.x;   
    double ny = expinfo.n_perp.y;
    double nz = expinfo.n_perp.z;


    Vector3D MirrorR(0.0, 0.0, mirror_radius); // Vector of center of rotation 
    Vector3D MirrorR2O(0.0,0.0, -mirror_radius); // Vector from center of rotation to Mirror center
    Vector3D rotMirrorR2O = MirrorR2O.rotateY(-mirror_rotangleY*TMath::Pi()/180.0).rotateX(mirror_rotangleX*TMath::Pi()/180.0);
    //double x0Mirror = (MirrorR + rotMirrorR2O).x;
    //double y0Mirror = (MirrorR + rotMirrorR2O).y;
    //double z0Mirror = (MirrorR + rotMirrorR2O).z;
    //std::cout << "Center of Mirror = " << (MirrorR+rotMirrorR2O) << std::endl;
    double x0Mirror = mirror_radius*std::sin(mirror_rotangleY*TMath::Pi()/180.0);//expinfo.mirror.x0;
    double y0Mirror = 0.0;//expinfo.mirror.y0;
    double z0Mirror = mirror_radius*(1.0-std::cos(mirror_rotangleY*TMath::Pi()/180.0));//expinfo.mirror.z0;
    std::cout << "Center of Mirror = (" << x0Mirror << ", " << y0Mirror << ", " << z0Mirror << ")"<< std::endl;
    
    double Radius = mirror_radius;//expinfo.mirror.Radius;
    double a0_Aero = expinfo.rindex.a0_Aero;
    double wlen0_Aero = expinfo.rindex.wlen0_Aero;
    double wlength_Aero = expinfo.rindex.wlength_Aero;
    double n_Air = expinfo.rindex.Air;

    TF1* refractive_index_func = new TF1("refractive_index_func", refractive_index, 200, 900, 2);
    refractive_index_func->SetParameters(a0_Aero, wlen0_Aero);
    double n_Aero = refractive_index_func->Eval(wlength_Aero);

    //double xD = -29.7523;
    //double yD = 1.18987;
    //double zD = 0.12986;
    double xD = 0.0;
    double yD = 0.0;
    double zD = 0.0;

    func -> SetParameter(0, zP);
    func -> SetParameter(1, xE);
    func -> SetParameter(2, yE);
    func -> SetParameter(3, zE);
    func -> SetParameter(4, nx);
    func -> SetParameter(5, ny);
    func -> SetParameter(6, nz);
    func -> SetParameter(7, x0Mirror);
    func -> SetParameter(8, y0Mirror);
    func -> SetParameter(9, z0Mirror);
    func -> SetParameter(10, Radius);
    func -> SetParameter(11, n_Aero);
    func -> SetParameter(12, n_Air);
    func -> SetParameter(13, xD);
    func -> SetParameter(14, yD);
    func -> SetParameter(15, zD);

    const char* outFile_name = expinfo.file.output_file.c_str();
    TFile* outFile = new TFile(outFile_name, "recreate");
    TTree* outTree = new TTree("tree", "tree");

    int nhit = 0;
    std::vector<Double_t>* wlen = 0;
    std::vector<Double_t>* qeff = 0;
    std::vector<Double_t>* radius = 0;
    std::vector<Double_t>* distanceO2D = 0;
    std::vector<Double_t>* xVD = 0;
    std::vector<Double_t>* yVD = 0;
    std::vector<Double_t>* xP = 0;
    std::vector<Double_t>* yP = 0;
    std::vector<Double_t>* minval = 0;
    std::vector<Double_t>* thetaCh = 0;
    std::vector<Double_t>* xP2 = 0;
    std::vector<Double_t>* yP2 = 0;
    std::vector<Double_t>* minval2 = 0;
    std::vector<Double_t>* thetaCh2 = 0;
    outTree -> Branch("wlen", &wlen);
    outTree -> Branch("qeff", &qeff);
    outTree -> Branch("radius", &radius);
    outTree -> Branch("distanceO2D", &distanceO2D);
    outTree -> Branch("xVD", &xVD);
    outTree -> Branch("yVD", &yVD);
    outTree -> Branch("xP", &xP);
    outTree -> Branch("yP", &yP);
    outTree -> Branch("minval", &minval);
    outTree -> Branch("thetaCh", &thetaCh);
    outTree -> Branch("xP2", &xP2);
    outTree -> Branch("yP2", &yP2);
    outTree -> Branch("minval2", &minval2);
    outTree -> Branch("thetaCh2", &thetaCh2);


    double x_minimum = 0.0;
    double y_minimum = 0.0;
    double x_minimum2 = 0.0;
    double y_minimum2 = 0.0;
    int printVal = 0;
    ROOT::Math::MinimizerOptions::SetDefaultTolerance(1e-2);
    ROOT::Math::MinimizerOptions::SetDefaultPrintLevel(0);

    int nhit_total = 0;
    int nhit_overThreshold = 0;
    int nhit_overThreshold2 = 0;
    double minval_threshold = 0.001;
    double minimum_value = 0.0;
    double minimum_value2 = 0.0;

    clock_t start_time = clock();


    double size_beamPM = 4.8; //cm
    double size_1ChbeamPM = 0.3; //cm
    // Variables for Tracker 1,2
    double x1, y1, z1;
    double x2, y2, z2;
    int hitChX, hitChY;
    double s1;
    int nsize0=0;
    Vector3D detec2center, detectionPoint;

    for(Int_t j=0; j<inputTree->GetEntries(); j++){
    //for(Int_t j=0; j<1; j++){
        inputTree -> GetEntry(j);
        nhit = xVD3 -> size();
        nhit_total += nhit;

        wlen -> assign(nhit, -2222.0);
        qeff -> assign(nhit, -2222.0);
        radius -> assign(nhit, -2222.0);
        distanceO2D -> assign(nhit, -2222.0);
        xVD -> assign(nhit, -2222.0);
        yVD -> assign(nhit, -2222.0);
        xP -> assign(nhit, -2222.0);
        yP -> assign(nhit, -2222.0);
        minval -> assign(nhit, -2222.0);
        thetaCh -> assign(nhit, -2222.0);
        xP2 -> assign(nhit, -2222.0);
        yP2 -> assign(nhit, -2222.0);
        minval2 -> assign(nhit, -2222.0);
        thetaCh2 -> assign(nhit, -2222.0);
        if(j%100==0) std::cout << "Vector Size (of " << j << " Event) = " << nhit << std::endl;

        // Tracking
        if(xVD1b->size()>0 && xVD4b->size()>0){

    
            for(Int_t i=0; i<nhit; i++){
            //for(Int_t i=0; i<10; i++){
                // tracking
                /*
                x1 = xVD1b->at(0);//0.0; //xVD1->at(0);
                y1 = yVD1b->at(0);//0.0; //yVD1->at(0);
                z1 = zVD1b->at(0);//-10.0; //zVD1->at(0);
                x2 = xVD4b->at(0);//0.0; //xVD2->at(0);
                y2 = yVD4b->at(0);//0.0; //yVD2->at(0);
                z2 = zVD4b->at(0);//150.0; //zVD2->at(0);

                s1 = (zE-z1)/(z2-z1);

                xE = x1 + s1*(x2-x1);
                yE = y1 + s1*(y2-y1);
                */

                x1 = upstream_trig_centerX;//0.0; //xVD1b->at(0);
                y1 = upstream_trig_centerY;//0.0; //yVD1b->at(0);
                z1 = upstream_trig_centerZ;//-70.0; //zVD1b->at(0);
                
                hitChX = std::distance(hitChXEdepVD4->begin(), std::max_element(hitChXEdepVD4->begin(), hitChXEdepVD4->end()));
                hitChY = std::distance(hitChYEdepVD4->begin(), std::max_element(hitChYEdepVD4->begin(), hitChYEdepVD4->end()));

                x2 = downstream_trig_centerX -0.5*size_beamPM + 0.5*size_1ChbeamPM + hitChX*size_1ChbeamPM; //xVD2b->at(0);
                y2 = downstream_trig_centerY -0.5*size_beamPM + 0.5*size_1ChbeamPM + hitChY*size_1ChbeamPM; //yVD2b->at(0);
                z2 = downstream_trig_centerZ;//130.0; //zVD2b->at(0);

                s1 = (zE-z1)/(z2-z1);

                xE = x1 + s1*(x2-x1);
                yE = y1 + s1*(y2-y1);

                func -> SetParameter(1, xE);
                func -> SetParameter(2, yE);

                qeff->at(i) = qeffVD3->at(i);
                wlen->at(i) = wlenVD3->at(i);
                xD = xVD3->at(i);
                yD = yVD3->at(i);
                zD = zVD3->at(i);
        
                xVD->at(i) = xD;
                yVD->at(i) = yD;

                func -> SetParameter(13, xD);
                func -> SetParameter(14, yD);
                func -> SetParameter(15, zD);           
                
                minimum_value = func->GetMinimumXY(x_minimum, y_minimum);
                minval->at(i) = minimum_value;
                if(minimum_value>minval_threshold){
                    nhit_overThreshold++;
                }
    
    
                xP->at(i) = x_minimum;
                yP->at(i) = y_minimum;
                //minval->at(i) = func->Eval(x_minimum, y_minimum);
                
                radius->at(i) = std::sqrt(std::pow(xD-detec_centerX,2.0) + std::pow(yD-detec_centerY,2.0));

                //std::cout << "Minimum: f(" << xs[0] << "," << xs[1] << "): " << min->MinValue()  << std::endl;
                Vector3D emit(x_minimum-xE, y_minimum-yE, zP - zE);
                Vector3D nbeam(x2-x1, y2-y1, z2-z1);
                thetaCh->at(i) = std::acos(Vector3D::dot(emit, nbeam.normalized())/emit.norm());

                detectionPoint = channel2position.GetChPos3D(parentIdVD3->at(i), hitChVD3->at(i));
                detec2center = Vector3D(xD,yD,zD) - detectionPoint;
                distanceO2D->at(i) = detec2center.norm();


                // Include all error
                
                func -> SetParameter(13, detectionPoint.x);
                func -> SetParameter(14, detectionPoint.y-0.1);
                func -> SetParameter(15, detectionPoint.z);           
                
                minimum_value2 = func->GetMinimumXY(x_minimum2, y_minimum2);
                minval2->at(i) = minimum_value2;
                if(minimum_value2>minval_threshold){
                    nhit_overThreshold2++;
                }                
                xP2->at(i) = x_minimum2;
                yP2->at(i) = y_minimum2;

                Vector3D emit2(x_minimum2-xE, y_minimum2-yE, zP - zE);
                Vector3D nbeam2(x2-x1, y2-y1, z2-z1);
                thetaCh2->at(i) = std::acos(Vector3D::dot(emit2, nbeam2.normalized())/emit2.norm());

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
    
            }
            outTree -> Fill();
        }
        else{
            nsize0++;
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
    std::cout << "Hit Number over threshold (= " << minval_threshold << ") for Pixel: " << nhit_overThreshold2 << std::endl;
    std::cout << "Ratio of over threshold to total for Pixel : " << 100.0*nhit_overThreshold2/nhit_total << " %" << std::endl;
    std::cout << "Beam Particle not detected : " << nsize0 << std::endl;
    std::cout << "Calculation time : " << (double)(end_time-start_time)/CLOCKS_PER_SEC << " sec" << std::endl;

    std::cout << "------------ Analysis Info ------------" << std::endl;
    std::cout << "Finish !!" << std::endl;

    return 0;    
}
