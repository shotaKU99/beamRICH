#include <time.h>

#include <algorithm>
#include <iostream>
#include <vector>

#include "Ch2Pos.hh"
#include "ExpConfig.hh"
#include "FunctionDef.hh"
#include "Math/Minimizer.h"
#include "ReadYAML.hh"
#include "Rtypes.h"
#include "TError.h"
#include "TF2.h"
#include "TFile.h"
#include "TRandom2.h"
#include "TTree.h"
#include "Vector3D.hh"
#include "GeoEff.hh"

int main(int argc, char** argv) {
    // Read info.yaml
    std::string configfile = "./info.yml";
    std::string settingname = "geant4_sim1";
    ExpConfig expconfig = ExpConfig(configfile);
    auto expinfo = expconfig.GetConstruction(settingname);

    // Read input root file & yaml file
    std::string inputfilename = "sample.root";
    if (argc == 1) {
        std::cout << "Using default inputfile name :" << inputfilename << std::endl;
    } else if (argc == 2) {
        inputfilename = argv[1];
    } else if (argc > 2) {
        std::cout << "=========================================" << std::endl;
        std::cout << "Too Many arguments !!" << std::endl;
        //std::cout << "Use default inputfile name : " << inputfilename << std::endl;
        std::cout << "=========================================" << std::endl;
        exit(EXIT_FAILURE);
    }

    ReadYAML readyaml = ReadYAML(inputfilename + ".yaml");
    auto configuration = readyaml.GetConfigration();

    TFile* inputFile = TFile::Open(inputfilename.c_str());
    TTree* inputTree = (TTree*)inputFile->Get("tree");

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
    int parentIDoffset = readyaml.GetParentIdOffset();
    int hitchoffset = readyaml.GetHitChOffset();
    
    std::cout << "test1" << std::endl;
    GeoEff geomeff = GeoEff(readyaml, expconfig);
    std::cout << "test2" << std::endl;

    Ch2Pos channel2position = Ch2Pos(readyaml);


    // Create output root file & TTree
    const char* outFile_name = expinfo.file.output_file.c_str();
    TFile* outFile = new TFile(outFile_name, "recreate");
    TTree* outTree = new TTree("tree", "tree");

    // definition of angle calculation function
    double xmin = -5.0;
    double xmax = 5.0;
    double ymin = -5.0;
    double ymax = 5.0;

    TF2* func = new TF2("func", distance_func, xmin, xmax, ymin, ymax, 16, 2);

    double zP = radiator_centerZ + radiator_sizeZ / 2.0;
    double xE = 0.0;                                            // configuration.emit.x;
    double yE = 0.0;                                            // configuration.emit.y;
    double zE = radiator_centerZ + radiator_sizeZ / 2.0 - 1.5;  // expinfo.emit.z;
    double nx = expinfo.n_perp.x;
    double ny = expinfo.n_perp.y;
    double nz = expinfo.n_perp.z;
    double x0Mirror =
        mirror_radius * std::sin(mirror_rotangleY * TMath::Pi() / 180.0);  // expinfo.mirror.x0;
    double y0Mirror = 0.0;                                                 // expinfo.mirror.y0;
    double z0Mirror =
        mirror_radius *
        (1.0 - std::cos(mirror_rotangleY * TMath::Pi() / 180.0));  // expinfo.mirror.z0;
    double Radius = mirror_radius;                                 // expinfo.mirror.Radius;
    double a0_Aero = expinfo.rindex.a0_Aero;
    double wlen0_Aero = expinfo.rindex.wlen0_Aero;
    double wlength_Aero = expinfo.rindex.wlength_Aero;
    double n_Air = expinfo.rindex.Air;

    TF1* refractive_index_func = new TF1("refractive_index_func", refractive_index, 200, 900, 2);
    refractive_index_func->SetParameters(a0_Aero, wlen0_Aero);
    double n_Aero = refractive_index_func->Eval(wlength_Aero);

    double xD = 0.0;
    double yD = 0.0;
    double zD = 0.0;

    func->SetParameter(0, zP);
    func->SetParameter(1, xE);
    func->SetParameter(2, yE);
    func->SetParameter(3, zE);
    func->SetParameter(4, nx);
    func->SetParameter(5, ny);
    func->SetParameter(6, nz);
    func->SetParameter(7, x0Mirror);
    func->SetParameter(8, y0Mirror);
    func->SetParameter(9, z0Mirror);
    func->SetParameter(10, Radius);
    func->SetParameter(11, n_Aero);
    func->SetParameter(12, n_Air);
    func->SetParameter(13, xD);
    func->SetParameter(14, yD);
    func->SetParameter(15, zD);

    // definition of tree branch of input root file
    std::vector<Double_t>* xVD1b = 0;
    std::vector<Double_t>* yVD1b = 0;
    std::vector<Double_t>* zVD1b = 0;
    std::vector<Double_t>* xVD4b = 0;
    std::vector<Double_t>* yVD4b = 0;
    std::vector<Double_t>* zVD4b = 0;
    std::vector<Int_t>* nhitchVD3 = 0;
    std::vector<Int_t>* nhitchAllVD3 = 0;
    std::vector<Double_t>* wlenVD3 = 0;
    std::vector<Double_t>* xVD3 = 0;
    std::vector<Double_t>* yVD3 = 0;
    std::vector<Double_t>* zVD3 = 0;
    std::vector<Double_t>* qeffVD3 = 0;
    std::vector<Double_t>* hitChVD3 = 0;
    std::vector<Double_t>* parentIdVD3 = 0;
    std::vector<Double_t>* hitChXEdepVD4 = 0;
    std::vector<Double_t>* hitChYEdepVD4 = 0;

    inputTree->SetBranchAddress("xVD1b", &xVD1b);
    inputTree->SetBranchAddress("yVD1b", &yVD1b);
    inputTree->SetBranchAddress("zVD1b", &zVD1b);
    inputTree->SetBranchAddress("xVD4b", &xVD4b);
    inputTree->SetBranchAddress("yVD4b", &yVD4b);
    inputTree->SetBranchAddress("zVD4b", &zVD4b);
    inputTree->SetBranchAddress("nhitchVD3", &nhitchVD3);
    inputTree->SetBranchAddress("nhitchAllVD3", &nhitchAllVD3);
    inputTree->SetBranchAddress("wlenVD3", &wlenVD3);
    inputTree->SetBranchAddress("xVD3", &xVD3);
    inputTree->SetBranchAddress("yVD3", &yVD3);
    inputTree->SetBranchAddress("zVD3", &zVD3);
    inputTree->SetBranchAddress("qeffVD3", &qeffVD3);
    inputTree->SetBranchAddress("hitChVD3", &hitChVD3);
    inputTree->SetBranchAddress("parentIdVD3", &parentIdVD3);
    inputTree->SetBranchAddress("hitChXEdepVD4", &hitChXEdepVD4);
    inputTree->SetBranchAddress("hitChYEdepVD4", &hitChYEdepVD4);

    // definition of output tree branch
    double LL_pi, LL_K, LL_pro, LL_ele;  // log-likelihood value

    outTree->Branch("LL_pi", &LL_pi);
    outTree->Branch("LL_K", &LL_K);
    outTree->Branch("LL_pro", &LL_pro);
    outTree->Branch("LL_ele", &LL_ele);

    // Variables for Tracker 1,2
    double size_beamPM = 4.8;     // cm
    double size_1ChbeamPM = 0.3;  // cm
    int hitChX, hitChY;
    double x1, y1, z1;
    double x2, y2, z2;
    double x_minimum = 0.0;
    double y_minimum = 0.0;
    double minimum_value = 0.0;
    double s1;
    int nsize0 = 0;
    Vector3D detectionPoint;
    int nChNumber;
    double thetaCh, thetaCh_theory;
    double momentum;
    int mass_hypo;
    double sigma_theta;
    double photon_num_all, photon_num_all_temp;
    double a_ij;

    Int_t mass_hypo_num = 4;  // 0: pi, 1: K, 2: proton, 3: electron
    double mass_list[4] = {139.57018, 493.677, 938.272081, 0.51099895};  // MeV/c^2
    int NumCh1MPPC = 16;
    double area_mppc = 0.3 * 0.3;  // cm^2
    double dark_current_rate = 0.5e6; // Hz
    double tdc_gatewidth = 10e-9; // sec

    // constant
    // for refractive index from lhcb
    //double a0 = 0.039562;
    //double wlen0 = 105.312;
    // for transmittance
    double A = 0.94;
    //double C = 0.0064;
    double d = 0.03;
    double reflectivity_mirror = 0.89;
    double geom_eff;// = 0.25;

    //TF1* refractive_index_func = new TF1("refractive_index_func", refractive_index, 200, 920, 2);
    //refractive_index_func->SetParameters(a0, wlen0);
    //TF1* transmittance_func = new TF1("transmittance_func", transmittance, 200, 920, 3);
    //transmittance_func->SetParameters(A, C, 10.0);
    TF1* transmittance_func2 = new TF1("transmittance_func",transmittance2, 200, 920,5);
    transmittance_func2 -> SetParameters(A, d, 10.0, a0_Aero, wlen0_Aero);

    TF1* qefficiency_func = new TF1("qefficiency_func", quantumEff, 200, 920, 0);

    clock_t start_time = clock();

    // calculate log-likelihood
    for (Int_t j = 0; j < inputTree->GetEntries(); j++) {
    //for (Int_t j = 0; j < 3; j++) {
        inputTree->GetEntry(j);
        nChNumber = nhitchAllVD3->size();  // 24 * 16
        if (j % 100 == 0) {
            std::cout << "Calculating " << j << " event..." << std::endl;
        }
        momentum = BeamMomentum;  // GeV/c

        LL_pi = 0.0;
        LL_K = 0.0;
        LL_pro = 0.0;
        LL_ele = 0.0;
        for (Int_t k = 0; k < mass_hypo_num; k++) {
            // for(Int_t k=0; k<1; k++){
            // std::cout << "vector length (k = " << k << " ) = " << nChNumber << std::endl;
            // thetaCh_theory = 0.199;
            thetaCh_theory =
                std::acos(std::sqrt(std::pow(momentum, 2.0) + std::pow(mass_list[k] * 1e-3, 2.0)) /
                          (momentum * refractive_index_func->Eval(wlength_Aero)));
            sigma_theta = 0.0035; //rad

            if (xVD1b->size() > 0 && xVD4b->size() > 0) {
                for (Int_t i = 0; i < nChNumber; i++) {
                    a_ij = 0.0;

                    if (nhitchAllVD3->at(i) != 0) {
                    //if (nhitchVD3->at(i) != 0) {
                        // std::cout << "test " << std::endl;
                        detectionPoint = channel2position.GetChPos3D(
                            i / NumCh1MPPC + parentIDoffset, i % NumCh1MPPC + hitchoffset);

                        // std::cout << "Detecpos =" << detectionPoint << std::endl;
                        func->SetParameter(13, detectionPoint.x);
                        func->SetParameter(14, detectionPoint.y);
                        func->SetParameter(15, detectionPoint.z);

                        // ここに tracking のループ
                        photon_num_all = 0.0;
                        // for(Int_t i=0; i<10; i++){
                        //  tracking
                        //x1 = xVD1b->at(0);  // 0.0; //xVD1->at(0);
                        //y1 = yVD1b->at(0);  // 0.0; //yVD1->at(0);
                        //z1 = zVD1b->at(0);  //-10.0; //zVD1->at(0);
                        //x2 = xVD4b->at(0);  // 0.0; //xVD2->at(0);
                        //y2 = yVD4b->at(0);  // 0.0; //yVD2->at(0);
                        //z2 = zVD4b->at(0);  // 150.0; //zVD2->at(0);

                        // Include all error
                        x1 = upstream_trig_centerX;//0.0; //xVD1b->at(0);
                        y1 = upstream_trig_centerY;//0.0; //yVD1b->at(0);
                        z1 = upstream_trig_centerZ;//30.0; //zVD1b->at(0);

                        hitChX =
                            std::distance(hitChXEdepVD4->begin(),
                                        std::max_element(hitChXEdepVD4->begin(), hitChXEdepVD4->end()));
                        hitChY =
                            std::distance(hitChYEdepVD4->begin(),
                                        std::max_element(hitChYEdepVD4->begin(), hitChYEdepVD4->end()));

                        x2 = downstream_trig_centerX - 0.5 * size_beamPM + 0.5 * size_1ChbeamPM +
                            hitChX * size_1ChbeamPM;  // xVD2b->at(0);
                        y2 = downstream_trig_centerY - 0.5 * size_beamPM + 0.5 * size_1ChbeamPM +
                            hitChY * size_1ChbeamPM;  // yVD2b->at(0);
                        z2 = downstream_trig_centerZ;  // 130.0; //zVD2b->at(0);


                        s1 = (zE - z1) / (z2 - z1);

                        xE = x1 + s1 * (x2 - x1);
                        yE = y1 + s1 * (y2 - y1);

                        func->SetParameter(1, xE);
                        func->SetParameter(2, yE);

                        minimum_value = func->GetMinimumXY(x_minimum, y_minimum);

                        Vector3D emit2(x_minimum - xE, y_minimum - yE, zP - zE);
                        Vector3D nbeam2(x2 - x1, y2 - y1, z2 - z1);
                        thetaCh =
                            std::acos(Vector3D::dot(emit2, nbeam2.normalized()) / emit2.norm());
                        // std::cout << "thetaCh = " << thetaCh << std::endl;

                        geom_eff = 0.5*geomeff.GetGeoEff(mass_list[k], momentum, Vector3D(xE, yE, zE), nbeam2);
                        //std::cout << "Geometrical efficiency (mass= " << mass_list[k] << ", j = " << j << " ) = " << geom_eff << std::endl;
                        // mu_j(h_j)
                        photon_num_all_temp = allDetectedPhoton(
                            mass_list[k], momentum, radiator_sizeZ, qefficiency_func,
                            transmittance_func2, refractive_index_func, geom_eff);

                        photon_num_all += photon_num_all_temp;
                        // a_ij(h_j)
                        if(geom_eff>0.001 && std::isfinite(thetaCh_theory)==1){
                            a_ij += photon_num_all_temp * 4.0 * area_mppc /
                                    (mirror_radius * mirror_radius * thetaCh) / geom_eff *
                                    angle_probDis(thetaCh, thetaCh_theory, sigma_theta);
                        }else{
                            a_ij += 0.0;
                        }
                        // トラッキングのループ終わり

                        // std::cout << "a_ij (k=" << k << ") = " << a_ij << std::endl;
                        if (k == 0) {
                            LL_pi += nhitchAllVD3->at(i) * std::log(a_ij + dark_current_rate*tdc_gatewidth);
                            //LL_pi += nhitchVD3->at(i) * std::log(a_ij);
                        } else if (k == 1) {
                            LL_K += nhitchAllVD3->at(i) * std::log(a_ij + dark_current_rate*tdc_gatewidth);
                            //LL_K += nhitchVD3->at(i) * std::log(a_ij);
                        } else if (k == 2) {
                            LL_pro += nhitchAllVD3->at(i) * std::log(a_ij + dark_current_rate*tdc_gatewidth);
                            //LL_pro += nhitchVD3->at(i) * std::log(a_ij);
                        } else if (k == 3) {
                            LL_ele += nhitchAllVD3->at(i) * std::log(a_ij + dark_current_rate*tdc_gatewidth);
                            //LL_ele += nhitchVD3->at(i) * std::log(a_ij);
                        }
                    } // if hit is not zero
                } // i MPPC Channel loop (384)
            } // if tracking exists
            // std::cout << "photon number (k= " << k << ") = " << photon_num_all << std::endl;
            if (k == 0) {
                LL_pi -= (photon_num_all + dark_current_rate*tdc_gatewidth*nChNumber);
            } else if (k == 1) {
                LL_K -= (photon_num_all + dark_current_rate*tdc_gatewidth*nChNumber);
            } else if (k == 2) {
                LL_pro -= (photon_num_all + dark_current_rate*tdc_gatewidth*nChNumber);
            } else if (k == 3) {
                LL_ele -= (photon_num_all + dark_current_rate*tdc_gatewidth*nChNumber);
            }
        } // k mass hypo loop
        outTree->Fill();
    } // j event loop
    clock_t end_time = clock();

    outTree->Write();
    outFile->Close();

    std::cout << "------------ Analysis Info ------------" << std::endl;

    std::cout << "Input Root File name : " << inputfilename << std::endl;
    // std::cout << "Config File name : " << configfile << std::endl;
    // std::cout << "Setting name : " << settingname << std::endl;
    std::cout << "Output Root File name :" << outFile_name << std::endl;

    std::cout << "Total Event Size : " << inputTree->GetEntries() << std::endl;
    // std::cout << "Total Hit Number : " << nhit_total << std::endl;
    // std::cout << "Hit Number over threshold (= " << minval_threshold << ") : " <<
    // nhit_overThreshold << std::endl; std::cout << "Ratio of over threshold to total : " <<
    // 100.0*nhit_overThreshold/nhit_total << " %" << std::endl;
    std::cout << "Calculation time : " << (double)(end_time - start_time) / CLOCKS_PER_SEC << " sec"
              << std::endl;

    std::cout << "------------ Analysis Info ------------" << std::endl;
    std::cout << "Finish !!" << std::endl;
    return 0;
}