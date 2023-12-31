#include "GeoEff.hh"

#include <vector>

#include "Ch2Pos.hh"
#include "FunctionDef.hh"
#include "ReadYAML.hh"

GeoEff::GeoEff(ReadYAML readyaml, ExpConfig expconfig) {
    auto configuration = readyaml.GetConfigration();
    auto expinfo = expconfig.GetConstruction("geant4_sim1");

    detec_sizeX = configuration.detector.size.x;  // cm
    detec_sizeY = configuration.detector.size.y;  // cm
    detec_sizeZ = configuration.detector.size.z;  // cm

    detec_centerX = configuration.detector.center.x;  // cm
    detec_centerY = configuration.detector.center.y;  // cm
    detec_centerZ = configuration.detector.center.z;  // cm

    detec_rotangleX = configuration.detector.rotangle.x;  // deg
    detec_rotangleY = configuration.detector.rotangle.y;  // deg
    detec_rotangleZ = configuration.detector.rotangle.z;  // deg

    radiator_centerX = configuration.radiator.center.x;
    radiator_centerY = configuration.radiator.center.y;
    radiator_centerZ = configuration.radiator.center.z;
    radiator_sizeX = configuration.radiator.size.x;
    radiator_sizeY = configuration.radiator.size.y;
    radiator_sizeZ = configuration.radiator.size.z;

    mirror_rotangleX = configuration.mirror.rotangle.x;
    mirror_rotangleY = configuration.mirror.rotangle.y;
    mirror_rotangleZ = configuration.mirror.rotangle.z;
    mirror_radius = configuration.mirror.radius;

    nx = expinfo.n_perp.x;
    ny = expinfo.n_perp.y;
    nz = expinfo.n_perp.z;

    a0_Aero = expinfo.rindex.a0_Aero;
    wlen0_Aero = expinfo.rindex.wlen0_Aero;
    wlength_Aero = expinfo.rindex.wlength_Aero;
    n_Air = expinfo.rindex.Air;

    TF1* refractive_index_func = new TF1("refractive_index_func", refractive_index, 200, 900, 2);
    refractive_index_func->SetParameters(a0_Aero, wlen0_Aero);
    n_Aero = refractive_index_func->Eval(wlength_Aero);

    parentIdOffset = readyaml.GetParentIdOffset();
    ChannelOffset = readyaml.GetHitChOffset();
    PosX_Mppc = readyaml.GetMppcPosX();   // unit cm
    PosY_Mppc = readyaml.GetMppcPosY();   // unit cm
    PosZ_Mppc = readyaml.GetMppcPosZ();   // unit cm
    rot_Z_mppc = readyaml.GetMppcZRot();  // unit deg
    NumMppc = PosX_Mppc.size();

    std::vector<Vector3D> Pos_MPPC3D;
    Pos_MPPC3D.reserve(NumMppc);
    for (int i = 0; i < NumMppc; i++) {
        Pos_MPPC3D.push_back(Vector3D(PosX_Mppc.at(i), PosY_Mppc.at(i), PosZ_Mppc.at(i)));
    }

    // Global position of PCB Board
    double gCenterX = detec_centerX + detec_sizeX * std::cos(detec_rotangleY * deg2rad) / 2.0;
    double gCenterY = detec_centerY;
    double gCenterZ = detec_centerZ + detec_sizeX * std::sin(detec_rotangleY * deg2rad) / 2.0;
    Vector3D globalCenter(gCenterX, gCenterY, gCenterZ);

    // MPPC Sensor
    double X_SensorSize = 0.30;  // unit of cm
    double Y_SensorSize = 0.30;
    double Z_SensorSize = 0.04;
    // MPPC Window
    double X_WindowSize = 1.30;  // unit of cm
    double Y_WindowSize = 1.30;
    double Z_WindowSize = 0.01;
    // Board base
    double X_BoardSize = 1.30;
    double Y_BoardSize = 1.30;
    double Z_BoardSize = 0.135;
    // Gap between channels
    double X_GapSize = (X_BoardSize - nChCol * X_SensorSize) / (nChCol + 1);
    double Y_GapSize = (Y_BoardSize - nChRow * Y_SensorSize) / (nChRow + 1);

    std::cout << "test1_1" << std::endl;

    std::vector<Vector3D> lChPos3D;
    lChPos3D.reserve(ChNumOf1Mppc);
    for (int i = 0; i < nChRow; i++) {
        for (int j = 0; j < nChCol; j++) {
            lChPos3D.push_back(
                Vector3D(-X_BoardSize / 2.0 + (j + 1) * X_GapSize + (0.5 + j) * X_SensorSize,
                         Y_BoardSize / 2.0 - (i + 1) * Y_GapSize - (0.5 + i) * Y_SensorSize,
                         Z_BoardSize / 2.0 - Z_WindowSize - Z_SensorSize / 2.0));
        }
    }
    std::cout << "test1_2" << std::endl;
    //std::cout << "lChPos3D size = " << lChPos3D.size() << std::endl;
    //std::cout << "rot_Z_mppc size = " << rot_Z_mppc.size() << std::endl;

    // Calculate global postion of each MPPC Channels
    Vector3D gBoard2MPPC3D, gMPPC2Ch3D;
    gChPos3D.reserve(GeoEff::NumMppc * GeoEff::ChNumOf1Mppc);
    for (int i = 0; i < GeoEff::NumMppc; i++) {
        //std::cout << "test1_3" << std::endl;
        gBoard2MPPC3D = Pos_MPPC3D.at(i).rotateY(-detec_rotangleY * deg2rad);
        //std::cout << "test1_4" << std::endl;

        for (int j = 0; j < GeoEff::ChNumOf1Mppc; j++) {
            gMPPC2Ch3D = lChPos3D.at(j).rotateZ(rot_Z_mppc.at(i) * deg2rad);
            //std::cout << "test1_5" << std::endl;
            gMPPC2Ch3D = gMPPC2Ch3D.rotateY(-detec_rotangleY * deg2rad);
            //std::cout << "test1_6" << std::endl;
            gChPos3D.push_back(globalCenter + gBoard2MPPC3D + gMPPC2Ch3D);
        }
    }
    std::cout << "test1_7" << std::endl;


}

GeoEff::~GeoEff() {}

Vector3D GeoEff::GetChPos3D(int parentId, int Channel) {
    return gChPos3D.at((parentId - parentIdOffset) * ChNumOf1Mppc + (Channel - ChannelOffset));
}

double GeoEff::GetChPosX(int parentId, int Channel) {
    return gChPos3D.at((parentId - parentIdOffset) * ChNumOf1Mppc + (Channel - ChannelOffset)).x;
}

double GeoEff::GetChPosY(int parentId, int Channel) {
    return gChPos3D.at((parentId - parentIdOffset) * ChNumOf1Mppc + (Channel - ChannelOffset)).y;
}

double GeoEff::GetChPosZ(int parentId, int Channel) {
    return gChPos3D.at((parentId - parentIdOffset) * ChNumOf1Mppc + (Channel - ChannelOffset)).z;
}

Vector3D GeoEff::detection_point(Vector3D rE, Vector3D nbeam, double theta, double phi) {
    double deg2rad = std::atan(1.0) * 4.0 / 180.0;

    // Emission Point 1
    double xE1 = rE.x;  // par[0];
    double yE1 = rE.y;  // par[1];
    double zE1 = rE.z;  // par[2];

    double zE2 = radiator_centerZ + radiator_sizeZ / 2.0;

    // Perpendicular of back of Radiator
    double nx = GeoEff::nx;
    double ny = GeoEff::ny;
    double nz = GeoEff::nz;

    // Center & Radius of Spherical Mirror
    double Radius = mirror_radius;
    double MirrorRotY = mirror_rotangleY;  // deg
    double x0Mirror = Radius * std::sin(MirrorRotY * deg2rad);
    double y0Mirror = 0.0;
    double z0Mirror = Radius * (1.0 - std::cos(MirrorRotY * deg2rad));

    // Refractive Index
    double n_Aero = GeoEff::n_Aero;
    double n_Air = GeoEff::n_Air;

    // Detection Point
    // double xD = par[12];
    // double yD = par[13];
    // double zD = par[14];

    double xD0 = GeoEff::detec_centerX;
    double yD0 = GeoEff::detec_centerY;
    double zD0 = GeoEff::detec_centerZ;

    double rotDY = GeoEff::detec_rotangleY;

    // Vector3D emit_dir(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi),
    //                   std::cos(theta));

    Vector3D nbeam_uni = nbeam.normalized();

    double theta_beam = std::acos(Vector3D::dot(nbeam_uni, Vector3D(0, 0, 1)));

    // 返り値 0-pi でも 0-2pi にしたい
    double phi_beam = std::acos(nbeam_uni.x / std::sin(theta_beam));
    if(nbeam_uni.y<0){
        phi_beam = 2.0*M_PI - phi_beam;
    }

    // (theta_b+theta,phi_b)方向のベクトル
    Vector3D emit_base(std::sin(theta_beam + theta) * std::cos(phi_beam),
                        std::sin(theta_beam + theta) * std::sin(phi_beam),
                        std::cos(theta_beam + theta));

    // ロドリゲスの回転公式でチェレンコフ光の角度を計算
    // -phi_beam + phi 回転
    double rotAngle_phi = -phi_beam + phi;
    Vector3D RotMat1(
        std::pow(nbeam_uni.x, 2.0) * (1 - std::cos(rotAngle_phi)) + std::cos(rotAngle_phi),
        nbeam_uni.x * nbeam_uni.y * (1 - std::cos(rotAngle_phi)) -
            nbeam_uni.z * std::sin(rotAngle_phi),
        nbeam_uni.x * nbeam_uni.z * (1 - std::cos(rotAngle_phi)) +
            nbeam_uni.y * std::sin(rotAngle_phi)
    );
    Vector3D RotMat2(
        nbeam_uni.x * nbeam_uni.y * (1 - std::cos(rotAngle_phi)) +
            nbeam_uni.z * std::sin(rotAngle_phi),
        std::pow(nbeam_uni.y, 2.0) * (1 - std::cos(rotAngle_phi)) + std::cos(rotAngle_phi),
        nbeam_uni.y * nbeam_uni.z * (1 - std::cos(rotAngle_phi)) -
            nbeam_uni.x * std::sin(rotAngle_phi)
    );
    Vector3D RotMat3(
        nbeam_uni.x * nbeam_uni.z * (1 - std::cos(rotAngle_phi)) -
            nbeam_uni.y * std::sin(rotAngle_phi),
        nbeam_uni.y * nbeam_uni.z * (1 - std::cos(rotAngle_phi)) +
            nbeam_uni.x * std::cos(rotAngle_phi),
        std::pow(nbeam_uni.z, 2.0) * (1 - std::cos(rotAngle_phi)) + std::cos(rotAngle_phi)
    );

    Vector3D emit_dir;
    emit_dir.x = Vector3D::dot(RotMat1, emit_base);
    emit_dir.y = Vector3D::dot(RotMat2, emit_base);
    emit_dir.z = Vector3D::dot(RotMat3, emit_base);

    Vector3D rE1(xE1, yE1, zE1);
    Vector3D n_perp(nx, ny, nz);
    Vector3D r0Mirror(x0Mirror, y0Mirror, z0Mirror);
    Vector3D rD0(xD0, yD0, zD0);

    // Refraction
    Vector3D d_refraction = Vector3D::refraction(emit_dir, n_perp, n_Aero, n_Air);

    Vector3D rP = rE1 + emit_dir * (zE2 - zE1) / emit_dir.z;

    // Cross-Point of refraction light and Spherical Mirror
    Vector3D rPM0 = rP - r0Mirror;
    double t_c = -1.0 * Vector3D::dot(rPM0, d_refraction) +
                 std::sqrt(Radius * Radius - rPM0.norm2() +
                           std::pow(Vector3D::dot(rPM0, d_refraction), 2.0));
    Vector3D rM = rP + t_c * d_refraction;

    // Reflect on Spherical Mirror
    Vector3D d_reflect = Vector3D::reflect(d_refraction, rM - r0Mirror);

    // calculate intersection of detection plane and reflected light
    Vector3D nperp_detec = Vector3D(0.0, 0.0, 1.0).rotateY(-rotDY * deg2rad);

    return rM + d_reflect * (Vector3D::dot(nperp_detec, rD0) - Vector3D::dot(nperp_detec, rM)) /
                    Vector3D::dot(nperp_detec, d_reflect);
}

double GeoEff::GetGeoEff(double mass, double momentum, Vector3D rE, Vector3D nbeam) {
    // 0: e, 1:pi, 2: K, 3: proton
    //double mass_list[4] = {0.51099895, 139.57018, 493.677, 938.272081};  // MeV/c^2

    double theta =
        std::acos(std::sqrt(std::pow(momentum, 2.0) + std::pow(mass * 1e-3, 2.0)) /
                  (momentum * GeoEff::n_Aero));

    if(std::isfinite(theta)==0){
        return 0.0;
    }

    // MPPC Sensor
    double X_SensorSize = 0.30;  // unit of cm
    double Y_SensorSize = 0.30;
    double Z_SensorSize = 0.04;
    // MPPC Window
    double X_WindowSize = 1.30;  // unit of cm
    double Y_WindowSize = 1.30;
    double Z_WindowSize = 0.01;
    // Board base
    double X_BoardSize = 1.30;
    double Y_BoardSize = 1.30;
    double Z_BoardSize = 0.135;
    // Gap between channels
    double X_GapSize = (X_BoardSize - nChCol * X_SensorSize) / (nChCol + 1);
    double Y_GapSize = (Y_BoardSize - nChRow * Y_SensorSize) / (nChRow + 1);
    

    double phi;
    Vector3D rD1;
    double ypos_mppc, ypos_ch, xpos_ch;
    int hitch_base;
    
    int yhit_no1;
    int yhit_no2;
    int xhit_no;


    int nstep = 200;
    int nhit = 0;
    for (int i = 0; i < nstep; i++) {
        // for(int i=0; i<2; i++){
        phi = -0.5 * M_PI + M_PI * i / nstep; // -pi/2 - pi/2

        rD1 = detection_point(rE, nbeam, theta, phi);
        //if(i%10==0) std::cout << "detection point (phi = " << phi/M_PI <<  ") = " << rD1 << std::endl;
        

        // Check y position
        yhit_no1 = -1;
        for(int j=0; j<NumMppc/2; j++){
            ypos_mppc = GeoEff::PosY_Mppc.at(2*j);

            if(std::abs(ypos_mppc-rD1.y)<Y_BoardSize/2.0){
                yhit_no1 = 2*j;
                //break;
            }else{
                continue;
            }
        }
        //if(i%20==0) std::cout << "yhit_no1 = " << yhit_no1 << std::endl;
        // if hit position is not hit to mppc
        if(yhit_no1==-1){
            continue;
        }

        yhit_no2 == -1;
        hitch_base = yhit_no1*ChNumOf1Mppc;
        for(int j=0; j<nChRow; j++){ // y hit
            ypos_ch = GeoEff::GetChPosY(yhit_no1+parentIdOffset, 4*j+ChannelOffset);
            if(std::abs(ypos_ch-rD1.y)<Y_SensorSize/2.0){
                yhit_no2 = 4*j;
                //break;
            }else{
                continue;
            }
        }
        //if(i%20==0) std::cout << "yhit_no2 = " << yhit_no2 << std::endl;
        if(yhit_no2==-1){
            continue;
        }

        xhit_no = -1;
        for(int j=0; j<2*nChCol; j++){
            if(j<4){
                xpos_ch = GeoEff::GetChPosX(yhit_no1+parentIdOffset, j+ChannelOffset);
            }else{
                xpos_ch = GeoEff::GetChPosX(yhit_no1+1+parentIdOffset, j-nChCol+ChannelOffset);
            }

            if(std::abs(xpos_ch-rD1.x)<X_SensorSize*std::cos(detec_rotangleY*deg2rad)/2.0){
                xhit_no = j;
                //break;
            }else{
                continue;
            }
        }
        //if(i%10==0) std::cout << "xhit_no = " << xhit_no << std::endl;
        if(xhit_no==-1){
            continue;
        }

        // xhit_no != -1 && yhit_no2 != -1
        nhit++;    
    }

    // return geometorical efficiency
    return (double)(nhit)/nstep;
}