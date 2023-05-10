#include "Ch2Pos.hh"
#include "Vector3D.hh"

Ch2Pos::Ch2Pos(ReadYAML readyaml) {
    auto configuration = readyaml.GetConfigration();
    double detec_sizeX = configuration.detector.size.x;  // cm
    double detec_sizeY = configuration.detector.size.y;  // cm
    double detec_sizeZ = configuration.detector.size.z;  // cm

    double detec_centerX = configuration.detector.center.x;  // cm
    double detec_centerY = configuration.detector.center.y;  // cm
    double detec_centerZ = configuration.detector.center.z;  // cm

    double detec_rotangleX = configuration.detector.rotangle.x;  // deg
    double detec_rotangleY = configuration.detector.rotangle.y;  // deg
    double detec_rotangleZ = configuration.detector.rotangle.z;  // deg

    std::vector<double> PosX_Mppc = readyaml.GetMppcPosX();   // unit cm
    std::vector<double> PosY_Mppc = readyaml.GetMppcPosY();   // unit cm
    std::vector<double> PosZ_Mppc = readyaml.GetMppcPosZ();   // unit cm
    std::vector<double> rot_Z_mppc = readyaml.GetMppcZRot();  // unit deg

    int NumMppc = PosX_Mppc.size();
    //std::cout << "Number of MPPC = " << NumMppc << std::endl;
    std::vector<Vector3D> Pos_MPPC3D;// = std::vector<Vector3D>(NumMppc, Vector3D(1.0,0.0,0.0));
    Pos_MPPC3D.reserve(NumMppc);
    for (int i = 0; i < NumMppc; i++) {
        //Pos_MPPC3D.at(i) = Vector3D(PosX_Mppc.at(i), PosY_Mppc.at(i), PosZ_Mppc.at(i));
        Pos_MPPC3D.push_back(Vector3D(PosX_Mppc.at(i), PosY_Mppc.at(i), PosZ_Mppc.at(i)));
    }

    int NumCh = ChNumOf1Mppc * NumMppc;
    //gChPosX = std::vector<double>(NumCh, -22222.0);
    //gChPosY = std::vector<double>(NumCh, -22222.0);
    //gChPosZ = std::vector<double>(NumCh, -22222.0);

    // Global position of PCB Board
    double gCenterX = detec_centerX + detec_sizeX * std::cos(detec_rotangleY * deg2rad) / 2.0;
    double gCenterY = detec_centerY;
    double gCenterZ = detec_centerZ + detec_sizeX * std::sin(detec_rotangleY * deg2rad) / 2.0;
    Vector3D globalCenter(gCenterX, gCenterY, gCenterZ);
    //std::cout << "globalCenter = " << globalCenter.x << "," << globalCenter.y << ","<< globalCenter.z <<std::endl;
            
    // Relative position of each channels of center of MPPC

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

    // std::vector<double> lChPosX = std::vector<double>(ChNumOf1Mppc, -22222.0);
    // std::vector<double> lChPosY = std::vector<double>(ChNumOf1Mppc, -22222.0);
    // std::vector<double> lChPosZ = std::vector<double>(ChNumOf1Mppc, -22222.0);
    //std::vector<Vector3D> lChPos3D = std::vector<Vector3D>(ChNumOf1Mppc, Vector3D(1.0,0.0,0.0));
    std::vector<Vector3D> lChPos3D;// = std::vector<Vector3D>(ChNumOf1Mppc, Vector3D(1.0,0.0,0.0));
    lChPos3D.reserve(ChNumOf1Mppc);
    for (int i = 0; i < nChRow; i++) {
        for (int j = 0; j < nChCol; j++) {
            //lChPos3D.at(i * nChCol + j) =
            //    Vector3D(-X_BoardSize / 2.0 + (j + 1) * X_GapSize + (0.5 + j) * X_SensorSize,
            //             Y_BoardSize / 2.0 - (i + 1) * Y_GapSize - (0.5 + i) * Y_SensorSize,
            //             Z_BoardSize / 2.0 - Z_WindowSize - Z_SensorSize / 2.0);
            lChPos3D.push_back(
                Vector3D(-X_BoardSize / 2.0 + (j + 1) * X_GapSize + (0.5 + j) * X_SensorSize,
                         Y_BoardSize / 2.0 - (i + 1) * Y_GapSize - (0.5 + i) * Y_SensorSize,
                         Z_BoardSize / 2.0 - Z_WindowSize - Z_SensorSize / 2.0)
            );
            // lChPosX.at(i * nChCol + j) = -X_BoardSize / 2.0 + (j + 1) * X_GapSize + (0.5 + j) *
            // X_SensorSize; lChPosY.at(i * nChCol + j) = Y_BoardSize / 2.0 - (i + 1) * Y_GapSize -
            // (0.5 + i) * Y_SensorSize; lChPosZ.at(i * nChCol + j) = Z_BoardSize / 2.0 -
            // Z_WindowSize - Z_SensorSize / 2.0;
        }
    }
    //std::cout << "lChPos3D = " << lChPos3D.at(10).x << "," << lChPos3D.at(10).y << ","<< lChPos3D.at(10).z <<std::endl;

    // Calculate global postion of each MPPC Channels
    //double gXBoard2MPPC, gYBoard2MPPC, gZBoard2MPPC;
    //double gXMPPC2Ch, gYMPPC2Ch, gZMPPC2Ch;
    Vector3D gBoard2MPPC3D, gMPPC2Ch3D;
    gChPos3D.reserve(NumMppc*ChNumOf1Mppc);// = std::vector<Vector3D>(NumMppc*ChNumOf1Mppc, Vector3D(1.0,0.0,0.0));
    for (int i = 0; i < NumMppc; i++) {
        /*gXBoard2MPPC = -PosZ_Mppc.at(i) * std::sin(detec_rotangleY * deg2rad) +
                       PosX_Mppc.at(i) * std::cos(detec_rotangleY * deg2rad);
        gYBoard2MPPC = PosY_Mppc.at(i);
        gZBoard2MPPC = PosZ_Mppc.at(i) * std::cos(detec_rotangleY * deg2rad) +
                       PosX_Mppc.at(i) * std::sin(detec_rotangleY * deg2rad);
        */
        gBoard2MPPC3D = Pos_MPPC3D.at(i).rotateY(-detec_rotangleY*deg2rad);
        //std::cout << Pos_MPPC3D.at(i).rotateY(-detec_rotangleY*deg2rad) << std::endl;
        //std::cout << "gBoard2MPPC3D = " << gBoard2MPPC3D.x << "," << gBoard2MPPC3D.y << ","<< gBoard2MPPC3D.z <<std::endl;

        for (int j = 0; j < ChNumOf1Mppc; j++) {
            /*gXMPPC2Ch = lChPosX.at(j) * std::cos(rot_Z_mppc.at(i) * deg2rad) -
                        lChPosY.at(j) * std::sin(rot_Z_mppc.at(i) * deg2rad);
            gYMPPC2Ch = -lChPosX.at(j) * std::sin(rot_Z_mppc.at(i) * deg2rad) +
                        lChPosY.at(j) * std::cos(rot_Z_mppc.at(i) * deg2rad);
            gZMPPC2Ch = lChPosZ.at(j);
            */
            gMPPC2Ch3D = lChPos3D.at(j).rotateZ(rot_Z_mppc.at(i)*deg2rad);
            gMPPC2Ch3D = gMPPC2Ch3D.rotateY(-detec_rotangleY*deg2rad);
            //std::cout << "gMPPC2Ch3D = " << gMPPC2Ch3D.x << "," << gMPPC2Ch3D.y << ","<< gMPPC2Ch3D.z <<std::endl;
            //gChPos3D.at(i * ChNumOf1Mppc + j) = globalCenter + gBoard2MPPC3D + gMPPC2Ch3D;
            gChPos3D.push_back(globalCenter + gBoard2MPPC3D + gMPPC2Ch3D);
        }
    }
    //std::cout << "gChPos3D size = " << gChPos3D.size() << std::endl;
    //std::cout << "gChPos3D = " << gChPos3D.at(10).x << "," << gChPos3D.at(10).y << ","<< gChPos3D.at(10).z <<std::endl;
    //std::cout << "gChPos3D = " << gChPos3D.at(100).x << "," << gChPos3D.at(100).y << ","<< gChPos3D.at(100).z <<std::endl;
}

Ch2Pos::~Ch2Pos() {}

Vector3D Ch2Pos::GetChPos3D(int parentId, int Channel){
    return gChPos3D.at((parentId-parentIdOffset)* ChNumOf1Mppc + (Channel-ChannelOffset));
}

double Ch2Pos::GetChPosX(int parentId, int Channel){
    return gChPos3D.at((parentId-parentIdOffset)* ChNumOf1Mppc + (Channel-ChannelOffset)).x;
}

double Ch2Pos::GetChPosY(int parentId, int Channel){
    return gChPos3D.at((parentId-parentIdOffset)* ChNumOf1Mppc + (Channel-ChannelOffset)).y;
}

double Ch2Pos::GetChPosZ(int parentId, int Channel){
    return gChPos3D.at((parentId-parentIdOffset)* ChNumOf1Mppc + (Channel-ChannelOffset)).z;
}

