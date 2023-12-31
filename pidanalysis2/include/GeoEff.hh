#ifndef GeoEff_h
#define GeoEff_h

#include <vector>
#include "ReadYAML.hh"
#include "Vector3D.hh"
#include "Ch2Pos.hh"
#include "ExpConfig.hh"

class ReadYAML;
class ExpConfig;

class GeoEff{
private:
    double deg2rad = std::atan(1.0)*4.0/180.0;

    const int nChRow = 4;  // y
    const int nChCol = 4;  // x
    int ChNumOf1Mppc = nChRow*nChCol;
    int NumMppc;
    int parentIdOffset;// = 100;
    int ChannelOffset;// = 0;
    //std::vector<double> gChPosX;
    //std::vector<double> gChPosY;
    //std::vector<double> gChPosZ;
    std::vector<Vector3D> gChPos3D;

    std::vector<double> PosX_Mppc;   // unit cm
    std::vector<double> PosY_Mppc;   // unit cm
    std::vector<double> PosZ_Mppc;   // unit cm
    std::vector<double> rot_Z_mppc;  // unit deg

    double detec_sizeX;  // cm
    double detec_sizeY;  // cm
    double detec_sizeZ;  // cm

    double detec_centerX;  // cm
    double detec_centerY;  // cm
    double detec_centerZ;  // cm

    double detec_rotangleX;  // deg
    double detec_rotangleY;  // deg
    double detec_rotangleZ;  // deg

    double radiator_centerX;
    double radiator_centerY;
    double radiator_centerZ;
    double radiator_sizeX;
    double radiator_sizeY;
    double radiator_sizeZ;

    double mirror_rotangleX;
    double mirror_rotangleY;
    double mirror_rotangleZ;
    double mirror_radius;

    double nx, ny, nz;
    double a0_Aero;
    double wlen0_Aero;
    double wlength_Aero;
    double n_Air;
    double n_Aero;

public:
    GeoEff(ReadYAML readyaml, ExpConfig expconfig);
    ~GeoEff();

    Vector3D detection_point(Vector3D rE, Vector3D nbeam, double theta, double phi);
    double GetGeoEff(double mass, double momentum, Vector3D rE, Vector3D nbeam);



    double GetChPosX(int parentId, int Channel);
    double GetChPosY(int parentId, int Channel);
    double GetChPosZ(int parentId, int Channel);
    Vector3D GetChPos3D(int parentId, int Channel);


};




#endif