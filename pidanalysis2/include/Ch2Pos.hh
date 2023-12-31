#ifndef Ch2Pos_h
#define Ch2Pos_h

#include <vector>
#include "ReadYAML.hh"
#include "Vector3D.hh"


class ReadYAML;

class Ch2Pos
{
private:
    double deg2rad = std::atan(1.0)*4.0/180.0;

    const int nChRow = 4;  // y
    const int nChCol = 4;  // x
    int ChNumOf1Mppc = nChRow*nChCol;
    int parentIdOffset;// = 100;
    int ChannelOffset;// = 0;
    //std::vector<double> gChPosX;
    //std::vector<double> gChPosY;
    //std::vector<double> gChPosZ;
    std::vector<Vector3D> gChPos3D;


public:
    Ch2Pos(ReadYAML readyaml);
    ~Ch2Pos();

    double GetChPosX(int parentId, int Channel);
    double GetChPosY(int parentId, int Channel);
    double GetChPosZ(int parentId, int Channel);
    Vector3D GetChPos3D(int parentId, int Channel);
};












#endif