#ifndef ReadYAML_h
#define ReadYAML_h

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>


class ReadYAML{

private:
    std::string configfilename;
    YAML::Node config;
    std::vector<std::string> keylist;
    std::string settingname;

    int Number_MPPC;
    int parentIdOffset, hitChOffset;
    std::vector<double> mppc_pos_x;// = std::vector<double>(50, -2222.0);
    std::vector<double> mppc_pos_y;// = std::vector<double>(50, -2222.0);
    std::vector<double> mppc_pos_z;// = std::vector<double>(50, -2222.0);
    std::vector<double> mppc_zrot;
    std::vector<double> mppc_darkcurrent;

public:
    ReadYAML(std::string fpath, std::string setting_name);
    std::string GetConfigFile() const {return configfilename;};
    bool GetConfigration();
    int GetNumMPPC() {return Number_MPPC;};
    int GetParentIdOffset() {return parentIdOffset;};
    int GetHitChOffset() {return hitChOffset;};
    std::vector<double> GetMppcPosX() {return mppc_pos_x;};
    std::vector<double> GetMppcPosY() {return mppc_pos_y;};
    std::vector<double> GetMppcPosZ() {return mppc_pos_z;};
    std::vector<double> GetMppcZRot() {return mppc_zrot;};

};



#endif