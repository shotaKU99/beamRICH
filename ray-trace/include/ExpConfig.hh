#ifndef ExpConfig_h
#define ExpConfig_h

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

struct File_t{
    std::string input_file;
    std::string output_file;
};

struct Emit_t{
    double x;
    double y;
    double z;
};

struct n_perp_t{
    double x;
    double y;
    double z;
};

struct Mirror_t{
    double x0;
    double y0;
    double z0;
    double Radius;
};

struct RIndex_t{
    double Aero;
    double Air;
};

struct construction_t{
    std::string setting_name;
    File_t file;
    Emit_t emit;
    n_perp_t n_perp;
    Mirror_t mirror;
    RIndex_t rindex;
};

class ExpConfig{

private:
    std::string configfilename;
    YAML::Node config;
    std::string info_setting;
    int setteing_number;
    std::vector<std::string> keylist;

public:
    ExpConfig(std::string fpath);

    std::string GetConfigFile() const {return configfilename;};

    const std::vector<std::string>& GetKeyList() const {return keylist;};
    construction_t GetConstruction(std::string settingname);





};



#endif