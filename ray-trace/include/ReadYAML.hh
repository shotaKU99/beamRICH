#ifndef ReadYAML_h
#define ReadYAML_h

#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

struct Setting_t {
    //G4String output_file;
    std::string output_file;
    int numberOfThread;
};

struct Momentum_t {
    double center;
    double width;
};

struct Xprime_t {
    double center;
    double width;
};

struct Yprime_t {
    double center;
    double width;
};

struct X_t {
    double center;
    double width;
};

struct Y_t {
    double center;
    double width;
};

struct Z_t {
    double center;
    double width;
};

struct Beam_t {
    int particle;
    int MomSpreadType;
    Momentum_t momentum;
    Xprime_t xprime;
    Yprime_t yprime;
    int PosSpreadType;
    X_t x;
    Y_t y;
    Z_t z;
};

struct Center_t
{
  double x;
  double y;
  double z;
};

struct RotAngle_t
{
  double x;
  double y;
  double z;
};


struct Detector_t
{
  Center_t center;
  RotAngle_t rotangle;
};

struct Mirror_t
{
  RotAngle_t rotangle;
  double radius;
};

struct BoxSize_t
{
  double x;
  double y;
  double z;
};

struct Radiator_t
{
  Center_t center;
  BoxSize_t size;
};

struct configuration_t {
    //G4String setting_name;
    std::string setting_name;
    Setting_t setting;
    Beam_t beam;
    Detector_t detector;
    Mirror_t mirror;
    Radiator_t radiator;
};


class ReadYAML{

private:
    std::string configfilename;
    YAML::Node config;
    std::vector<std::string> keylist;
    std::string settingname;

public:
    ReadYAML(std::string fpath);

    std::string GetConfigFile() const {return configfilename;};

    configuration_t GetConfigration();
};



#endif