#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "ReadYAML.hh"

ReadYAML::ReadYAML(std::string fpath, std::string setting_name){
    configfilename = fpath;
    config = YAML::LoadFile(fpath);
    settingname = setting_name;

    int setting_number = 0;
    //auto constructionInfo = config["config"];
    for(auto it=config.begin(); it!=config.end(); ++it){
        std::string key = it->first.as<std::string>();
        keylist.push_back(key);
        setting_number++;
    }
}

bool ReadYAML::GetConfigration(){
      
    //std::cout << "test readyaml1" << std::endl;
    parentIdOffset = config["config"][settingname]["Detector"]["parentIdOffset"].as<int>();
    hitChOffset = config["config"][settingname]["Detector"]["hitChOffset"].as<int>();
    auto configurationInfo_mppc_pos = config["config"][settingname]["Detector"]["MPPC"];
    Number_MPPC = configurationInfo_mppc_pos.size();
    std::cout << "number of mppc = " << Number_MPPC << std::endl;
    //std::cout << "test readyaml2" << std::endl;

    //std::cout << "Number of MPPC = " << Number_MPPC << std::endl;

    mppc_pos_x = std::vector<double>(Number_MPPC, -2222.0);
    mppc_pos_y = std::vector<double>(Number_MPPC, -2222.0);
    mppc_pos_z = std::vector<double>(Number_MPPC, -2222.0);
    mppc_zrot = std::vector<double>(Number_MPPC, -2222.0);
    mppc_darkcurrent = std::vector<double>(Number_MPPC, -2222.0);
    
    for(int i=0; i<Number_MPPC; i++){
        mppc_pos_x.at(i) = configurationInfo_mppc_pos[std::to_string(i)]["x"].as<double>();
        mppc_pos_y.at(i) = configurationInfo_mppc_pos[std::to_string(i)]["y"].as<double>();
        mppc_pos_z.at(i) = configurationInfo_mppc_pos[std::to_string(i)]["z"].as<double>();
        mppc_zrot.at(i) = configurationInfo_mppc_pos[std::to_string(i)]["zrot"].as<double>();
        mppc_darkcurrent.at(i) = configurationInfo_mppc_pos[std::to_string(i)]["darkcurrent"].as<double>();
    }    



    return true;
}