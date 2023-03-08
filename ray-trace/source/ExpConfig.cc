#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "ExpConfig.hh"

ExpConfig::ExpConfig(std::string fpath){
    configfilename = fpath;
    config = YAML::LoadFile(fpath);

    setteing_number = 0;
    auto constructionInfo = config["config"];
    for(auto it=constructionInfo.begin(); it!=constructionInfo.end(); ++it){
        std::string key = it->first.as<std::string>();
        keylist.push_back(key);
        setteing_number++;
    }
}

construction_t ExpConfig::GetConstruction(std::string settingname){
    construction_t construction;

    auto findresult = std::find(keylist.begin(), keylist.end(), settingname);

    if(findresult == keylist.end()){
        std::cout << "Setting Not Found" << std::endl;
    }else{
        construction.setting_name = settingname;
    
        construction.file.input_file = config["config"][settingname]["input_file"].as<std::string>();
        construction.file.output_file = config["config"][settingname]["output_file"].as<std::string>();

        construction.emit.x = config["config"][settingname]["Emit"]["x"].as<double>();
        construction.emit.y = config["config"][settingname]["Emit"]["y"].as<double>();
        construction.emit.z = config["config"][settingname]["Emit"]["z"].as<double>();

        construction.n_perp.x = config["config"][settingname]["n_perp"]["x"].as<double>();
        construction.n_perp.y = config["config"][settingname]["n_perp"]["y"].as<double>();
        construction.n_perp.z = config["config"][settingname]["n_perp"]["z"].as<double>();

        construction.mirror.x0 = config["config"][settingname]["Mirror"]["x0"].as<double>();
        construction.mirror.y0 = config["config"][settingname]["Mirror"]["y0"].as<double>();
        construction.mirror.z0 = config["config"][settingname]["Mirror"]["z0"].as<double>();
        construction.mirror.Radius = config["config"][settingname]["Mirror"]["Radius"].as<double>();

        construction.rindex.Aero = config["config"][settingname]["RIndex"]["Aero"].as<double>();
        construction.rindex.Air = config["config"][settingname]["RIndex"]["Air"].as<double>();
    }
    return construction;
}