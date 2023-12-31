#include <iostream>
#include <string>
#include <vector>
#include <algorithm>

#include "ReadYAML.hh"

ReadYAML::ReadYAML(std::string fpath){
    configfilename = fpath;
    config = YAML::LoadFile(fpath);

    int setting_number = 0;
    //auto constructionInfo = config["config"];
    for(auto it=config.begin(); it!=config.end(); ++it){
        std::string key = it->first.as<std::string>();
        keylist.push_back(key);
        setting_number++;
    }

    if(setting_number!=1){
        std::cout << "Setting name is not one" << std::endl;
        exit(EXIT_FAILURE);
    }else{
        settingname = keylist.at(0);
    }
}

configuration_t ReadYAML::GetConfigration(){
    configuration_t configuration;

    configuration.setting_name = settingname;
    
    //configuration.setting.output_file = config[settingname]["output_file"].as<G4String>();
    
    //std::cout << "test 2" << std::endl;
    
    configuration.setting.output_file = config[settingname]["output_file"].as<std::string>();
    configuration.setting.numberOfThread = config[settingname]["numberOfThread"].as<int>();
    //std::cout << "test 3" << std::endl;
    
    configuration.beam.particle = config[settingname]["Beam"]["particle"].as<int>();
    //std::cout << "test 3.5" << std::endl;
    
    configuration.beam.MomSpreadType = config[settingname]["Beam"]["MomSpreadType"].as<int>();
    configuration.beam.momentum.center = config[settingname]["Beam"]["momentum"]["center"].as<double>();
    configuration.beam.momentum.width = config[settingname]["Beam"]["momentum"]["width"].as<double>();
    configuration.beam.xprime.center = config[settingname]["Beam"]["xprime"]["center"].as<double>();
    configuration.beam.xprime.width = config[settingname]["Beam"]["xprime"]["width"].as<double>();
    configuration.beam.yprime.center = config[settingname]["Beam"]["yprime"]["center"].as<double>();
    configuration.beam.yprime.width = config[settingname]["Beam"]["yprime"]["width"].as<double>();
    //std::cout << "test 4" << std::endl;
    
    configuration.beam.PosSpreadType = config[settingname]["Beam"]["PosSpreadType"].as<int>();
    configuration.beam.x.center = config[settingname]["Beam"]["x"]["center"].as<double>();
    configuration.beam.x.width = config[settingname]["Beam"]["x"]["width"].as<double>();
    configuration.beam.y.center = config[settingname]["Beam"]["y"]["center"].as<double>();
    configuration.beam.y.width = config[settingname]["Beam"]["y"]["width"].as<double>();
    configuration.beam.z.center = config[settingname]["Beam"]["z"]["center"].as<double>();
    configuration.beam.z.width = config[settingname]["Beam"]["z"]["width"].as<double>();

    configuration.detector.size.x = config[settingname]["Detector"]["size"]["x"].as<double>();
    configuration.detector.size.y = config[settingname]["Detector"]["size"]["y"].as<double>();
    configuration.detector.size.z = config[settingname]["Detector"]["size"]["z"].as<double>();

    configuration.detector.center.x = config[settingname]["Detector"]["center"]["x"].as<double>();
    configuration.detector.center.y = config[settingname]["Detector"]["center"]["y"].as<double>();
    configuration.detector.center.z = config[settingname]["Detector"]["center"]["z"].as<double>();

    configuration.detector.rotangle.x = config[settingname]["Detector"]["rotangle"]["x"].as<double>();
    configuration.detector.rotangle.y = config[settingname]["Detector"]["rotangle"]["y"].as<double>();
    configuration.detector.rotangle.z = config[settingname]["Detector"]["rotangle"]["z"].as<double>();

    configuration.mirror.rotangle.x = config[settingname]["Mirror"]["rotangle"]["x"].as<double>();
    configuration.mirror.rotangle.y = config[settingname]["Mirror"]["rotangle"]["y"].as<double>();
    configuration.mirror.rotangle.z = config[settingname]["Mirror"]["rotangle"]["z"].as<double>();
    configuration.mirror.radius = config[settingname]["Mirror"]["radius"].as<double>();

    configuration.radiator.center.x = config[settingname]["Radiator"]["center"]["x"].as<double>();
    configuration.radiator.center.y = config[settingname]["Radiator"]["center"]["y"].as<double>();
    configuration.radiator.center.z = config[settingname]["Radiator"]["center"]["z"].as<double>();
    configuration.radiator.size.x = config[settingname]["Radiator"]["size"]["x"].as<double>();
    configuration.radiator.size.y = config[settingname]["Radiator"]["size"]["y"].as<double>();
    configuration.radiator.size.z = config[settingname]["Radiator"]["size"]["z"].as<double>();

    configuration.uptrig.center.x = config[settingname]["TriggerScinti"]["upstream"]["center"]["x"].as<double>();
    configuration.uptrig.center.y = config[settingname]["TriggerScinti"]["upstream"]["center"]["y"].as<double>();
    configuration.uptrig.center.z = config[settingname]["TriggerScinti"]["upstream"]["center"]["z"].as<double>();

    configuration.downtrig.center.x = config[settingname]["TriggerScinti"]["downstream"]["center"]["x"].as<double>();
    configuration.downtrig.center.y = config[settingname]["TriggerScinti"]["downstream"]["center"]["y"].as<double>();
    configuration.downtrig.center.z = config[settingname]["TriggerScinti"]["downstream"]["center"]["z"].as<double>();
   

    //std::cout << "test readyaml1" << std::endl;
    parentIdOffset = config[settingname]["PhotonDetector"]["parentIdOffset"].as<int>();
    hitChOffset = config[settingname]["PhotonDetector"]["hitChOffset"].as<int>();
    auto configurationInfo_mppc_pos = config[settingname]["PhotonDetector"]["MPPC"];
    Number_MPPC = configurationInfo_mppc_pos.size();
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



    return configuration;
}